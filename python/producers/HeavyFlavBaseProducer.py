import os
import itertools
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.utils import deltaR, closest, polarP4, sumP4, get_subjets, configLogger
from ..helpers.nnHelper import convert_prob
from ..helpers.jetmetCorrector import JetMETCorrector, rndSeed, JetIdCalculator, AK15patcher
from .leptonSFProducer import MuonSFProducer
from .puWeightProducer import PileupWeightProducer
from .EventVetoMapProducer import EventVetoMapProducer

import logging
logger = logging.getLogger('nano')
configLogger('nano', loglevel=logging.INFO)

lumi_dict = {2015: 19.52, 2016: 16.81, 2017: 41.48, 2018: 59.83, 2021: 7.98, 2022: 26.67, 2023:17.794, 2024:109.95}

class _NullObject:
    '''An null object which does not store anything, and does not raise exception.'''

    def __bool__(self):
        return False

    def __nonzero__(self):
        return False

    def __getattr__(self, name):
        pass

    def __setattr__(self, name, value):
        pass


class METObject(Object):

    def p4(self):
        return polarP4(self, eta=None, mass=None)


class HeavyFlavBaseProducer(Module, object):

    def __init__(self, channel, **kwargs):
        self._channel = channel  # 'muon'
        self.year = int(kwargs['year'])
        self.jetType = kwargs.get('jetType', 'ak15').lower()
        self._jmeSysts = {'jec': False, 'jes': None, 'jes_source': '', 'jes_uncertainty_file_prefix': '',
                          'jer': 'nominal', 'met_unclustered': None, 'applyHEMUnc': False,
                          'smearMET': False}
        self._opts = {'fillSystWeights': True, 'nanoV15': False}
        for k in kwargs:
            if k in self._jmeSysts:
                self._jmeSysts[k] = kwargs[k]
            else:
                self._opts[k] = kwargs[k]
        self._needsJMECorr = any([self._jmeSysts['jec'], self._jmeSysts['jes'],
                                  self._jmeSysts['jer'],
                                  self._jmeSysts['met_unclustered'], self._jmeSysts['applyHEMUnc']])
        self._doJetCleaning = True

        logger.info('Running %s channel for %s jets with JME systematics %s, other options %s',
                    self._channel, self.jetType, str(self._jmeSysts), str(self._opts))

        if self.jetType != 'ak15':
            raise RuntimeError('Jet type %s is not supported (only ak15)!' % self.jetType)
        self._jetConeSize = 1.5
        self._fj_name = 'AK15Puppi'
        self._fj_gen_name = 'GenJetAK15'

        self.jetidcalc = JetIdCalculator(year=self.year, jetType="AK4PUPPI", nanov15=self._opts['nanoV15'])
        self.fatjetidcalc = JetIdCalculator(year=self.year, jetType="AK8PUPPI", nanov15=self._opts['nanoV15'])
        self.ak15patcher = AK15patcher(year=self.year)

        if self._needsJMECorr:
            self.jetmetCorr = JetMETCorrector(year=self.year, jetType="AK4PFPuppi", **self._jmeSysts)
            self.fatjetCorr = JetMETCorrector(year=self.year, jetType="AK8PFPuppi", **self._jmeSysts)
            self.subjetCorr = JetMETCorrector(year=self.year, jetType="AK4PFPuppi", **self._jmeSysts)

        # AK4 b-tagging, medium WP
        # NanoAODv15 (2024) has no PNet b-tag for AK4 jets: use UParTAK4 (BTV Summer24 json, M = 0.1272)
        self.btag_name, self.btag_WP_M = {
            2021: ('btagPNetB', 0.245),
            2022: ('btagPNetB', 0.2605),
            2023: ('btagPNetB', 0.1917),
            2024: ('btagUParTAK4B', 0.1272),
        }[self.year]

        # same modules as in NanoTrees vhTreeProducer, 1L channel; the muon iso SF matches miniIso < 0.10
        # no trigger SF (as in the original HRT muon channel: HLT_Mu50 is only stored as passMuTrig)
        self._modules = {
            'muonSF': MuonSFProducer,
            'puWeight': PileupWeightProducer,
            'eventJetVeto': EventVetoMapProducer,
        }
        for k, cls in self._modules.items():
            self._modules[k] = cls(self.year, fillSystWeights=self._opts['fillSystWeights'])
        self._modules['muonSF'].SetIDWP('Tight')
        self._modules['muonSF'].SetIsoWP('TightMini')

    def beginJob(self):
        self.jetidcalc.beginJob()
        self.fatjetidcalc.beginJob()
        if self._needsJMECorr:
            self.jetmetCorr.beginJob()
            self.fatjetCorr.beginJob()
            self.subjetCorr.beginJob()
        for mod in self._modules.values():
            mod.beginJob()

    def endJob(self):
        for mod in self._modules.values():
            mod.endJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.rho_branch_name = 'Rho_fixedGridRhoFastjetAll' if bool(
            inputTree.GetBranch('Rho_fixedGridRhoFastjetAll')) else 'fixedGridRhoFastjetAll'

        self.out = wrappedOutputTree

        # NOTE: branch names must start with a lower case letter
        # check keep_and_drop_output.txt
        self.out.branch("year", "I")
        self.out.branch("lumiwgt", "F")
        self.out.branch("jetR", "F")
        self.out.branch("passmetfilters", "O")
        self.out.branch("l1PreFiringWeight", "F")
        self.out.branch("l1PreFiringWeightUp", "F")
        self.out.branch("l1PreFiringWeightDown", "F")
        self.out.branch("nlep", "I")
        self.out.branch("ht", "F")
        self.out.branch("met", "F")
        self.out.branch("metphi", "F")

        # Large-R jets
        for idx in [1]:
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")
            self.out.branch(prefix + "nsubjets", "I")

            # masses, same definitions as in NanoTrees vhTreeProducer
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "regressed_mass", "F")
            self.out.branch(prefix + "ParT_resonanceMass", "F")
            self.out.branch(prefix + "ParT_genericMass", "F")
            self.out.branch(prefix + "ParT_X2pMass", "F")

            self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
            self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
            self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
            self.out.branch(prefix + "ParticleNetMD_QCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccOrXqqVsQCD", "F")

            self.out.branch(prefix + "ParT_HbbVsQCD", "F")
            self.out.branch(prefix + "ParT_HccVsQCD", "F")
            self.out.branch(prefix + "ParT_TVsQCD", "F")
            self.out.branch(prefix + "ParT_QCD", "F")
            for score in self._part_scores:
                self.out.branch(prefix + "ParT_" + score, "F")

            # matching variables
            if self.isMC:
                self.out.branch(prefix + "nbhadrons", "I")
                self.out.branch(prefix + "nchadrons", "I")
                # info of the closest hadGenH
                self.out.branch(prefix + "dr_H", "F")
                self.out.branch(prefix + "dr_H_daus", "F")
                self.out.branch(prefix + "H_pt", "F")
                self.out.branch(prefix + "H_decay", "I")

                # info of the closest hadGenZ
                self.out.branch(prefix + "dr_Z", "F")
                self.out.branch(prefix + "dr_Z_daus", "F")
                self.out.branch(prefix + "Z_pt", "F")
                self.out.branch(prefix + "Z_decay", "I")

                # info of the closest hadGenW
                self.out.branch(prefix + "dr_W", "F")
                self.out.branch(prefix + "dr_W_daus", "F")
                self.out.branch(prefix + "W_pt", "F")
                self.out.branch(prefix + "W_decay", "I")

                # info of the closest hadGenTop
                self.out.branch(prefix + "dr_T", "F")
                self.out.branch(prefix + "dr_T_b", "F")
                self.out.branch(prefix + "dr_T_Wq_max", "F")
                self.out.branch(prefix + "dr_T_Wq_min", "F")
                self.out.branch(prefix + "T_Wq_max_pdgId", "I")
                self.out.branch(prefix + "T_Wq_min_pdgId", "I")
                self.out.branch(prefix + "T_pt", "F")

        for mod in self._modules.values():
            mod.beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        for mod in self._modules.values():
            mod.endFile(inputFile, outputFile, inputTree, wrappedOutputTree)

    # ParTv3 raw scores stored in the output (NanoAODv15 names)
    _part_scores = ['Xbb', 'Xcc', 'Xcs', 'Xqq', 'TopbWqq', 'TopbWq']

    def selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []  # used for jet lepton cleaning & lepton counting

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            # follow the main analysis:
            if el.pt > 20 and abs(el.eta) < 2.5 and el.mvaIso_WP90:
                event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            # follow the main analysis, looseId :
            if mu.pt > 20 and abs(mu.eta) < 2.4 and mu.looseId and mu.pfRelIso04_all < 0.25:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

    def correctJetsAndMET(self, event):
        # same as NanoTrees vhTreeProducer: analyze() + _correctJetAndMET()
        event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)
        event._allAK15jets = Collection(event, "AK15Puppi")
        event._allAK8jets = Collection(event, "FatJet")
        event._allJets = Collection(event, "Jet")
        event.ak15Subjets = Collection(event, "AK15PuppiSubJet")  # do not sort after updating!!
        event.met = METObject(event, "PuppiMET")

        self.ak15patcher.CalculateMissingVars(event._allJets, event._allAK8jets, event._allAK15jets)
        self.jetidcalc.CalculateJetId(event._allJets)
        self.fatjetidcalc.CalculateJetId(event._allAK15jets)

        if self._needsJMECorr:
            rho = getattr(event, self.rho_branch_name)
            # correct AK4 jets and MET
            self.jetmetCorr.setSeed(rndSeed(event, event._allJets))
            self.jetmetCorr.correctJetAndMET(
                jets=event._allJets,
                lowPtJets=Collection(event, "CorrT1METJet"),
                met=event.met,
                rawMET=METObject(event, "RawPuppiMET"),
                rho=rho, genjets=Collection(event, 'GenJet') if self.isMC else None,
                isMC=self.isMC, runNumber=event.run)
            event._allJets = sorted(event._allJets, key=lambda x: x.pt, reverse=True)  # sort by pt after updating
            # correct AK15 fatjets
            self.fatjetCorr.setSeed(rndSeed(event, event._allAK15jets))
            self.fatjetCorr.correctJetAndMET(jets=event._allAK15jets, met=None, rho=rho,
                                             genjets=Collection(event, 'GenJetAK15') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            # correct AK15 subjets: do not sort after updating!!
            self.subjetCorr.setSeed(rndSeed(event, event.ak15Subjets))
            self.subjetCorr.correctJetAndMET(jets=event.ak15Subjets, met=None, rho=rho,
                                             genjets=Collection(event, 'GenSubJetAK15') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

        # link fatjet to (updated) subjets and recompute softdrop mass
        for idx, fj in enumerate(event._allAK15jets):
            fj.idx = idx
            fj.subjets = get_subjets(fj, event.ak15Subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = sumP4(*fj.subjets).M()
        event._allAK15jets = sorted(event._allAK15jets, key=lambda x: x.pt, reverse=True)  # sort by pt

        # select lepton-cleaned jets
        if self._doJetCleaning:
            event.fatjets = [fj for fj in event._allAK15jets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
                fj.jetId & 2) and closest(fj, event.looseLeptons)[1] >= self._jetConeSize]
            event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
                j.jetId & 4) and closest(j, event.looseLeptons)[1] >= 0.4]
        else:
            event.fatjets = [fj for fj in event._allAK15jets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
                fj.jetId & 2)]
            event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
                j.jetId & 4)]
        event.ht = sum([j.pt for j in event.ak4jets])

    def runModules(self, event):
        # needs event.selectedLeptons (muon SFs) and event._allJets (jet veto map)
        for mod in self._modules.values():
            assert mod.analyze(event)

    def loadGenHistory(self, event, fatjets):
        # gen matching
        if not self.isMC:
            return

        try:
            genparts = event.genparts
        except RuntimeError as e:
            genparts = Collection(event, "GenPart")
            for idx, gp in enumerate(genparts):
                if 'dauIdx' not in gp.__dict__:
                    gp.dauIdx = []
                if gp.genPartIdxMother >= 0:
                    mom = genparts[gp.genPartIdxMother]
                    if 'dauIdx' not in mom.__dict__:
                        mom.dauIdx = [idx]
                    else:
                        mom.dauIdx.append(idx)
            event.genparts = genparts

        def isHadronic(gp):
            if len(gp.dauIdx) == 0:
                return False
                # raise ValueError('Particle has no daughters!')
            for idx in gp.dauIdx:
                if abs(genparts[idx].pdgId) < 6:
                    return True
            return False

        def getFinal(gp):
            for idx in gp.dauIdx:
                dau = genparts[idx]
                if dau.pdgId == gp.pdgId:
                    return getFinal(dau)
            return gp

        lepGenTops = []
        hadGenTops = []
        hadGenWs = []
        hadGenZs = []
        hadGenHs = []

        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                continue
            if abs(gp.pdgId) == 6:
                for idx in gp.dauIdx:
                    dau = genparts[idx]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            hadGenTops.append(gp)
                        else:
                            lepGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                if isHadronic(gp):
                    hadGenWs.append(gp)
            elif abs(gp.pdgId) == 23:
                if isHadronic(gp):
                    hadGenZs.append(gp)
            elif abs(gp.pdgId) == 25:
                if isHadronic(gp):
                    hadGenHs.append(gp)

        for parton in itertools.chain(lepGenTops, hadGenTops):
            parton.daus = (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            parton.genW.daus = parton.daus[1:]
        for parton in itertools.chain(hadGenWs, hadGenZs, hadGenHs):
            parton.daus = (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])

        for fj in fatjets:
            fj.genH, fj.dr_H = closest(fj, hadGenHs)
            fj.genZ, fj.dr_Z = closest(fj, hadGenZs)
            fj.genW, fj.dr_W = closest(fj, hadGenWs)
            fj.genT, fj.dr_T = closest(fj, hadGenTops)
            fj.genLepT, fj.dr_LepT = closest(fj, lepGenTops)

    def evalTagger(self, jets):
        # same as NanoTrees vhTreeProducer._evalTagger (nanoV15)
        for j in jets:
            j.pn_Xbb = j.ParticleNetMD_probXbb
            j.pn_Xcc = j.ParticleNetMD_probXcc
            j.pn_Xqq = j.ParticleNetMD_probXqq
            j.pn_QCD = convert_prob(j, None, prefix='ParticleNetMD_prob')
            j.pn_XbbVsQCD = convert_prob(j, ['Xbb'], ['QCD'], prefix='pn_')
            j.pn_XccVsQCD = convert_prob(j, ['Xcc'], ['QCD'], prefix='pn_')
            j.pn_XccOrXqqVsQCD = convert_prob(j, ['Xcc', 'Xqq'], ['QCD'], prefix='pn_')
            j.part_HbbVsQCD = convert_prob(j, ['Xbb'], ['QCD'], prefix='ParTv3_prob')
            j.part_HccVsQCD = convert_prob(j, ['Xcc'], ['QCD'], prefix='ParTv3_prob')
            j.part_TVsQCD = convert_prob(j, ['TopbWqq', 'TopbWq'], ['QCD'], prefix='ParTv3_prob')

    def evalMassRegression(self, event, jets):
        # same as NanoTrees vhTreeProducer._evalMassRegression (nanoV15)
        for j in jets:
            j.masses = {
                'sdmass': j.msoftdrop,
                'regressed_mass': j.ParticleNetMD_mass,
                'ParT_resonanceMass': j.ParTv3_massCorrResonance * j.mass * (1. - j.rawFactor),
                'ParT_genericMass': j.ParTv3_massCorrGeneric * j.mass * (1. - j.rawFactor),
                'ParT_X2pMass': j.ParTv3_massCorrX2p * j.mass * (1. - j.rawFactor),
            }

    def fillBaseEventInfo(self, event):
        self.out.fillBranch("jetR", self._jetConeSize)
        self.out.fillBranch("year", self.year)
        self.out.fillBranch("lumiwgt", lumi_dict[self.year])

        met_filters = bool(
            event.Flag_goodVertices and
            event.Flag_globalSuperTightHalo2016Filter and
            event.Flag_HBHENoiseFilter and
            event.Flag_HBHENoiseIsoFilter and
            event.Flag_EcalDeadCellTriggerPrimitiveFilter and
            event.Flag_BadPFMuonFilter and
            event.Flag_BadPFMuonDzFilter and
            event.Flag_eeBadScFilter
        )
        # TODO: check which is used in the main analysis
        if self.year >= 2017:
            #if self.year in (2017, 2018):
            met_filters = met_filters and event.Flag_ecalBadCalibFilter
        self.out.fillBranch("passmetfilters", met_filters)

        # L1 prefire weights
        if self.year <= 2017:
            self.out.fillBranch("l1PreFiringWeight", event.L1PreFiringWeight_Nom)
            self.out.fillBranch("l1PreFiringWeightUp", event.L1PreFiringWeight_Up)
            self.out.fillBranch("l1PreFiringWeightDown", event.L1PreFiringWeight_Dn)
        else:
            self.out.fillBranch("l1PreFiringWeight", 1.0)
            self.out.fillBranch("l1PreFiringWeightUp", 1.0)
            self.out.fillBranch("l1PreFiringWeightDown", 1.0)

        self.out.fillBranch("nlep", len(event.looseLeptons))
        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("met", event.met.pt)
        self.out.fillBranch("metphi", event.met.phi)

    def fillFatJetInfo(self, event, fatjets):
        for idx in [1]:
            prefix = 'fj_%d_' % idx

            if len(fatjets) <= idx - 1:
                # fill zeros if fatjet fails probe selection
                for b in self.out._branches.keys():
                    if b.startswith(prefix):
                        self.out.fillBranch(b, 0)
                continue

            fj = fatjets[idx - 1]

            # fatjet kinematics
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            self.out.fillBranch(prefix + "tau32", fj.tau3 / fj.tau2 if fj.tau2 > 0 else 99)
            self.out.fillBranch(prefix + "nsubjets", len(fj.subjets))

            for name, val in fj.masses.items():
                self.out.fillBranch(prefix + name, val)

            # ParticleNet-MD
            self.out.fillBranch(prefix + "ParticleNetMD_Xbb", fj.pn_Xbb)
            self.out.fillBranch(prefix + "ParticleNetMD_Xcc", fj.pn_Xcc)
            self.out.fillBranch(prefix + "ParticleNetMD_Xqq", fj.pn_Xqq)
            self.out.fillBranch(prefix + "ParticleNetMD_QCD", fj.pn_QCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", fj.pn_XbbVsQCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", fj.pn_XccVsQCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XccOrXqqVsQCD", fj.pn_XccOrXqqVsQCD)

            # ParTv3
            self.out.fillBranch(prefix + "ParT_HbbVsQCD", fj.part_HbbVsQCD)
            self.out.fillBranch(prefix + "ParT_HccVsQCD", fj.part_HccVsQCD)
            self.out.fillBranch(prefix + "ParT_TVsQCD", fj.part_TVsQCD)
            self.out.fillBranch(prefix + "ParT_QCD", fj.ParTv3_probQCD)
            for score in self._part_scores:
                self.out.fillBranch(prefix + "ParT_" + score, getattr(fj, 'ParTv3_prob' + score))

            # matching variables
            if self.isMC:
                self.out.fillBranch(prefix + "nbhadrons", fj.nBHadrons)
                self.out.fillBranch(prefix + "nchadrons", fj.nCHadrons)

                # info of the closest hadGenH
                self.out.fillBranch(prefix + "dr_H", fj.dr_H)
                self.out.fillBranch(prefix + "dr_H_daus",
                                    max([deltaR(fj, dau) for dau in fj.genH.daus]) if fj.genH else 99)
                self.out.fillBranch(prefix + "H_pt", fj.genH.pt if fj.genH else -1)
                self.out.fillBranch(prefix + "H_decay", abs(fj.genH.daus[0].pdgId) if fj.genH else 0)

                # info of the closest hadGenZ
                self.out.fillBranch(prefix + "dr_Z", fj.dr_Z)
                self.out.fillBranch(prefix + "dr_Z_daus",
                                    max([deltaR(fj, dau) for dau in fj.genZ.daus]) if fj.genZ else 99)
                self.out.fillBranch(prefix + "Z_pt", fj.genZ.pt if fj.genZ else -1)
                self.out.fillBranch(prefix + "Z_decay", abs(fj.genZ.daus[0].pdgId) if fj.genZ else 0)

                # info of the closest hadGenW
                self.out.fillBranch(prefix + "dr_W", fj.dr_W)
                self.out.fillBranch(prefix + "dr_W_daus",
                                    max([deltaR(fj, dau) for dau in fj.genW.daus]) if fj.genW else 99)
                self.out.fillBranch(prefix + "W_pt", fj.genW.pt if fj.genW else -1)
                self.out.fillBranch(prefix + "W_decay", max([abs(d.pdgId) for d in fj.genW.daus]) if fj.genW else 0)

                # info of the closest hadGenTop
                drwq1, drwq2 = [deltaR(fj, dau) for dau in fj.genT.genW.daus] if fj.genT else [99, 99]
                wq1_pdgId, wq2_pdgId = [dau.pdgId for dau in fj.genT.genW.daus] if fj.genT else [0, 0]
                if drwq1 < drwq2:
                    drwq1, drwq2 = drwq2, drwq1
                    wq1_pdgId, wq2_pdgId = wq2_pdgId, wq1_pdgId
                self.out.fillBranch(prefix + "dr_T", fj.dr_T)
                self.out.fillBranch(prefix + "dr_T_b", deltaR(fj, fj.genT.genB) if fj.genT else 99)
                self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
                self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
                self.out.fillBranch(prefix + "T_Wq_max_pdgId", wq1_pdgId)
                self.out.fillBranch(prefix + "T_Wq_min_pdgId", wq2_pdgId)
                self.out.fillBranch(prefix + "T_pt", fj.genT.pt if fj.genT else -1)