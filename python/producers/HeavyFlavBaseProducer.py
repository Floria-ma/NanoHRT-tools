import os
import itertools
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.utils import deltaR, closest, polarP4, sumP4, get_subjets, corrected_svmass, configLogger
from ..helpers.xgbHelper import XGBEnsemble
from ..helpers.nnHelper import convert_prob, ensemble
from ..helpers.jetmetCorrector import JetMETCorrector, rndSeed, JetIdCalculator, AK15patcher
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
        self._channel = channel  # 'qcd', 'photon', 'inclusive', 'muon'
        self.year = int(kwargs['year'])
        self.jetType = kwargs.get('jetType', 'ak15').lower()
        self._jmeSysts = {'jec': False, 'jes': None, 'jes_source': '', 'jes_uncertainty_file_prefix': '',
                          'jer': None, 'jmr': None, 'met_unclustered': None, 'smearMET': True, 'applyHEMUnc': False,
                          'jesr_extra_br': True}
        self._opts = {'sfbdt_threshold': -99,
                      'run_tagger': False, 'tagger_versions': ['V02b', 'V02c', 'V02d'],
                      'run_mass_regression': False, 'mass_regression_versions': ['V01a', 'V01b', 'V01c'],
                      'WRITE_CACHE_FILE': False}
        for k in kwargs:
            if k in self._jmeSysts:
                self._jmeSysts[k] = kwargs[k]
            else:
                self._opts[k] = kwargs[k]
        self._needsJMECorr = any([self._jmeSysts['jec'], self._jmeSysts['jes'],
                                  self._jmeSysts['jer'], self._jmeSysts['jmr'],
                                  self._jmeSysts['met_unclustered'], self._jmeSysts['applyHEMUnc']])
        self._doJetCleaning = True

        logger.info('Running %s channel for %s jets with JME systematics %s, other options %s',
                    self._channel, self.jetType, str(self._jmeSysts), str(self._opts))
        
        if self.jetType == 'ak8':
            self._jetConeSize = 0.8
            self._fj_name = 'FatJet'
            #self._sj_name = 'SubJet'
            self._fj_gen_name = 'GenJetAK8'
        elif self.jetType == 'ak15':
            self._jetConeSize = 1.5
            self._fj_name = 'AK15Puppi'
            #self._sj_name = 'AK15PuppiSubJet'
            self._fj_gen_name = 'GenJetAK15'
        else:
            raise RuntimeError('Jet type %s is not recognized!' % self.jetType)

        self._fill_sv = self._channel in ('qcd', 'photon', 'higgs', 'inclusive') and self._opts['sfbdt_threshold'] > -99
        
        self.jetidcalc = JetIdCalculator(year=self._year, jetType="AK4PUPPI")
        self.fatjetidcalc = JetIdCalculator(year=self._year, jetType="AK8PUPPI") # does not work because some branches are missing
        self.ak15patcher = AK15patcher(year=self._year)

        if self._needsJMECorr:
            self.jetmetCorr = JetMETCorrector(year=self.year, jetType="AK4PFchs", **self._jmeSysts)
            self.fatjetCorr = JetMETCorrector(year=self.year, jetType="AK8PFPuppi", **self._jmeSysts)
            self.subjetCorr = JetMETCorrector(year=self.year, jetType="AK4PFPuppi", **self._jmeSysts)

        #check the pnet WP for 2024
        self.PNet_WP_M = {
            2015: 0,
            2016: 0,
            2017: 0,
            2018: 0,
            2021: 0.245,
            2022: 0.2605,
            2023: 0.1917,
            2024: 0.1919
        }[self.year]

    def beginJob(self):
        self.jetidcalc.beginJob()
        self.fatjetidcalc.beginJob()
        if self._needsJMECorr:
            self.jetmetCorr.beginJob()
            self.fatjetCorr.beginJob()
            self.subjetCorr.beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.rho_branch_name = 'Rho_fixedGridRhoFastjetAll' if bool(
            inputTree.GetBranch('Rho_fixedGridRhoFastjetAll')) else 'fixedGridRhoFastjetAll'
        #self.hasParticleNetProb = bool(inputTree.GetBranch(self._fj_name + '_ParticleNet_probHbb'))

        # remove all possible h5 cache files
        for f in os.listdir('.'):
            if f.endswith('.h5'):
                os.remove(f)

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
        if self.isMC and self._jmeSysts['jesr_extra_br']:
            # HT with JES/JER correction
            self.out.branch("ht_jesUncFactorUp", "F")
            self.out.branch("ht_jesUncFactorDn", "F")
            self.out.branch("ht_jerSmearFactorUp", "F")
            self.out.branch("ht_jerSmearFactorDn", "F")

        # Large-R jets
        for idx in ([1, 2] if self._channel in ['qcd', 'mutagged'] else [1]):
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            #self.out.branch(prefix + "is_qualified", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "regressed_mass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")
            
            self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
            self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
            self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
            #self.out.branch(prefix + "ParticleNetMD_QCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccOrXqqVsQCD", "F")

            self.out.branch(prefix + "ParT_Hbb", "F")
            self.out.branch(prefix + "ParT_Hcc", "F")
            self.out.branch(prefix + "ParT_Hqq", "F")
            self.out.branch(prefix + "ParT_HbbVsQCD", "F")
            self.out.branch(prefix + "ParT_HccVsQCD", "F")
            self.out.branch(prefix + "ParT_TVsQCD", "F")

            self.out.branch(prefix + "ParticleNet_mass", "F")
            self.out.branch(prefix + "btagDDBvLV2", "F")
            self.out.branch(prefix + "btagDDCvBV2", "F")
            self.out.branch(prefix + "btagDDCvLV2", "F")
            self.out.branch(prefix + "btagDeepB", "F")
            self.out.branch(prefix + "btagHbb", "F")

            if self._opts['run_tagger']:
                self.out.branch(prefix + "origParticleNetMD_XccVsQCD", "F")
                self.out.branch(prefix + "origParticleNetMD_XbbVsQCD", "F")

            # matching variables
            if self.isMC:
                self.out.branch(prefix + "nbhadrons", "I")
                self.out.branch(prefix + "nchadrons", "I")
                self.out.branch(prefix + "partonflavour", "I")
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
 

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        #for mod in self._modules.values():
        #    mod.endFile(inputFile, outputFile, inputTree, wrappedOutputTree)
        
        if self._opts['run_tagger'] and self._opts['WRITE_CACHE_FILE']:
            for p in self.pnTaggers:
                p.update_cache()

        if self._opts['run_mass_regression'] and self._opts['WRITE_CACHE_FILE']:
            for p in self.pnMassRegressions:
                p.update_cache()

        # remove all h5 cache files
        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            for f in os.listdir('.'):
                if f.endswith('.h5'):
                    os.remove(f)

    def selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []  # used for jet lepton cleaning & lepton counting

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            # follow the main analysis:
            if el.pt > 20 and abs(el.eta) < 2.5 and el.mvaIso_WP90: #el.mvaFall17V2noIso_WP90 and el.miniPFRelIso_all < 0.4: #THESE ID IS NOT SAVED IN OUR CUSTOM NANOS
                event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            # follow the main analysis, looseId :
            if mu.pt > 20 and abs(mu.eta) < 2.4 and mu.looseId and mu.pfRelIso04_all < 0.25:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

    def correctJetsAndMET(self, event):
        # correct Jets and MET
        event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)
        event._allAK15jets = Collection(event, "AK15Puppi")
        event._allAK8jets = Collection(event, "FatJet")        
        event._allJets = Collection(event, "Jet")
        event.ak15Subjets = Collection(event, "AK15PuppiSubJet")  # do not sort after updating!!
        event.met = METObject(event, "PuppiMET")
        event.trigobjs = Collection(event, "TrigObj")

        self.ak15patcher.CalculateMissingVars(event._allJets, event._allAK8jets, event._allAK15jets)
        self.jetidcalc.CalculateJetId(event._allJets)
        self.jetidcalc.CalculateJetId(event._allAK15jets) 

        if self._needsJMECorr:
            rho = getattr(event, self.rho_branch_name)
            #rho = event.fixedGridRhoFastjetAll
            # correct AK4 jets and MET
            self.jetmetCorr.setSeed(rndSeed(event, event._allJets))
            self.jetmetCorr.correctJetAndMET(jets=event._allJets, lowPtJets=Collection(event, "CorrT1METJet"),
                                             met=event.met, rawMET=METObject(event, "RawMET"),
                                             defaultMET=METObject(event, "MET"),
                                             rho=rho, genjets=Collection(event, 'GenJet') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            event._allJets = sorted(event._allJets, key=lambda x: x.pt, reverse=True)  # sort by pt after updating

            # correct fatjets
            self.fatjetCorr.setSeed(rndSeed(event, event._allFatJets))
            self.fatjetCorr.correctJetAndMET(jets=event._allFatJets, met=None, rho=rho,
                                             genjets=Collection(event, self._fj_gen_name) if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            # correct subjets
            self.subjetCorr.setSeed(rndSeed(event, event.subjets))
            self.subjetCorr.correctJetAndMET(jets=event.subjets, met=None, rho=rho,
                                             #genjets=Collection(event, self._sj_gen_name) if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

        # jet mass resolution smearing
        if self.isMC and self._jmeSysts['jmr']:
            raise NotImplementedError

        # link fatjet to subjets and recompute softdrop mass
        for idx, fj in enumerate(event._allFatJets):
            fj.idx = idx
            #fj.is_qualified = True
            fj.subjets = get_subjets(fj, event.subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = sumP4(*fj.subjets).M()
        event._allFatJets = sorted(event._allFatJets, key=lambda x: x.pt, reverse=True)  # sort by pt

        # select lepton-cleaned jets
        if self._doJetCleaning:
            event.fatjets = [fj for fj in event._allFatJets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
                fj.jetId & 2) and closest(fj, event.looseLeptons)[1] >= self._jetConeSize]
            event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
                j.jetId & 4) and closest(j, event.looseLeptons)[1] >= 0.4]
        else:
            event.fatjets = [fj for fj in event._allFatJets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
                fj.jetId & 2)]
            event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
                j.jetId & 4)]
        event.ht = sum([j.pt for j in event.ak4jets])
        if self.isMC and self._jmeSysts['jesr_extra_br']:
            event.ht_jesUncFactorUp = sum([j.pt * j.jesUncFactorUp for j in event.ak4jets])
            event.ht_jesUncFactorDn = sum([j.pt * j.jesUncFactorDn for j in event.ak4jets])
            event.ht_jerSmearFactorUp = sum([j.pt * j.jerSmearFactorUp for j in event.ak4jets])
            event.ht_jerSmearFactorDn = sum([j.pt * j.jerSmearFactorDn for j in event.ak4jets])
    
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

    #TODO: check which mass variable is used in the main analysis
    def evalMassRegression(self, event, jets):
        for j in jets:
            if self._opts['run_mass_regression']:
                outputs = [p.predict_with_cache(self.tagInfoMaker, event.idx, j.idx, j) for p in self.pnMassRegressions]
                j.regressed_mass = ensemble(outputs, np.median)['mass']
            else:
                try:
                    j.regressed_mass = j.ParticleNetMD_mass
                except RuntimeError:
                    j.regressed_mass = 0

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

    def _get_filler(self, obj):

        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)

        return filler

    def fillFatJetInfo(self, event, fatjets):
        for idx in ([1, 2] if self._channel in ['qcd', 'mutagged'] else [1]):
            prefix = 'fj_%d_' % idx

            if len(fatjets) <= idx - 1 or not fatjets[idx - 1].is_qualified:
                # fill zeros if fatjet fails probe selection
                for b in self.out._branches.keys():
                    if b.startswith(prefix):
                        self.out.fillBranch(b, 0)
                continue

            fj = fatjets[idx - 1]

            # fatjet kinematics
            #self.out.fillBranch(prefix + "is_qualified", fj.is_qualified)
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "regressed_mass", fj.regressed_mass)
            #self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            #self.out.fillBranch(prefix + "tau32", fj.tau3 / fj.tau2 if fj.tau2 > 0 else 99)
            try:
                # need to be synchronized withh the main analysis
                self.out.fillBranch(prefix + "ParT_HccVsQCD", convert_prob(fj, ['Hcc'], None, prefix='ParT_prob'))
                self.out.fillBranch(prefix + "ParT_HbbVsQCD", convert_prob(fj, ['Hbb'], None, prefix='ParT_prob'))
                self.out.fillBranch(prefix + "ParT_TVsQCD", convert_prob(fj, ['TopbWcs', 'TopbWqq', 'TopbWc', 'TopbWs', 'TopbWq','TopWcs', 'TopWqq'], None, prefix='ParT_prob'))

            except RuntimeError:
                self.out.fillBranch(prefix + "ParT_HccVsQCD", -1)
                self.out.fillBranch(prefix + "ParT_HbbVsQCD", -1)
                self.out.fillBranch(prefix + "ParT_TopVsQCD", -1)

            # ParticleNet-MD
            self.out.fillBranch(prefix + "ParticleNetMD_Xbb", fj.pn_Xbb)
            self.out.fillBranch(prefix + "ParticleNetMD_Xcc", fj.pn_Xcc)
            self.out.fillBranch(prefix + "ParticleNetMD_Xqq", fj.pn_Xqq)
            #self.out.fillBranch(prefix + "ParticleNetMD_QCD", fj.pn_QCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", fj.pn_XbbVsQCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", fj.pn_XccVsQCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XccOrXqqVsQCD", fj.pn_XccOrXqqVsQCD)

            #ParT
            self.out.fillBranch(prefix + "ParT_Hbb", fj.pt_Hbb)
            self.out.fillBranch(prefix + "ParT_Hcc", fj.pt_Hcc)
            self.out.fillBranch(prefix + "ParT_Hqq", fj.pt_Hqq)

            if self._opts['run_tagger']:
                self.out.fillBranch(prefix + "origParticleNetMD_XccVsQCD",
                                    convert_prob(fj, ['Xcc'], None, prefix='ParticleNetMD_prob'))
                self.out.fillBranch(prefix + "origParticleNetMD_XbbVsQCD",
                                    convert_prob(fj, ['Xbb'], None, prefix='ParticleNetMD_prob'))

            try:
                self.out.fillBranch(prefix + "ParticleNet_mass", fj.particleNet_mass)
            except RuntimeError:
                self.out.fillBranch(prefix + "ParticleNet_mass", -1)
            try:
                self.out.fillBranch(prefix + "btagDDBvLV2", fj.btagDDBvLV2)
                self.out.fillBranch(prefix + "btagDDCvBV2", fj.btagDDCvBV2)
                self.out.fillBranch(prefix + "btagDDCvLV2", fj.btagDDCvLV2)
                self.out.fillBranch(prefix + "btagDeepB", fj.btagDeepB)
                self.out.fillBranch(prefix + "btagHbb", fj.btagHbb)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagDDBvLV2", -1)
                self.out.fillBranch(prefix + "btagDDCvBV2", -1)
                self.out.fillBranch(prefix + "btagDDCvLV2", -1)
                self.out.fillBranch(prefix + "btagDeepB", -1)
                self.out.fillBranch(prefix + "btagHbb", -1)

            # matching variables
            if self.isMC:

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

                if self._jmeSysts['jesr_extra_br']:
                    self.out.fillBranch(prefix + "jesUncFactorUp", fj.jesUncFactorUp)
                    self.out.fillBranch(prefix + "jesUncFactorDn", fj.jesUncFactorDn)
                    self.out.fillBranch(prefix + "jerSmearFactorUp", fj.jerSmearFactorUp)
                    self.out.fillBranch(prefix + "jerSmearFactorDn", fj.jerSmearFactorDn)

                if self.isMC and idx==(2 if self._channel == 'qcd' else 1) and self._channel != 'higgs':
                    for hadtype in ['b', 'c']:
                        for hadidx in [1, 2]:
                            gp = fj.genBhadron[hadidx - 1] if hadtype=='b' else fj.genChadron[hadidx - 1]
                            fill_gp = self._get_filler(gp)  # wrapper, fill default value if sv=None
                            fill_gp(prefix + "gen{}hadron{}_pt".format(hadtype, hadidx), gp.pt)
                            fill_gp(prefix + "gen{}hadron{}_eta".format(hadtype, hadidx), gp.eta)
                            fill_gp(prefix + "gen{}hadron{}_phi".format(hadtype, hadidx), gp.phi)
                            fill_gp(prefix + "gen{}hadron{}_mass".format(hadtype, hadidx), gp.mass)
                            fill_gp(prefix + "gen{}hadron{}_pdgId".format(hadtype, hadidx), gp.pdgId)

                if self.isMC and self._channel != 'higgs':
                    self.out.fillBranch(prefix + "npart", fj.npart)
                    self.out.fillBranch(prefix + "nbpart", fj.nbpart)
                    self.out.fillBranch(prefix + "ncpart", fj.ncpart)
                    self.out.fillBranch(prefix + "ngpart", fj.ngpart)
                    self.out.fillBranch(prefix + "part_sumpt", fj.part_sumpt)
                    self.out.fillBranch(prefix + "bpart_sumpt", fj.bpart_sumpt)
                    self.out.fillBranch(prefix + "cpart_sumpt", fj.cpart_sumpt)
                    self.out.fillBranch(prefix + "gpart_sumpt", fj.gpart_sumpt)
                    self.out.fillBranch(prefix + "npart50", fj.npart50)
                    self.out.fillBranch(prefix + "nbpart50", fj.nbpart50)
                    self.out.fillBranch(prefix + "ncpart50", fj.ncpart50)
                    self.out.fillBranch(prefix + "ngpart50", fj.ngpart50)
                    self.out.fillBranch(prefix + "part50_sumpt", fj.part50_sumpt)
                    self.out.fillBranch(prefix + "bpart50_sumpt", fj.bpart50_sumpt)
                    self.out.fillBranch(prefix + "cpart50_sumpt", fj.cpart50_sumpt)
                    self.out.fillBranch(prefix + "gpart50_sumpt", fj.gpart50_sumpt)