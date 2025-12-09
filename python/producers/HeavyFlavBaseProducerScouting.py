import os
import logging
import numpy as np
import itertools
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.nnHelper import convert_prob, ensemble
from ..helpers.jetmetCorrector import JetMETCorrector, rndSeed
from ..helpers.triggerHelper import passTrigger
from ..helpers.utils import deltaR, closest, polarP4, configLogger, deltaR2, deltaPhi

logger = logging.getLogger("nano")
configLogger("nano", loglevel=logging.INFO)


# define luminosity per year.
# note: for now, this is only meant for 2024 scouting samples.
lumi_dict = {2024: 109.08}


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


class HeavyFlavBaseProducerScouting(Module, object):

    def __init__(self, channel, **kwargs):
        super(HeavyFlavBaseProducerScouting, self).__init__()

        # basic settings
        self._channel = channel
        # note: for now, this is only meant for channel = "muon"
        self.year = int(kwargs["year"])
        # note: for now, this is only meant for 2024 samples
        self.jetType = kwargs.get("jetType", "scouting").lower()
        self._jmeSysts = {'jec': False,
                          'jes': None,
                          'jes_source': '',
                          'jes_uncertainty_file_prefix': '',
                          'jer': None,
                          'jmr': None,
                          'met_unclustered': None,
                          'smearMET': True,
                          'applyHEMUnc': False,
                          'jesr_extra_br': True}
        
        #self._opts = {
        #    "WRITE_CACHE_FILE": False,
        #    'mass_range': (40, 250),
        #}
        self._opts = {'sfbdt_threshold': -99,
                      'run_tagger': False, 'tagger_versions': ['V02b', 'V02c', 'V02d'],
                      'run_mass_regression': False, 'mass_regression_versions': ['V01a', 'V01b', 'V01c'],
                      'WRITE_CACHE_FILE': False}

        for k in kwargs:
            if k in self._jmeSysts:
                self._jmeSysts[k] = kwargs[k]
            else:
                self._opts[k] = kwargs[k]

        #self._needsJMECorr = False  
        self._needsJMECorr = any([self._jmeSysts['jec'], self._jmeSysts['jes'],
                                  self._jmeSysts['jer'], self._jmeSysts['jmr'],
                                  self._jmeSysts['met_unclustered'], self._jmeSysts['applyHEMUnc']])
        #logger.info('Running %s channel for year %s with JME systematics %s, other options %s',
        #            self._channel, str(self._year), str(self._jmeSysts), str(self._opts))

        self._doJetCleaning = True

        #if self.jetType == 'ak8':
        # pfjet collection
        self._ak4_name = "ScoutingPFJet"
        # fatjet collection
        self._fatjet_name = "ScoutingFatPFJetRecluster"
        self._jetConeSize = 0.8 # we use ak8 for now

        logger.info(
            "HeavyFlavBaseProducerScouting: channel=%s, year=%s, ak4=%s, fatjets=%s",
            self._channel, self.year, self._ak4_name, self._fatjet_name
        )

    def beginJob(self):
        # nothing heavy to initialize for scouting
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # MC or data
        self.isMC = bool(inputTree.GetBranch("genWeight"))

        # try to get an available rho branch
        if inputTree.GetBranch("ScoutingRho_fixedGridRhoFastjetAll"):
            self.rho_branch_name = "ScoutingRho_fixedGridRhoFastjetAll"
        else:
            self.rho_branch_name = None   # no rho in scouting

        self.out = wrappedOutputTree

        # NOTE: branch names must start with a lower case letter
        # check keep_and_drop_output.txt
        self.out.branch("year", "I")
        self.out.branch("lumiwgt", "F")
        self.out.branch("jetR", "F")
        self.out.branch("passmetfilters", "O")
        self.out.branch("nlep", "I")
        self.out.branch("ht", "F")
        self.out.branch("met", "F")
        self.out.branch("met_phi", "F")
        #self.out.branch("passTrigPFHTScouting", "O")

        # prefiring weight branches
        # (not really defined for 2024 and/or scouting,
        # kept only for consistency with other years)
        self.out.branch("l1PreFiringWeight", "F")
        self.out.branch("l1PreFiringWeightUp", "F")
        self.out.branch("l1PreFiringWeightDown", "F")
        
        # Large-R jets
        for idx in ([1, 2] if self._channel in ['qcd', 'mutagged'] else [1]):
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            #self.out.branch(prefix + "is_qualified", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "mass", "F")
            #self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "ParT_resonanceMass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")
            #self.out.branch(prefix + "btagcsvv2", "F")
            #self.out.branch(prefix + "btagjp", "F")

            # ParticleNet (scouting)
            self.out.branch(prefix + "PN_Hbb", "F")
            self.out.branch(prefix + "PN_Hcc", "F")
            self.out.branch(prefix + "PN_Hqq", "F")
            self.out.branch(prefix + "PN_QCD", "F")
            self.out.branch(prefix + "PN_mass", "F")

            # ParT (scoutGlobalParT)
            self.out.branch(prefix + "scoutGloParT_Xbb", "F")
            self.out.branch(prefix + "scoutGloParT_Xcc", "F")
            self.out.branch(prefix + "scoutGloParT_Xqq", "F")
            self.out.branch(prefix + "scoutGloParT_Xgg", "F")
            self.out.branch(prefix + "scoutGloParT_Xbc", "F")
            self.out.branch(prefix + "scoutGloParT_Xbs", "F")
            self.out.branch(prefix + "scoutGloParT_Xcs", "F")
            self.out.branch(prefix + "scoutGloParT_Xss", "F")
            self.out.branch(prefix + "scoutGloParT_tauhtaue", "F")
            self.out.branch(prefix + "scoutGloParT_tauhtauh", "F")
            self.out.branch(prefix + "scoutGloParT_tauhtaum", "F")
            self.out.branch(prefix + "scoutGloParT_QCD", "F")

            # VSQCD scores
            self.out.branch(prefix + "PN_HbbVsQCD", "F")
            self.out.branch(prefix + "PN_HccVsQCD", "F")
            self.out.branch(prefix + "PN_HqqVsQCD", "F")

            self.out.branch(prefix + "scoutGloParT_HbbVsQCD", "F")
            self.out.branch(prefix + "scoutGloParT_HccVsQCD", "F")
            self.out.branch(prefix + "scoutGloParT_HqqVsQCD", "F")
            self.out.branch(prefix + "scoutGloParT_twoprong", "F")
            self.out.branch(prefix + "scoutGloParT_HbqVsQCD",  "F")
            self.out.branch(prefix + "scoutGloParT_HcsVsQCD",  "F")
            self.out.branch(prefix + "scoutGloParT_HbbVsHcc", "F")
            self.out.branch(prefix + "scoutGloParT_HVsQCD",  "F")
            self.out.branch(prefix + "scoutGloParT_ZVsQCD",  "F")
            self.out.branch(prefix + "scoutGloParT_WVsQCD",  "F")

            # massCorr branches
            self.out.branch(prefix + "scoutGloParT_massCorrGeneric", "F")
            self.out.branch(prefix + "scoutGloParT_massCorrGenericW2p", "F")
            self.out.branch(prefix + "scoutGloParT_massCorrGenericX2p", "F")
            self.out.branch(prefix + "scoutGloParT_massCorrResonance", "F")

            # gen-matching variables
            if self.isMC:
                # info of the closest hadGenTop
                self.out.branch(prefix + "dr_T", "F")
                self.out.branch(prefix + "dr_T_b", "F")
                self.out.branch(prefix + "dr_T_Wq_max", "F")
                self.out.branch(prefix + "dr_T_Wq_min", "F")
                self.out.branch(prefix + "T_Wq_max_pdgId", "I")
                self.out.branch(prefix + "T_Wq_min_pdgId", "I")
                #self.out.branch(prefix + "T_pt", "F")
                self.out.branch(prefix + "dr_W_daus", "F")
                self.out.branch(prefix + "dr_W", "F")
                #self.out.branch(prefix + "W_pt", "F")
                self.out.branch(prefix + "W_decay", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # nothing to clean up
        pass

    def selectLeptons(self, event):
        # do loose lepton selection
        # (typically used for jet lepton cleaning & lepton counting)
        
        event.looseLeptons = []

        # electrons
        electrons = Collection(event, "ScoutingElectron")
        for el in electrons:
            el.etaSC = el.eta #+ el.deltaEtaSC
            if el.pt > 10 and abs(el.eta) < 2.5:# and el.mvaIso_WP90:
                event.looseLeptons.append(el)

        # muons
        muons = Collection(event, "ScoutingMuonVtx")
        for mu in muons:
            if mu.pt > 10 and abs(mu.eta) < 2.4:
                event.looseLeptons.append(mu)

        # sort by pt
        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)
    
    def correctJetAndMET(self, event):
        # initialize jets and MET
       
        #event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)

        # scoutingPFJet
        try:
            event._allJets = Collection(event, self._ak4_name)
        except RuntimeError:
            logger.error("ScoutingPFJet collection not found in event")
            event._allJets = []

        # scoutingFatPFJetRecluster
        try:
            event._allFatJets = Collection(event, self._fatjet_name)
        except RuntimeError:
            logger.error("ScoutingFatPFJetRecluster collection not found in event")
            event._allFatJets = []

        # MET
        try:
            event.met = METObject(event, "ScoutingMET")
        except RuntimeError:
            class _DummyMET:
                pt = 0.0
                phi = 0.0

                def p4(self):
                    return polarP4(self, eta=None, mass=None)

            logger.error("ScoutingMET not found in event, using dummy MET instead.")
            event.met = _DummyMET()

        # sort jets by pt
        try:
            event._allJets = sorted(event._allJets, key=lambda x: x.pt, reverse=True)
        except Exception: pass

        # sort fat jets by pt
        try:
            for idx, fj in enumerate(event._allFatJets):
                fj.idx = idx
                fj.is_qualified = True
            event._allFatJets = sorted(event._allFatJets, key=lambda x: x.pt, reverse=True)
        except Exception: pass
        
        # select jets
        if self._doJetCleaning:
            event.ak4jets = [
                j for j in event._allJets
                if j.pt > 25
                and abs(j.eta) < 2.4
                and closest(j, getattr(event, "looseLeptons", []))[1] >= 0.4
            ]
        else:
            event.ak4jets = [
                j for j in event._allJets
                if j.pt > 25 and abs(j.eta) < 2.4
            ]

        # select fatjets
        if self._doJetCleaning:
            event.fatjets = [
                fj for fj in event._allFatJets
                if fj.pt > 200
                and abs(fj.eta) < 2.4
                and closest(fj, getattr(event, "looseLeptons", []))[1] >= self._jetConeSize
            ]
        else:
            event.fatjets = [
                fj for fj in event._allFatJets
                if fj.pt > 200 and abs(fj.eta) < 2.4
            ]

        # HT = scalar sum of selected AK4 jet pT
        event.ht = float(sum(j.pt for j in event.ak4jets))

    def evalMassRegression(self, jets):
        for j in jets:
            j.ParT_resonanceMass = j.scoutGlobalParT_massCorrResonance * j.mass
            j.sdmass = j.msoftdrop
            #j.masses = {
            #    'sdmass': j.msoftdrop,
                #'regressed_mass': j.particleNet_mass,
            #    'ParT_resonanceMass': j.scoutGlobalParT_massCorrResonance * j.mass,# * (1. - j.rawFactor),
                #'ParT_visibleMass': j.ParT_visiableMassCorr * j.mass * (1. - j.rawFactor),
            #}

    def fillBaseEventInfo(self, event):
        # jet radius (use reclustered fatjet R)
        self.out.fillBranch("jetR", self._jetConeSize)

        # year and lumi weight
        self.out.fillBranch("year", self.year)
        self.out.fillBranch("lumiwgt", lumi_dict.get(self.year, 1.0))

        # MET filters
        def _safe_flag(name, default=True):
            try:
                return bool(getattr(event, name))
            except AttributeError:
                return default

        met_filters = (
            _safe_flag("Flag_goodVertices") and
            _safe_flag("Flag_globalSuperTightHalo2016Filter") and
            _safe_flag("Flag_HBHENoiseFilter") and
            _safe_flag("Flag_HBHENoiseIsoFilter") and
            _safe_flag("Flag_EcalDeadCellTriggerPrimitiveFilter") and
            _safe_flag("Flag_BadPFMuonFilter") and
            _safe_flag("Flag_BadPFMuonDzFilter") and
            _safe_flag("Flag_eeBadScFilter")
        )

        if self.year >= 2017:
            met_filters = met_filters and _safe_flag("Flag_ecalBadCalibFilter")
        
        self.out.fillBranch("passmetfilters", met_filters)
        # trigger
        #self.out.fillBranch("passTrigPFHTScouting", passTrigger(event,["DST_PFScouting_JetHT"]))

        # fill prefire weight branches with dummy values
        self.out.fillBranch("l1PreFiringWeight", 1.0)
        self.out.fillBranch("l1PreFiringWeightUp", 1.0)
        self.out.fillBranch("l1PreFiringWeightDown", 1.0)

        self.out.fillBranch("nlep", len(event.looseLeptons))
        self.out.fillBranch("ht", event.ht)

        try:
            self.out.fillBranch("met", event.met.pt)
            self.out.fillBranch("met_phi", event.met.phi)
        except Exception:
            self.out.fillBranch("met", 0.0)
            self.out.fillBranch("met_phi", 0.0)

    def loadGenHistory(self, event, fatjets):
        if not self.isMC:
            return

        # use cleaned fatjets if available
        #if fatjets is None or len(fatjets) < 1:
        #    return False

        ak8  = fatjets[0]
        #ak8_2 = fatjets[1] if len(fatjets) > 1 else None

        # build genparts with dauIdx
        try:
            genparts = event.genparts
        except RuntimeError:
            genparts = Collection(event, "GenPart")
            for idx, gp in enumerate(genparts):
                if "dauIdx" not in gp.__dict__:
                    gp.dauIdx = []
                if gp.genPartIdxMother >= 0:
                    mom = genparts[gp.genPartIdxMother]
                    if "dauIdx" not in mom.__dict__:
                        mom.dauIdx = [idx]
                    else:
                        mom.dauIdx.append(idx)
            event.genparts = genparts

        def isHadronic(gp):
            return any(abs(genparts[i].pdgId) < 6 for i in getattr(gp, "dauIdx", []))

        def getFinal(gp):
            for i in getattr(gp, "dauIdx", []):
                dau = genparts[i]
                if dau.pdgId == gp.pdgId:
                    return getFinal(dau)
            return gp

        lepGenTops, hadGenTops, hadGenWs = [], [], []

        for gp in genparts:
            # maybe relax this completely for scouting
            if hasattr(gp, "statusFlags") and gp.statusFlags & (1 << 13) == 0:
                continue

            if abs(gp.pdgId) == 6:
                for i in gp.dauIdx:
                    dau = genparts[i]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            hadGenTops.append(gp)
                        else:
                            lepGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24 and isHadronic(gp):
                hadGenWs.append(gp)

        def get_daughters(parton):
            if parton is None:
                return []
            if abs(parton.pdgId) == 6 and hasattr(parton, "genW") and hasattr(parton, "genB"):
                W = parton.genW
                if getattr(W, "dauIdx", []):
                    return [parton.genB] + [genparts[i] for i in W.dauIdx[:2]]
            elif abs(parton.pdgId) in (23, 24, 25) and getattr(parton, "dauIdx", []):
                return [genparts[i] for i in parton.dauIdx[:2]]
            return []
        
        prefix = "fj_1_"  # assuming ak8 is fj_1
        
        genW, dr_W = closest(ak8, hadGenWs)

        # info of the closest hadGenW
        wdecay_ = max([abs(d.pdgId) for d in get_daughters(genW)], default=0) if genW else 0
        self.out.fillBranch(prefix + "dr_W", dr_W)
        self.out.fillBranch(prefix + "dr_W_daus",
                            max([deltaR(ak8, dau) for dau in get_daughters(genW)]) if genW else 99)
        self.out.fillBranch(prefix + "W_decay", wdecay_)

        # sort tops by deltaR to ak8
        genT, dr_T = closest(ak8, hadGenTops)
        #hadGenTops.sort(key=lambda x: deltaR2(x, ak8))
        #t = genT if len(genT) else None
        self.out.fillBranch(prefix + "dr_T_b", deltaR(ak8, genT.genB) if len(hadGenTops) else 99)

        drwq1, drwq2 = [deltaR(ak8, dau) for dau in get_daughters(
                genT.genW)] if len(hadGenTops) else [99, 99]
        wq1_pdgId, wq2_pdgId = [dau.pdgId for dau in get_daughters(genT.genW)] if len(hadGenTops) else [0, 0]
        if drwq1 < drwq2:
                drwq1, drwq2 = drwq2, drwq1
                wq1_pdgId, wq2_pdgId = wq2_pdgId, wq1_pdgId
        
        self.out.fillBranch(prefix + "dr_T", dr_T)
        self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
        self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
        self.out.fillBranch(prefix + "T_Wq_max_pdgId", wq1_pdgId)
        self.out.fillBranch(prefix + "T_Wq_min_pdgId", wq2_pdgId)
        
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
            #genparts = Collection(event, "GenPart")

            # fatjet kinematics
            #self.out.fillBranch(prefix + "is_qualified", fj.is_qualified)
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "mass", getattr(fj, "mass", 0.0))
            #self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "ParT_resonanceMass", fj.ParT_resonanceMass)
            tau1 = getattr(fj, "tau1", -1)
            tau2 = getattr(fj, "tau2", -1)
            tau3 = getattr(fj, "tau3", -1)
            self.out.fillBranch(prefix + "tau21", tau2 / tau1 if tau1 > 0 else 99)
            self.out.fillBranch(prefix + "tau32", tau3 / tau2 if tau2 > 0 else 99)

            # ParticleNet
            self.out.fillBranch(prefix + "PN_Hbb", getattr(fj, "particleNet_prob_Hbb", -1))
            self.out.fillBranch(prefix + "PN_Hcc", getattr(fj, "particleNet_prob_Hcc", -1))
            self.out.fillBranch(prefix + "PN_Hqq", getattr(fj, "particleNet_prob_Hqq", -1))
            self.out.fillBranch(prefix + "PN_QCD", getattr(fj, "particleNet_prob_QCD", -1))
            self.out.fillBranch(prefix + "PN_mass", getattr(fj, "particleNet_mass", -1))

            p_qcd = getattr(fj, "particleNet_prob_QCD", 1.0)
            self.out.fillBranch(prefix + "PN_HccVsQCD", convert_prob(fj, ['Hcc'], None, prefix='particleNet_prob_'))
            self.out.fillBranch(prefix + "PN_HbbVsQCD", convert_prob(fj, ['Hbb'], None, prefix='particleNet_prob_'))
            #self.out.fillBranch(prefix + "PN_HbbVsQCD", getattr(fj, "particleNet_prob_Hbb", 0) / p_qcd if p_qcd>0 else 0)
            #self.out.fillBranch(prefix + "PN_HccVsQCD", getattr(fj, "particleNet_prob_Hcc", 0) / p_qcd if p_qcd>0 else 0)
            self.out.fillBranch(prefix + "PN_HqqVsQCD", convert_prob(fj, ['Hqq'], None, prefix='particleNet_prob_'))

            # ParT main class scores
            self.out.fillBranch(prefix + "scoutGloParT_Xbb", getattr(fj, "scoutGlobalParT_prob_Xbb", -1))
            self.out.fillBranch(prefix + "scoutGloParT_Xcc", getattr(fj, "scoutGlobalParT_prob_Xcc", -1))
            self.out.fillBranch(prefix + "scoutGloParT_Xqq", getattr(fj, "scoutGlobalParT_prob_Xqq", -1))
            self.out.fillBranch(prefix + "scoutGloParT_Xgg", getattr(fj, "scoutGlobalParT_prob_Xgg", -1))

            self.out.fillBranch(prefix + "scoutGloParT_Xbc", getattr(fj, "scoutGlobalParT_prob_Xbc", -1))
            self.out.fillBranch(prefix + "scoutGloParT_Xbs", getattr(fj, "scoutGlobalParT_prob_Xbs", -1))
            self.out.fillBranch(prefix + "scoutGloParT_Xcs", getattr(fj, "scoutGlobalParT_prob_Xcs", -1))
            self.out.fillBranch(prefix + "scoutGloParT_Xss", getattr(fj, "scoutGlobalParT_prob_Xss", -1))

            self.out.fillBranch(prefix + "scoutGloParT_tauhtaue", getattr(fj, "scoutGlobalParT_prob_Xtauhtaue", -1))
            self.out.fillBranch(prefix + "scoutGloParT_tauhtauh", getattr(fj, "scoutGlobalParT_prob_Xtauhtauh", -1))
            self.out.fillBranch(prefix + "scoutGloParT_tauhtaum", getattr(fj, "scoutGlobalParT_prob_Xtauhtaum", -1))

            self.out.fillBranch(prefix + "scoutGloParT_QCD", getattr(fj, "scoutGlobalParT_prob_QCD", -1))

            qcd = getattr(fj, "scoutGlobalParT_prob_QCD", 1.0)
            self.out.fillBranch(prefix + "scoutGloParT_HbbVsQCD", convert_prob(fj, ['Xbb'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_HccVsQCD", convert_prob(fj, ['Xcc'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_HqqVsQCD", convert_prob(fj, ['Xqq'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_twoprong", convert_prob(fj, ['Xbb', 'Xcc', 'Xcs', 'Xqq', 'Xss'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            #self.out.fillBranch(prefix + "scoutGloParT_HqqVsQCD", convert_prob(fj, ['Xqq','Xss'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_HbqVsQCD", convert_prob(fj, ['Xbc','Xbs'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_HcsVsQCD", convert_prob(fj, ['Xcs'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_HbbVsHcc", convert_prob(fj, ['Xbb'], ['Xcc'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_HVsQCD", convert_prob(fj, ['Xbb','Xcc'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_ZVsQCD", convert_prob(fj, ['Xbb','Xcc','Xss','Xqq'], ['QCD'], prefix='scoutGlobalParT_prob_'))
            self.out.fillBranch(prefix + "scoutGloParT_WVsQCD", convert_prob(fj, ['Xcs','Xqq'], ['QCD'], prefix='scoutGlobalParT_prob_'))

            self.out.fillBranch(prefix + "scoutGloParT_massCorrGeneric", getattr(fj, "scoutGlobalParT_massCorrGeneric", -1))
            self.out.fillBranch(prefix + "scoutGloParT_massCorrGenericW2p", getattr(fj, "scoutGlobalParT_massCorrGenericW2p", -1))
            self.out.fillBranch(prefix + "scoutGloParT_massCorrGenericX2p", getattr(fj, "scoutGlobalParT_massCorrGenericX2p", -1))
            self.out.fillBranch(prefix + "scoutGloParT_massCorrResonance", getattr(fj, "scoutGlobalParT_massCorrResonance", -1))
            
            
