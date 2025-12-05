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

lumi_dict = {2015: 19.52, 2016: 16.81, 2017: 41.48, 2018: 59.83, 
             2021: 7.98, 2022: 26.67, 2023:17.794, 2024:9.451}

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

        self._channel = channel    # 'qcd', 'photon', 'inclusive', 'muon'
        self.year = int(kwargs["year"])
        self.jetType = kwargs.get("jetType", "scouting").lower()
        self._jmeSysts = {'jec': False, 'jes': None, 'jes_source': '', 'jes_uncertainty_file_prefix': '',
                          'jer': None, 'jmr': None, 'met_unclustered': None, 'smearMET': True, 'applyHEMUnc': False,
                          'jesr_extra_br': True}
        
        self._opts = {
            "WRITE_CACHE_FILE": False,
            'mass_range': (40, 250),
        }

        for k in kwargs:
            if k in self._jmeSysts:
                self._jmeSysts[k] = kwargs[k]
            else:
                self._opts[k] = kwargs[k]

        #self._needsJMECorr = False  
        #self._needsJMECorr = any([self._jmeSysts['jec'], self._jmeSysts['jes'],
        #                          self._jmeSysts['jer'], self._jmeSysts['jmr'],
        #                          self._jmeSysts['met_unclustered'], self._jmeSysts['applyHEMUnc']])
        #logger.info('Running %s channel for year %s with JME systematics %s, other options %s',
        #            self._channel, str(self._year), str(self._jmeSysts), str(self._opts))


        self._doJetCleaning = True

        # pfjet collection
        self._ak4_name = "ScoutingPFJet"

        # fatjet collection
        self._fatjet_name = "ScoutingFatPFJetRecluster"
        self._jetConeSize = 0.8  # ak8

        logger.info(
            "HeavyFlavBaseProducerScouting: channel=%s, year=%s, ak4=%s, fatjets=%s",
            self._channel, self.year, self._ak4_name, self._fatjet_name
        )

        # DeepJet WPs 
        self.DeepJet_WP_L = {2015: 0.0508, 2016: 0.0480, 2017: 0.0532, 2018: 0.0490, 2021:0.0583, 2022: 0.0614, 2023: 0.0479, 2024: 0.048}[self.year]
        self.DeepJet_WP_M = {2015: 0.2598, 2016: 0.2489, 2017: 0.3040, 2018: 0.2783, 2021:0.3086, 2022: 0.3196, 2023: 0.2431, 2024: 0.2435}[self.year]
        self.DeepJet_WP_T = {2015: 0.6502, 2016: 0.6377, 2017: 0.7476, 2018: 0.7100, 2021:0.7183, 2022: 0.7300, 2023: 0.6553, 2024: 0.6563}[self.year]

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
        # nothing heavy to initialize for scouting
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # MC or data
        self.isMC = bool(inputTree.GetBranch("genWeight"))

        # try to detect an available rho-like branch
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

        for idx in ([1, 2] if self._channel in ['qcd', 'mutagged'] else [1]):
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            #self.out.branch(prefix + "is_qualified", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "mass", "F")
            #self.out.branch(prefix + "rawmass", "F")
            #self.out.branch(prefix + "sdmass", "F")
            #self.out.branch(prefix + "regressed_mass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")
            #self.out.branch(prefix + "btagcsvv2", "F")
            #self.out.branch(prefix + "btagjp", "F")

            #scoutGlobalParT
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

            # VS-QCD scores (auto-computed)
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

            #self.out.branch(prefix + "scoutGloParT_HbbVsQCD", "F")
            #self.out.branch(prefix + "scoutGloParT_HccVsQCD", "F")

            # massCorr branches
            self.out.branch(prefix + "scoutGloParT_massCorrGeneric", "F")
            self.out.branch(prefix + "scoutGloParT_massCorrGenericW2p", "F")
            self.out.branch(prefix + "scoutGloParT_massCorrGenericX2p", "F")
            self.out.branch(prefix + "scoutGloParT_massCorrResonance", "F")

            # matching variables
            if self.isMC:
                # info of the closest hadGenTop
                #self.out.branch(prefix + "dr_T", "F")
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


            if self.year <= 2017:
                self.out.branch("l1PreFiringWeight", "F")
                self.out.branch("l1PreFiringWeightUp", "F")
                self.out.branch("l1PreFiringWeightDown", "F")
            else:
                self.out.branch("l1PreFiringWeight", "F")
                self.out.branch("l1PreFiringWeightUp", "F")
                self.out.branch("l1PreFiringWeightDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # nothing to clean up
        pass

    def selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []  # used for jet lepton cleaning & lepton counting

        electrons = Collection(event, "ScoutingElectron")
        for el in electrons:
            el.etaSC = el.eta #+ el.deltaEtaSC
            if el.pt > 20 and abs(el.eta) < 2.5:# and el.mvaIso_WP90:
                event.looseLeptons.append(el)

        muons = Collection(event, "ScoutingMuonVtx")
        for mu in muons:
            if mu.pt > 20 and abs(mu.eta) < 2.4:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)
    
    def correctJetAndMET(self, event):
       
        event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)

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

            event.met = _DummyMET()

        # sort jets by pt
        try:
            event._allJets = sorted(event._allJets, key=lambda x: x.pt, reverse=True)
        except Exception:
            pass

        try:
            for idx, fj in enumerate(event._allFatJets):
                fj.idx = idx
                fj.is_qualified = True
            event._allFatJets = sorted(event._allFatJets, key=lambda x: x.pt, reverse=True)
        except Exception:
            pass
        
        # select _allJets jets
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

        # select fatjets (lepton-cleaned)
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

    def _evalMassRegression(self, jets):
        for j in jets:
            j.masses = {
                'sdmass': j.msoftdrop,
                #'regressed_mass': j.particleNet_mass,
                'ParT_resonanceMass': j.scoutGlobalParT_massCorrResonance * j.mass,# * (1. - j.rawFactor),
                #'ParT_visibleMass': j.ParT_visiableMassCorr * j.mass * (1. - j.rawFactor),
            }

    def fillBaseEventInfo(self, event):
        """
        Fill basic event-level info: year, lumi weight, MET filters, nlep, ht, MET.
        """
        # jet radius (use reclustered fatjet R)
        self.out.fillBranch("jetR", self._jetConeSize)
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
        
        #self.out.fillBranch("met", event.met.pt)
        #self.out.fillBranch("met_phi", event.met.phi)
        #if self._channel == '0L' or self._channel == 'had':
        #self.out.fillBranch("dphi_met_tkmet", event._dphi_met_tkmet)
        #self.out.fillBranch("min_dphi_met_jet", event._min_dphi_met_jet)

        self.out.fillBranch("passmetfilters", met_filters)
        # trigger
        #self.out.fillBranch("passTrigPFHTScouting", passTrigger(event,["DST_PFScouting_JetHT"]))

        # L1 prefire weights: not really defined for scouting; default to 1 if missing
        if self.year <= 2017:
            self.out.fillBranch(
                "l1PreFiringWeight",
                getattr(event, "L1PreFiringWeight_Nom", 1.0),
            )
            self.out.fillBranch(
                "l1PreFiringWeightUp",
                getattr(event, "L1PreFiringWeight_Up", 1.0),
            )
            self.out.fillBranch(
                "l1PreFiringWeightDown",
                getattr(event, "L1PreFiringWeight_Dn", 1.0),
            )
        else:
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

    '''
    def _selectEvent(self, event, fatjets):
            # logger.debug('processing event %d' % event.event)
            #event.Vboson = None
            #event.higgs = None
            #event.tag = None  
            #event.ak8jets = fatjets

            # sort by pt
            #event.ak8jets = sorted(fatjets, key=lambda x: x.pt, reverse=True)
            fatjets = sorted(fatjets, key=lambda x: x.pt, reverse=True)
            #event.ak8jets = sorted(event.ak8jets, key=lambda x: x.pt, reverse=True)
            if fatjets is None or len(fatjets) < 2:
                return False

            mass0 = fatjets[0].scoutGlobalParT_massCorrResonance * fatjets[0].mass
            mass1 = fatjets[1].scoutGlobalParT_massCorrResonance * fatjets[1].mass

            if (mass0 < mass1):
                H_idx=1
                V_idx=0
            else:
                H_idx=0
                V_idx=1

            # #these are the actual jet objects, mass ordered, but this assignment is not used in the final analysis
            event.recoV = fatjets[V_idx]
            event.recoH = fatjets[H_idx]
            #event.H_idx = H_idx

            #these have just the p4
            event.Vboson = fatjets[V_idx].p4()
            event.higgs = fatjets[H_idx].p4()

            #these are the jet objects pt ordered
            if len(fatjets)<2:
                return False
            event.ak8 = fatjets[0]
            event.ak8_2 = fatjets[1]

            if abs(deltaPhi(event.ak8.phi, event.ak8_2.phi)) < 2.2:
                return False

            #if abs(deltaPhi(event.ak8.p4().Phi(), event.ak8_2.p4().Phi())) < 2.2:
            #    return False
        
            self._evalMassRegression([event.ak8_2])
            
          
            #if self._opts['mass_range'] is not None:
                # if any mass is in the window, keep the event
            #    inMassWindow = False
            #    for mass in event.ak8_2.masses.values():
            #        if (self._opts['mass_range'][0] <= mass <= self._opts['mass_range'][1]):
            #                inMassWindow = True
            #                break
            #    if not inMassWindow:
         
            # channel specific selections
            #if self._channel == 'had':
            ht = sum([j.pt for j in getattr(event, "_allJets", [])])
            #ht = sum([j.pt for j in event._allJets])
            if (event.met.pt > 150) :
                    return False
                
            # return True if passes selection
            return True
    '''

    def loadGenHistory(self, event, fatjets):
        if not self.isMC:
            return

        # use cleaned fatjets if available
        if fatjets is None or len(fatjets) < 1:
            return False

        ak8  = fatjets[0]
        ak8_2 = fatjets[1] if len(fatjets) > 1 else None

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

        genW1, genWdr1 = closest(ak8, hadGenWs)
        if ak8_2 is not None:
            genW2, genWdr2 = closest(ak8_2, hadGenWs)
        else:
            genW2, genWdr2 = (None, 999)

        genW = genW1 if genWdr1<genWdr2 else genW2

        # info of the closest hadGenW
        wdecay_ = max([abs(d.pdgId) for d in get_daughters(genW)], default=0) if genW else 0
        #self.out.fillBranch(prefix + "dr_W", fj.dr_W)
        self.out.fillBranch(prefix + "dr_W_daus",
                            max([deltaR(ak8, dau) for dau in get_daughters(genW)]) if genW else 99)
        self.out.fillBranch(prefix + "W_decay", wdecay_)

        # sort tops by deltaR to ak8
        hadGenTops.sort(key=lambda x: deltaR2(x, ak8))
        t = hadGenTops[0]
        self.out.fillBranch(prefix + "dr_T_b", deltaR(ak8, t.genB) if len(hadGenTops) else 99)
        dW = get_daughters(t.genW) if len(hadGenTops) else []
        print(dW)
        drwq1, drwq2 = [deltaR(ak8, dau) for dau in get_daughters(
                hadGenTops[0].genW)] if len(hadGenTops) else [99, 99]
        wq1_pdgId, wq2_pdgId = [dau.pdgId for dau in get_daughters(hadGenTops[0].genW)] if len(hadGenTops) else [0, 0]
        if drwq1 < drwq2:
                drwq1, drwq2 = drwq2, drwq1
                wq1_pdgId, wq2_pdgId = wq2_pdgId, wq1_pdgId
        self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
        self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
        self.out.fillBranch(prefix + "T_Wq_max_pdgId", wq1_pdgId)
        self.out.fillBranch(prefix + "T_Wq_min_pdgId", wq2_pdgId)
        
        '''
        lepGenTops.sort(key=lambda x: deltaR2(x, event.ak8))  # sort by deltaR
        self.out.fillBranch("dr_ak8_leptop1_b", deltaR(event.ak8, lepGenTops[0].genB) if len(lepGenTops) > 0 else 99)
        wlep = [dau for dau in get_daughters(lepGenTops[0].genW) if abs(
                dau.pdgId) in (11, 13, 15)][0] if len(lepGenTops) > 0 else None
        self.out.fillBranch("dr_ak8_leptop1_wlep", deltaR(event.ak8, wlep) if wlep else 99)
        self.out.fillBranch("leptop1_wlep_pdgId", wlep.pdgId if wlep else 0)

        self.out.fillBranch("dr_ak8_leptop2_b", deltaR(event.ak8, lepGenTops[1].genB) if len(lepGenTops) > 1 else 99)
        wlep = [dau for dau in get_daughters(lepGenTops[1].genW) if abs(
                dau.pdgId) in (11, 13, 15)][0] if len(lepGenTops) > 1 else None
        self.out.fillBranch("dr_ak8_leptop2_wlep", deltaR(event.ak8, wlep) if wlep else 99)
        self.out.fillBranch("leptop2_wlep_pdgId", wlep.pdgId if wlep else 0)
        '''
        '''
        if len(dW) >= 2:
                q1, q2 = dW[0], dW[1]
                drwq1, drwq2 = deltaR(ak8, q1), deltaR(ak8, q2)
                if drwq1 < drwq2:
                    drwq1, drwq2 = drwq2, drwq1
                    q1, q2 = q2, q1
                self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
                self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
                self.out.fillBranch(prefix + "T_Wq_max_pdgId", q1.pdgId)
                self.out.fillBranch(prefix + "T_Wq_min_pdgId", q2.pdgId)
        else:
                self.out.fillBranch(prefix + "dr_T_Wq_max", 99)
                self.out.fillBranch(prefix + "dr_T_Wq_min", 99)
                self.out.fillBranch(prefix + "T_Wq_max_pdgId", 0)
                self.out.fillBranch(prefix + "T_Wq_min_pdgId", 0)
        '''

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
            #self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            #self.out.fillBranch(prefix + "regressed_mass", fj.regressed_mass)
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
            
            '''
            if self.isMC:
                fj = fatjets[idx - 1]
                #genparts = Collection(event, "GenPart")

                
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
                
                # find b closest to fj
                drs_b = [(deltaR(fj, gp), gp) for gp in genparts if abs(gp.pdgId) == 5]
                if drs_b:
                    dr_b, genb = min(drs_b, key=lambda x: x[0])
                    self.out.fillBranch(prefix + "dr_T_b", dr_b)
                else:
                    self.out.fillBranch(prefix + "dr_T_b", 99)

                # find W daughters
                Wd = [gp for gp in genparts if abs(gp.pdgId) in [1,2,3,4]]   # light quarks
                drs_Wq = [(deltaR(fj, gp), gp) for gp in Wd]
                if drs_Wq:
                    dr_Wq_max, gpmax = min(drs_Wq, key=lambda x: x[0])
                    self.out.fillBranch(prefix + "dr_T_Wq_max", dr_Wq_max)
                    self.out.fillBranch(prefix + "T_Wq_max_pdgId", gpmax.pdgId)
                else:
                    self.out.fillBranch(prefix + "dr_T_Wq_max", 99)
                    self.out.fillBranch(prefix + "T_Wq_max_pdgId", 0)

                # ΔR(W daughters)
                if len(Wd) >= 2:
                    # choose the first two daughters arbitrarily (you can improve logic)
                    dr_W_daus = deltaR(Wd[0], Wd[1])
                    self.out.fillBranch(prefix + "dr_W_daus", dr_W_daus)
                else:
                    self.out.fillBranch(prefix + "dr_W_daus", 99)

                self.out.fillBranch(prefix + "dr_T", fj.dr_T)
                self.out.fillBranch(prefix + "dr_T_b", deltaR(fj, fj.genT.genB) if fj.genT else 99)
                self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
                self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
                self.out.fillBranch(prefix + "T_Wq_max_pdgId", wq1_pdgId)
                self.out.fillBranch(prefix + "T_Wq_min_pdgId", wq2_pdgId)
                self.out.fillBranch(prefix + "T_pt", fj.genT.pt if fj.genT else -1)
            '''
        

