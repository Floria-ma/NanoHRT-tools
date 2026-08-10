import logging
import numpy as np
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ..helpers.utils import deltaPhi, polarP4, configLogger
from ..helpers.triggerHelper import passTrigger
from ..helpers.utils import deltaR, closest

from .HeavyFlavBaseProducerScouting import HeavyFlavBaseProducerScouting, METObject

logger = logging.getLogger("nano")
configLogger("nano", loglevel=logging.INFO)

class MuonSampleProducerScouting(HeavyFlavBaseProducerScouting):

    def __init__(self, **kwargs):
        super(MuonSampleProducerScouting, self).__init__(channel="muon", **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(MuonSampleProducerScouting, self).beginFile(
            inputFile, outputFile, inputTree, wrappedOutputTree
        )

        # Trigger information
        self.out.branch("passMuTrig", "O")

        # Delta Phi between jets and MET
        self.out.branch("n_ak4jets", "I")
        self.out.branch("deltaPhi_jets_met", "F", lenVar="n_ak4jets")

        # Selected lepton (muon or electron) matching the ttCR region
        self.out.branch("n_lep", "I")
        self.out.branch("lep_pdgId", "I")

        # Muon variables
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_phi", "F")
        self.out.branch("muon_trk_dxy", "F")
        self.out.branch("muon_trk_dz", "F")
        self.out.branch("muon_trackIso", "F")
        self.out.branch("muon_normchi2", "F")
        self.out.branch("muon_nValidRecoMuonHits", "I")
        self.out.branch("muon_nRecoMuonMatchedStations", "I")
        self.out.branch("muon_nTrackerLayersWithMeasurement", "I")

        # Leptonic W
        self.out.branch("leptonicW_pt", "F")


    def analyze(self, event):

        # print event being analyzed
        # (only for debugging, disable for production)        
        #print(
        #    f">>> run={event.run} lumi={event.luminosityBlock} event={event.event}",
        #    flush=True
        #)

        # muon selection: select events with at least 1 good lepton (muon or electron),
        # matching the vhTreeProducer ttCR region.
        event._allJets = Collection(event, self._ak4_name)
        event._allFatJets = Collection(event, self._fatjet_name)
        event.met = METObject(event, "ScoutingMET")
        event._allMuons = Collection(event, "ScoutingMuonVtx")
        event._allElectrons = Collection(event, "ScoutingElectron")
        event.primaryVertices = Collection(event, "ScoutingPrimaryVertex")

        # muon ID matching the vhTreeProducer ttCR region (lep_id_muon)
        # dxy/dz in scouting are not relative to the PV and are corrected here:
        # https://github.com/cms-sw/cmssw/blob/dc1c73bb9d1735dd1f3a192def64782b4e189f79/DQM/HLTEvF/plugins/ScoutingCollectionMonitor.cc#L85-L103
        def lep_id_muon(mu, pvs):
            # Tight ID
            pass_normchi2 = mu.normchi2 < 3.0  # offline ID is <10 but muon scouting contact recommended 3
            pass_nValidMuonHits = mu.nValidRecoMuonHits > 0
            pass_nRecoMuonMatchedStations = mu.nRecoMuonMatchedStations > 1
            pass_nTrackerLayers = mu.nTrackerLayersWithMeasurement > 5

            pv = pvs[0]
            px = mu.pt * np.cos(mu.phi)
            py = mu.pt * np.sin(mu.phi)
            pz = mu.pt * np.sinh(mu.phi)
            pt2 = mu.pt * mu.pt

            dx = mu.trk_vx - pv.x
            dy = mu.trk_vy - pv.y
            dz = mu.trk_vz - pv.z

            tk_dxyPV = abs((-dx * py + dy * px) / mu.pt)
            tk_dzPV = abs(dz - (dx * px + dy * py) * pz / pt2)

            pass_dxy = tk_dxyPV < 0.2
            pass_dz = tk_dzPV < 0.5

            pass_nValidPixelHits = mu.nValidPixelHits > 0
            pass_nTrackerLayersWithMeasurement = mu.nTrackerLayersWithMeasurement > 5

            pass_noiso = (pass_normchi2 and pass_nValidMuonHits and pass_nRecoMuonMatchedStations and
                          pass_nTrackerLayers and pass_dxy and pass_dz) and \
                         (pass_nValidPixelHits and pass_nTrackerLayersWithMeasurement)

            pass_trackIso = (mu.trackIso / mu.pt) < 0.15

            return pass_noiso and pass_trackIso

        # electron ID matching the vhTreeProducer ttCR region (lep_id_wp80_electron)
        def lep_id_wp80_electron(ele):
            barrel = abs(ele.eta) <= 1.479
            # WP80 ID (no isolation)
            pass_dEta = abs(ele.dEtaIn) < (0.1207 if barrel else 0.1818)
            pass_dPhi = abs(ele.dPhiIn) < (0.0811 if barrel else 0.1221)
            pass_sigmaIEtaIEta = ele.sigmaIetaIeta < (0.0158 if barrel else 0.0289)
            pass_ooEMOop = abs(ele.ooEMOop) < (0.0158 if barrel else 0.0311)
            pass_hOverE = ele.hOverE < (0.0823 if barrel else 0.0986)
            pass_noiso = (pass_dEta and pass_dPhi and pass_sigmaIEtaIEta and
                          pass_ooEMOop and pass_hOverE)
            # WP80 isolation
            pass_trackerIso = (ele.trackIso / ele.pt) < (0.1863 if barrel else 0.0406)
            pass_ecalIso = (ele.ecalIso / ele.pt) < (0.1779 if barrel else 0.2005)
            pass_hcalIso = (ele.hcalIso / ele.pt) < (0.4463 if barrel else 0.1660)
            return pass_noiso and pass_trackerIso and pass_ecalIso and pass_hcalIso

        # build selected leptons (muons + electrons), each tagged with pdgId
        event.leptons = []
        event.leptons_pdgId = []
        event.muons = []
        for mu in event._allMuons:
            if mu.pt > 20 and abs(mu.eta) < 2.4 and lep_id_muon(mu, event.primaryVertices):
                event.muons.append(mu)
                event.leptons.append(mu)
                event.leptons_pdgId.append(13)
        for el in event._allElectrons:
            if el.pt > 20 and abs(el.eta) < 2.5 and lep_id_wp80_electron(el):
                event.leptons.append(el)
                event.leptons_pdgId.append(11)

        # ttCR: require at least one selected lepton (muon or electron)
        if len(event.leptons) == 0: return False
        if len(event.muons) == 0: return False

        # select the leading (highest-pt) lepton; keep its pdgId aligned
        order = sorted(range(len(event.leptons)), key=lambda i: event.leptons[i].pt, reverse=True)
        event.muons = [lep for lep, pid in zip(event.leptons, event.leptons_pdgId) if pid == 13]
        event.leptons = [event.leptons[i] for i in order]
        event.leptons_pdgId = [event.leptons_pdgId[i] for i in order]
        event.mu = event.muons[0] #once sorted by pt, is the leading lepton still muon (? 
        event.lep_pdgId = event.leptons_pdgId[0]

        # select leptons for jet cleaning
        self.selectLeptons(event)

        # get jet & MET collections (no corrections in scouting)
        self.correctJetAndMET(event)

        # ttCR region: no MET or leptonic-W pt selection.
        # leptonic W still computed for the output branch.
        event.mu._mass = 0.1057
        event.leptonicW = polarP4(event.mu, mass="_mass") + event.met.p4()

        # ------------------------------------------------------------------
        # Object cleaning following the vhTreeProducer ttCR _cleanObjects:
        #   - AK8: pt>170, |eta|<2.4, jet_id, ΔR(fj, selectedLeptons) >= 0.8
        #   - AK4: pt>25,  |eta|<2.4, jet_id, ΔR(j, leading AK8) >= 1.2,
        #                                       ΔR(j, selectedLeptons) >= 0.4
        # We build these locally (event.ak8jets_clean, event.ak4jets) and drop
        # the base-class event.ak4jets so the two don't disagree.
        # ------------------------------------------------------------------

        # local jet_id matching vhTreeProducer.jet_id (ScoutingJetID 13.6 TeV).
        # guarded so a missing field doesn't crash if the scouting AK4/AK8
        # collection here doesn't carry the ID inputs.
        def jet_id(jet):
            try:
                return (
                    (abs(jet.eta) < 2.6 and jet.nConstituents > 1 and jet.nCh > 0
                     and jet.neEmEF < 0.9 and jet.muEF < 0.8 and jet.neHEF < 0.99)
                    or (abs(jet.eta) >= 2.6 and abs(jet.eta) < 2.7
                        and jet.nNh > 1 and jet.muEF < 0.8)
                    or (abs(jet.eta) >= 2.7 and abs(jet.eta) < 3.0
                        and jet.nNh > 0 and jet.neEmEF < 0.99)
                )
            except (AttributeError, RuntimeError):
                # ID inputs not available on this collection -> don't reject on ID
                return True

        # ak4 jet truth info 
        # in case genjet is not available (data)
        if self.isMC:
            genjets = Collection(event, "GenJet")
        else:
            genjets = None

        event.fatjets = [fj for fj in event.fatjets  if fj.pt > 170 and abs(fj.eta) < 2.4 and jet_id(fj)]
        event.ak4jets = [j for j in event.ak4jets if jet_id(j)]
        event.ak4jets = sorted(event.ak4jets, key=lambda j: j.pt, reverse=True)
        DeltaPhi_jets_met = []
        for j in event.ak4jets:
            # compute DeltaPhi between jets and MET
            DeltaPhi = abs(deltaPhi(j.phi, event.met.phi))
            DeltaPhi_jets_met.append(DeltaPhi)

            fj, dr = closest(j, event.fatjets)
            j.closestFatJet = fj
            j.closestFatJet_dr = dr

            if self.isMC:
                genjetIdx = getattr(j, "genJetIdx", -1)
                j.genJetIdx = genjetIdx
            else:
                j.genJetIdx = None

            # find the hadron flavour of the closest genjet
            if j.genJetIdx is not None and j.genJetIdx < len(genjets) and j.genJetIdx >= 0:
                genjet = genjets[j.genJetIdx]
                j.flavor = getattr(genjet, "hadronFlavour", -1)
            else:
                j.flavor = -1

        self.out.fillBranch("n_ak4jets", len(event.ak4jets))
        self.out.fillBranch("deltaPhi_jets_met", DeltaPhi_jets_met)

        # fat-jet (probe) selection matching the ttCR region: pt > 170, |eta| < 2.4,
        # and separated from the selected lepton. Require at least one probe jet.
        probe_jets = [fj for fj in event.fatjets
                      if abs(deltaPhi(fj, event.mu)) > 2]
        if len(probe_jets) == 0: return False

        # selection is done, now calculate event properties
        # (in particular of the probe jet)
        probe_jets = probe_jets[:1]
        self.loadGenHistory(event, probe_jets)
        self.evalMassRegression(probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)
        self.out.fillBranch("passMuTrig", passTrigger(event, ["DST_PFScouting_SingleMuon"]))
        self.out.fillBranch("n_lep", len(event.leptons))
        self.out.fillBranch("lep_pdgId", event.lep_pdgId)
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_phi", event.mu.phi)
        self.out.fillBranch("muon_trk_dxy", event.mu.trk_dxy)
        self.out.fillBranch("muon_trk_dz", event.mu.trk_dz)
        self.out.fillBranch("muon_trackIso", event.mu.trackIso)
        self.out.fillBranch("muon_normchi2", event.mu.normchi2)
        self.out.fillBranch("muon_nValidRecoMuonHits", event.mu.nValidRecoMuonHits)
        self.out.fillBranch("muon_nRecoMuonMatchedStations", event.mu.nRecoMuonMatchedStations)
        self.out.fillBranch("muon_nTrackerLayersWithMeasurement", event.mu.nTrackerLayersWithMeasurement)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())

        for mod in self._modules.values():
            if not mod.analyze(event):
                return False

        return True

def MuonTree_2016(): return MuonSampleProducerScouting(year=2016)
def MuonTree_2017(): return MuonSampleProducerScouting(year=2017)
def MuonTree_2018(): return MuonSampleProducerScouting(year=2018)
def MuonTree_2021(): return MuonSampleProducerScouting(year=2021)
def MuonTree_2022(): return MuonSampleProducerScouting(year=2022)
def MuonTree_2023(): return MuonSampleProducerScouting(year=2023)
def MuonTree_2024(): return MuonSampleProducerScouting(year=2024)