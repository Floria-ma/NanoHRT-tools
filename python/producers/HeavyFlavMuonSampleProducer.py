import logging
import numpy as np
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ..helpers.utils import deltaPhi, polarP4, configLogger
from ..helpers.triggerHelper import passTrigger
from ..helpers.utils import deltaR, closest

from .HeavyFlavBaseProducerScouting import HeavyFlavBaseProducerScouting, METObject

logger = logging.getLogger("nano")
configLogger("nano", loglevel=logging.INFO)


# ----------------------------------------------------------------------
# ID / jet_id helpers, copied verbatim from vhTreeProducer so the two
# stay bit-for-bit identical. Keeping them module-level (as in vhTree)
# rather than nested avoids accidental drift.
# ----------------------------------------------------------------------
def lep_id_wp80_electron_noiso(ele):
    barrel = abs(ele.eta) <= 1.479
    pass_dEta = abs(ele.dEtaIn) < (0.1207 if barrel else 0.1818)
    pass_dPhi = abs(ele.dPhiIn) < (0.0811 if barrel else 0.1221)
    pass_sigmaIEtaIEta = ele.sigmaIetaIeta < (0.0158 if barrel else 0.0289)
    pass_ooEMOop = abs(ele.ooEMOop) < (0.0158 if barrel else 0.0311)
    pass_hOverE = ele.hOverE < (0.0823 if barrel else 0.0986)
    return (pass_dEta and pass_dPhi and pass_sigmaIEtaIEta and pass_ooEMOop and pass_hOverE)


def lep_id_wp80_electron(ele):
    barrel = abs(ele.eta) <= 1.479
    pass_trackerIso = (ele.trackIso / ele.pt) < (0.1863 if barrel else 0.0406)
    pass_ecalIso = (ele.ecalIso / ele.pt) < (0.1779 if barrel else 0.2005)
    pass_hcalIso = (ele.hcalIso / ele.pt) < (0.4463 if barrel else 0.1660)
    return (lep_id_wp80_electron_noiso(ele) and pass_trackerIso and pass_ecalIso and pass_hcalIso)


def lep_id_muon_noiso(muo, pvs):
    # Tight ID
    pass_normchi2 = muo.normchi2 < 3.0  # offline ID is <10 but muon scouting contact recommended 3
    pass_nValidMuonHits = muo.nValidRecoMuonHits > 0
    pass_nRecoMuonMatchedStations = muo.nRecoMuonMatchedStations > 1
    pass_nTrackerLayers = muo.nTrackerLayersWithMeasurement > 5

    # dxy/dz in scouting are not relative to the PV and are corrected here:
    # https://github.com/cms-sw/cmssw/blob/dc1c73bb9d1735dd1f3a192def64782b4e189f79/DQM/HLTEvF/plugins/ScoutingCollectionMonitor.cc#L85-L103
    pv = pvs[0]
    px = muo.pt * np.cos(muo.phi)
    py = muo.pt * np.sin(muo.phi)
    pz = muo.pt * np.sinh(muo.phi)
    pt2 = muo.pt * muo.pt

    dx = muo.trk_vx - pv.x
    dy = muo.trk_vy - pv.y
    dz = muo.trk_vz - pv.z

    tk_dxyPV = abs((-dx * py + dy * px) / muo.pt)
    tk_dzPV = abs(dz - (dx * px + dy * py) * pz / pt2)

    pass_dxy = tk_dxyPV < 0.2
    pass_dz = tk_dzPV < 0.5

    pass_nValidPixelHits = muo.nValidPixelHits > 0
    pass_nTrackerLayersWithMeasurement = muo.nTrackerLayersWithMeasurement > 5

    return (pass_normchi2 and pass_nValidMuonHits and pass_nRecoMuonMatchedStations and
            pass_nTrackerLayers and pass_dxy and pass_dz) and (pass_nValidPixelHits and pass_nTrackerLayersWithMeasurement)


def lep_id_muon(muo, pvs):
    pass_trackIso = (muo.trackIso / muo.pt) < 0.15
    return lep_id_muon_noiso(muo, pvs) and pass_trackIso


def jet_id(jet):
    # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/ScoutingJetID13p6TeV
    return (
        (abs(jet.eta) < 2.6 and jet.nConstituents > 1 and jet.nCh > 0
         and jet.neEmEF < 0.9 and jet.muEF < 0.8 and jet.neHEF < 0.99)
        or (abs(jet.eta) >= 2.6 and abs(jet.eta) < 2.7 and jet.nNh > 1 and jet.muEF < 0.8)
        or (abs(jet.eta) >= 2.7 and abs(jet.eta) < 3.0 and jet.nNh > 0 and jet.neEmEF < 0.99)
    )


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

    # ------------------------------------------------------------------
    # Selection helpers mirroring vhTreeProducer's ttCR path:
    #   _selectLeptons -> _preSelect -> _cleanObjects -> _selectEvent
    # ------------------------------------------------------------------

    def _selectLeptons(self, event):
        # vhTreeProducer._selectLeptons: build looseLeptons (electrons + muons),
        # used for lepton counting and AK8 preselection cleaning.
        event.looseLeptons = []

        electrons = Collection(event, "ScoutingElectron")
        for el in electrons:
            el.etaSC = el.eta
            if el.pt > 20 and abs(el.eta) < 2.5:
                event.looseLeptons.append(el)

        muons = Collection(event, "ScoutingMuonVtx")
        for mu in muons:
            if mu.pt > 20 and abs(mu.eta) < 2.4:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

    def _preSelect(self, event):
        # vhTreeProducer._preSelect (ttCR branch only).
        event.primaryVertices = Collection(event, "ScoutingPrimaryVertex")

        event.selectedLeptons = []       # used to reconstruct the vector boson
        event.selectedLeptons_pdgId = []
        event.nonIsoLeptons = []
        event.nonIsoLeptons_pdgId = []

        # loose AK8: |eta|<2.4, dR(fj, looseLeptons) >= 0.8, jet_id
        ak8tmp = [fj for fj in event._allAK8jets
                  if (abs(fj.eta) < 2.4 and closest(fj, event.looseLeptons)[1] >= 0.8) and jet_id(fj)]
        if len(ak8tmp) == 0:
            return False

        for lep in event.looseLeptons:
            try:
                isElectron = (lep.sigmaIetaIeta >= 0)
            except Exception:
                isElectron = False

            if not isElectron:
                if lep_id_muon(lep, event.primaryVertices):
                    event.selectedLeptons_pdgId.append(13)
                    event.selectedLeptons.append(lep)
                elif lep_id_muon_noiso(lep, event.primaryVertices):
                    event.nonIsoLeptons_pdgId.append(13)
                    event.nonIsoLeptons.append(lep)
            else:
                if lep_id_wp80_electron(lep):
                    event.selectedLeptons_pdgId.append(11)
                    event.selectedLeptons.append(lep)
                elif lep_id_wp80_electron_noiso(lep):
                    event.nonIsoLeptons_pdgId.append(11)
                    event.nonIsoLeptons.append(lep)

        # ttCR: require at least one selected (isolated) lepton
        if len(event.selectedLeptons) == 0:
            return False

        # at least one of the ak8 candidates must pass a loose tagger bb/cc score cut
        self.evalTagger(ak8tmp)
        if self._opts['tagger_threshold_cc'] is not None:
            ak8tmp = [fj for fj in ak8tmp if (
                (fj.pn_HccVsQCD >= self._opts['tagger_threshold_cc']) or (fj.pn_HbbVsQCD >= self._opts['tagger_threshold_bb']))]
            if len(ak8tmp) == 0:
                return False

        return True

    def _cleanObjects(self, event):
        # vhTreeProducer._cleanObjects — AK4 leg only.
        # event.ak8jets is set by analyze() from the DECORATED event.fatjets
        # probe BEFORE this is called, so we do not rebuild it here (rebuilding
        # from raw event._allAK8jets would give undecorated jets that crash
        # fillFatJetInfo on data). We only build the cleaned AK4 collection,
        # referencing event.ak8jets[:1] for the ΔR>1.2 cut.
        event.central_jets = []
        event.ak4jets = []

        for j in event._allJets:
            if not (j.pt > 25 and abs(j.eta) < 2.4 and jet_id(j)):
                continue
            event.central_jets.append(j)
            # ttCR: require ΔR(jet, leading AK8) > 1.2
            if closest(j, event.ak8jets[:1])[1] < 1.2:
                continue
            if closest(j, event.selectedLeptons)[1] < 0.4:
                continue
            event.ak4jets.append(j)

    def _selectEvent(self, event):
        # vhTreeProducer._selectEvent (ttCR branch).
        if len(event.ak8jets) < 1:
            return False

        event.ak8 = event.ak8jets[0]

        # ttCR-specific bookkeeping (lepton-jet ΔR), matching vhTree.
        from ..helpers.utils import minValue
        event._min_dr_lep1_jet = minValue([deltaR(event.selectedLeptons[0], j) for j in event.ak4jets])
        if len(event.selectedLeptons) > 1:
            event._min_dr_lep2_jet = minValue([deltaR(event.selectedLeptons[1], j) for j in event.ak4jets])
        else:
            event._min_dr_lep2_jet = -99.
        if len(event.nonIsoLeptons) > 0:
            event._min_dr_lepnoniso_jet = minValue([deltaR(event.nonIsoLeptons[0], j) for j in event.ak4jets])
        else:
            event._min_dr_lepnoniso_jet = -99.

        return True

    def analyze(self, event):

        # collections
        event._allJets = Collection(event, self._ak4_name)
        event._allFatJets = Collection(event, self._fatjet_name)
        # aliases used by the vhTree-style helpers
        event._allAK8jets = event._allFatJets
        event.met = METObject(event, "ScoutingMET")

        # ---- vhTree ttCR flow: _selectLeptons -> _preSelect -> correct -> _cleanObjects -> _selectEvent
        self._selectLeptons(event)
        if self._preSelect(event) is False:
            return False

        # JEC/MET corrections. Base class needs event._allMuons / _allElectrons
        # for its own selectLeptons; provide them, but note that the *cleaning*
        # below is done by our own _cleanObjects, not the base class.
        event._allMuons = Collection(event, "ScoutingMuonVtx")
        event._allElectrons = Collection(event, "ScoutingElectron")
        self.correctJetAndMET(event)

        # After correctJetAndMET the fatjets/jets are pt-sorted in place on
        # event._allFatJets / event._allJets; re-alias so _cleanObjects sees
        # the corrected, sorted collections.
        event._allAK8jets = event._allFatJets

        # ---- ttCR probe jet from the DECORATED base-class fatjets.
        # correctJetAndMET populates event.fatjets with jets that carry .idx,
        # .is_qualified, and (via the msoftdrop selection cut) a read msoftdrop.
        # We must feed THESE to evalMassRegression/loadGenHistory/fillFatJetInfo,
        # NOT the raw event.ak8jets built in _cleanObjects — those are undecorated
        # and lack msoftdrop/ParT_resonanceMass, which fillFatJetInfo reads
        # unguarded and would null-pointer on data.
        # Apply the ttCR probe requirement (pt>170, |eta|<2.4, jet_id, dR>=0.8
        # vs the selected leptons) on top of the base-class fatjet selection.
        probe_jets = [fj for fj in event.fatjets
                      if fj.pt > 170 and abs(fj.eta) < 2.4 and jet_id(fj)
                      and closest(fj, event.selectedLeptons)[1] >= 0.8]
        if len(probe_jets) == 0:
            return False
        probe_jets = probe_jets[:1]

        # expose the decorated probe as event.ak8jets so _cleanObjects' AK4
        # ΔR(jet, leading AK8) >= 1.2 cut references the same probe jet.
        event.ak8jets = probe_jets

        self._cleanObjects(event)   # AK4-only now (see method); leaves event.ak8jets intact

        # leading selected lepton (muon OR electron), exactly as vhTree
        event.lep = event.selectedLeptons[0]
        event.lep_pdgId = event.selectedLeptons_pdgId[0]

        # leptonic W (built from the leading selected lepton, whatever flavour).
        # mass assigned by pdgId so an electron event isn't given the muon mass.
        event.lep._mass = 0.000511 if abs(event.lep_pdgId) == 11 else 0.1057
        event.leptonicW = polarP4(event.lep, mass="_mass") + event.met.p4()

        # ---- deltaPhi(jets, MET) and flavour bookkeeping over the CLEANED ak4jets
        DeltaPhi_jets_met = []
        if self.isMC:
            genjets = Collection(event, "GenJet")
        else:
            genjets = None

        for j in event.ak4jets:
            DeltaPhi_jets_met.append(abs(deltaPhi(j.phi, event.met.phi)))

            fj, dr = closest(j, event.ak8jets)
            j.closestFatJet = fj
            j.closestFatJet_dr = dr

            if self.isMC:
                genjetIdx = getattr(j, "genJetIdx", -1)
                j.genJetIdx = genjetIdx
            else:
                j.genJetIdx = None

            if j.genJetIdx is not None and 0 <= j.genJetIdx < len(genjets):
                genjet = genjets[j.genJetIdx]
                j.flavor = getattr(genjet, "hadronFlavour", -1)
            else:
                j.flavor = -1

        # ---- fatjet info on the DECORATED probe (event.ak8jets == probe_jets).
        # evalMassRegression sets ParT_resonanceMass/sdmass BEFORE fillFatJetInfo
        # reads them, exactly as in the reference muon producer.
        self.evalMassRegression(event.ak8jets)
        if self.isMC:
            self.loadGenHistory(event, event.ak8jets)
        self.fillFatJetInfo(event, event.ak8jets)

        # ---- fill output
        self.fillBaseEventInfo(event)

        self.out.fillBranch("passMuTrig", passTrigger(event, ["DST_PFScouting_SingleMuon"]))
        self.out.fillBranch("n_ak4jets", len(event.ak4jets))
        self.out.fillBranch("deltaPhi_jets_met", DeltaPhi_jets_met)
        self.out.fillBranch("n_lep", len(event.selectedLeptons))
        self.out.fillBranch("lep_pdgId", event.lep_pdgId)

        # kinematics exist for both muons and electrons
        self.out.fillBranch("muon_pt", event.lep.pt)
        self.out.fillBranch("muon_eta", event.lep.eta)
        self.out.fillBranch("muon_phi", event.lep.phi)

        # muon-specific branches: only valid when the leading selected lepton
        # is a muon. The datamodel raises RuntimeError (not AttributeError) on a
        # missing branch, so getattr(..., default) does NOT protect us here —
        # an electron lep would throw on ScoutingElectron_trk_dxy etc.
        if abs(event.lep_pdgId) == 13:
            mu = event.lep
            self.out.fillBranch("muon_trk_dxy", mu.trk_dxy)
            self.out.fillBranch("muon_trk_dz", mu.trk_dz)
            self.out.fillBranch("muon_trackIso", mu.trackIso)
            self.out.fillBranch("muon_normchi2", mu.normchi2)
            self.out.fillBranch("muon_nValidRecoMuonHits", mu.nValidRecoMuonHits)
            self.out.fillBranch("muon_nRecoMuonMatchedStations", mu.nRecoMuonMatchedStations)
            self.out.fillBranch("muon_nTrackerLayersWithMeasurement", mu.nTrackerLayersWithMeasurement)
        else:
            self.out.fillBranch("muon_trk_dxy", -99.)
            self.out.fillBranch("muon_trk_dz", -99.)
            self.out.fillBranch("muon_trackIso", -99.)
            self.out.fillBranch("muon_normchi2", -99.)
            self.out.fillBranch("muon_nValidRecoMuonHits", -99)
            self.out.fillBranch("muon_nRecoMuonMatchedStations", -99)
            self.out.fillBranch("muon_nTrackerLayersWithMeasurement", -99)
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