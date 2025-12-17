import logging
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ..helpers.utils import deltaPhi, polarP4, configLogger
from ..helpers.triggerHelper import passTrigger

from .HeavyFlavBaseProducerScouting import HeavyFlavBaseProducerScouting

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

        # Muon variables
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_phi", "F")
        #self.out.branch("muon_miniIso", "F")

        # Leptonic W
        self.out.branch("leptonicW_pt", "F")

    def analyze(self, event):

        # print event being analyzed
        # (only for debugging, disable for production)        
        #print(
        #    f">>> run={event.run} lumi={event.luminosityBlock} event={event.event}",
        #    flush=True
        #)

        # muon selection: select events with exactly 1 good muon, reject all others
        event._allMuons = Collection(event, "ScoutingMuonVtx")
        event.muons = [mu for mu in event._allMuons if mu.pt > 55 and abs(mu.eta) < 2.4]
        if len(event.muons) != 1: return False
        event.mu = event.muons[0]

        # select leptons for jet cleaning
        self.selectLeptons(event)

        # get jet & MET collections (no corrections in scouting)
        self.correctJetAndMET(event)

        # MET selection
        if event.met.pt < 50: return False

        # leptonic W pt selection
        event.mu._mass = 0.1057
        event.leptonicW = polarP4(event.mu, mass="_mass") + event.met.p4()
        if event.leptonicW.Pt() < 100: return False

        # b-jet selection: select events where there is at least one b-tagged jet
        # relatively close to the selected muon.
        # note: to update for scouting, right now just use all jets.
        bjets = [j for j in event.ak4jets if abs(deltaPhi(j, event.mu)) < 2]
        if len(bjets) == 0: return False
        event.bjets = bjets

        # fat-jet selection: select events where there is at least one fat jet
        # relatively far from the selected muon.
        if len(event.fatjets) == 0: return False
        probe_jets = [fj for fj in event.fatjets if abs(deltaPhi(fj, event.mu)) > 2]
        if len(probe_jets) == 0: return False

        # selection is done, now calculate event properties
        probe_jets = probe_jets[:1]
        self.loadGenHistory(event, probe_jets)
        self.evalMassRegression(probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)
        self.out.fillBranch("passMuTrig", passTrigger(event, ["DST_PFScouting_SingleMuon"]))
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_phi", event.mu.phi)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())

        return True

def MuonTree_2016(): return MuonSampleProducerScouting(year=2016)
def MuonTree_2017(): return MuonSampleProducerScouting(year=2017)
def MuonTree_2018(): return MuonSampleProducerScouting(year=2018)
def MuonTree_2021(): return MuonSampleProducerScouting(year=2021)
def MuonTree_2022(): return MuonSampleProducerScouting(year=2022)
def MuonTree_2023(): return MuonSampleProducerScouting(year=2023)
def MuonTree_2024(): return MuonSampleProducerScouting(year=2024)
