import logging
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ..helpers.utils import deltaPhi, polarP4, configLogger
from ..helpers.triggerHelper import passTrigger
from ..helpers.utils import deltaR, closest
from ..helpers.nnHelper import convert_prob

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
        self.out.branch("muon_trk_dxy", "F")
        self.out.branch("muon_trk_dz", "F")
        self.out.branch("muon_trackIso", "F")

        # Leptonic W
        self.out.branch("leptonicW_pt", "F")

        # B-tagged jet close to muon
        self.out.branch("bjet_pt", "F")
        self.out.branch("bjet_eta", "F")
        self.out.branch("bjet_phi", "F")
        self.out.branch("bjet_particleNet_prob_b", "F")
        self.out.branch("bjet_particleNet_prob_bb", "F")
        self.out.branch("bjet_particleNet_prob_c", "F")
        self.out.branch("bjet_particleNet_prob_cc", "F")
        self.out.branch("bjet_particleNet_prob_uds", "F")
        self.out.branch("bjet_particleNet_prob_g", "F")

        # closest fat jet to the b-jet scores
        prefix = "bjet_closestFatJet_"
        self.out.branch("n_bjets", "I")
        self.out.branch(prefix + "scoutGloParT_Xbb", "F", lenVar="n_bjets")
        self.out.branch(prefix + "scoutGloParT_Xbc", "F", lenVar="n_bjets")
        self.out.branch(prefix + "scoutGloParT_Xbs", "F", lenVar="n_bjets")
        self.out.branch(prefix + "scoutGloParT_HbbVsQCD", "F", lenVar="n_bjets")
        self.out.branch(prefix + "scoutGloParT_HbqVsQCD",  "F", lenVar="n_bjets")
        self.out.branch(prefix + "scoutGloParT_HcsVsQCD",  "F", lenVar="n_bjets")
        self.out.branch(prefix + "scoutGloParT_HbbVsHcc", "F")
        self.out.branch(prefix + "hadronFlavour",  "I", lenVar="n_bjets")
        self.out.branch(prefix + "ScoutingPFJetRecluster_genJetIdx", "I", lenVar="n_bjets")
        self.out.branch(prefix + "ak4_ak8_dR", "F")

    def analyze(self, event):

        # print event being analyzed
        # (only for debugging, disable for production)        
        #print(
        #    f">>> run={event.run} lumi={event.luminosityBlock} event={event.event}",
        #    flush=True
        #)

        # muon selection: select events with exactly 1 good muon, reject all others
        event._allMuons = Collection(event, "ScoutingMuonVtx")
        event.muons = [mu for mu in event._allMuons
                       if mu.pt > 55 and abs(mu.eta) < 2.4
                       and abs(mu.trk_dxy) < 0.15
                       #and abs(mu.trk_dz) < 1.0
                       and mu.trackIso < 0.1
                      ]
        if len(event.muons) != 1: return False

        # select the (unique) muon found above
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

        # for each ak4 jet, find the closest ak8 jet from the uncleaned collection
        # some baseline cuts on the fatjets
        allFatJetsSelected = [fj for fj in event._allFatJets
                              if fj.pt > 200 and abs(fj.eta) < 2.4]

        # Fatjet truth info 
        # enable branches explicitly
        #event._tree.SetBranchStatus("nGenJet", 1)
        #event._tree.SetBranchStatus("GenJet_*", 1)
        genjets = Collection(event, "GenJet") 

        for j in event.ak4jets:
            prefix = "bjet_closestFatJet_"
            fj, dr = closest(j, allFatJetsSelected)
            j.closestFatJet = fj
            j.closestFatJet_dr = dr

            genjetIdx = getattr(j, "genJetIdx", -1)
            # find the hadron flavour of the closest genjet
            if genjetIdx is not None and genjetIdx < len(genjets) and genjetIdx >= 0:
                genjet = genjets[genjetIdx]
                j.flavor = getattr(genjet, "hadronFlavour", -1)
            else:
                j.flavor = -1

            #if fj and dr < 0.8:
                # define the scouting ParT b tagging scores for the closest fatjet
            j.closestak8_scoutGlobalParT_prob_Xbb = fj.scoutGlobalParT_prob_Xbb if fj else -1
            j.closestak8_scoutGlobalParT_prob_Xbc = fj.scoutGlobalParT_prob_Xbc if fj else -1
            j.closestak8_scoutGlobalParT_prob_Xbs = fj.scoutGlobalParT_prob_Xbs if fj else -1
            j.closestak8_scoutGlobalParT_HbbVsQCD = convert_prob(fj, ['Xbb'], ['QCD'], prefix='scoutGlobalParT_prob_') if fj else -1
            j.closestak8_scoutGlobalParT_HbqVsQCD = convert_prob(fj, ['Xbc','Xbs'], ['QCD'], prefix='scoutGlobalParT_prob_') if fj else -1
            j.closestak8_scoutGlobalParT_HcsVsQCD = convert_prob(fj, ['Xcs'], ['QCD'], prefix='scoutGlobalParT_prob_') if fj else -1
            j.closestak8_scoutGlobalParT_HbbVsHcc = convert_prob(fj, ['Xbb'], ['Xcc'], prefix='scoutGlobalParT_prob_') if fj else -1

        # b-jet selection: select events where there is at least one b-tagged jet
        # relatively close to the selected muon.
        # note: preliminary; to find out which scores and threshold to use for 2024 scouting.
        bjets = [j for j in event.ak4jets if abs(deltaPhi(j, event.mu)) < 2
                 and j.closestFatJet is not None
                 and j.closestFatJet_dr < 0.8
                 #and j.closestak8_scoutGlobalParT_HbqVsQCD > self.scouting_ak4_PNet_WP_M 
                ]
        if len(bjets) == 0: return False

        HbqvsQCD = []
        flavors = []
        Xbb = []
        Xbc = []
        Xbs = []
        HbbvsQCD = []
        HcsvsQCD = []
        genidx = []

        # fill bjet closest fatjet scores
        prefix = "bjet_closestFatJet_"  
        for j in bjets:
            HbqvsQCD.append(j.closestak8_scoutGlobalParT_HbqVsQCD)
            Xbb.append(j.closestak8_scoutGlobalParT_prob_Xbb)
            Xbc.append(j.closestak8_scoutGlobalParT_prob_Xbc)
            Xbs.append(j.closestak8_scoutGlobalParT_prob_Xbs)
            HbbvsQCD.append(j.closestak8_scoutGlobalParT_HbbVsQCD)
            HcsvsQCD.append(j.closestak8_scoutGlobalParT_HcsVsQCD)
            flavors.append(j.flavor)
            genidx.append(getattr(j, "genJetIdx", -1))
            #flavor = j.flavor
        self.out.fillBranch("n_bjets", len(HbqvsQCD))
            #self.out.fillBranch(prefix + "ak4_ak8_dR", j.closestFatJet_dr)
        self.out.fillBranch(prefix + "hadronFlavour", flavors)
        self.out.fillBranch(prefix + "ScoutingPFJetRecluster_genJetIdx", genidx)
        self.out.fillBranch(prefix + "scoutGloParT_Xbb", Xbb)
        self.out.fillBranch(prefix + "scoutGloParT_Xbc", Xbc)
        self.out.fillBranch(prefix + "scoutGloParT_Xbs", Xbs)
        self.out.fillBranch(prefix + "scoutGloParT_HbbVsQCD", HbbvsQCD)
        self.out.fillBranch(prefix + "scoutGloParT_HbqVsQCD", HbqvsQCD)
        self.out.fillBranch(prefix + "scoutGloParT_HcsVsQCD", HcsvsQCD)
            #self.out.fillBranch(prefix + "scoutGloParT_HbbVsHcc", j.closestak8_scoutGlobalParT_HbbVsHcc)
        
        # select the b-jet with the highest b-tagging score close to the muon
        bscores = [j.closestak8_scoutGlobalParT_HbqVsQCD for j in bjets]
        maxindex = bscores.index(max(bscores))
        event.bjet = bjets[maxindex]

        # fat-jet selection: select events where there is at least one fat jet
        # relatively far from the selected muon.
        if len(event.fatjets) == 0: return False
        probe_jets = [fj for fj in event.fatjets if abs(deltaPhi(fj, event.mu)) > 2]
        if len(probe_jets) == 0: return False

        # selection is done, now calculate event properties
        # (in particular of the probe jet)
        probe_jets = probe_jets[:1]
        self.loadGenHistory(event, probe_jets)
        self.evalMassRegression(probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)
        self.out.fillBranch("passMuTrig", passTrigger(event, ["HLT_Mu50"]))
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_phi", event.mu.phi)
        self.out.fillBranch("muon_trk_dxy", event.mu.trk_dxy)
        self.out.fillBranch("muon_trk_dz", event.mu.trk_dz)
        self.out.fillBranch("muon_trackIso", event.mu.trackIso)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())
        self.out.fillBranch("bjet_pt", event.bjet.pt)
        self.out.fillBranch("bjet_eta", event.bjet.eta)
        self.out.fillBranch("bjet_phi", event.bjet.phi)
        self.out.fillBranch("bjet_particleNet_prob_b", event.bjet.particleNet_prob_b)
        self.out.fillBranch("bjet_particleNet_prob_bb", event.bjet.particleNet_prob_bb)
        self.out.fillBranch("bjet_particleNet_prob_c", event.bjet.particleNet_prob_c)
        self.out.fillBranch("bjet_particleNet_prob_cc", event.bjet.particleNet_prob_cc)
        self.out.fillBranch("bjet_particleNet_prob_uds", event.bjet.particleNet_prob_uds)
        self.out.fillBranch("bjet_particleNet_prob_g", event.bjet.particleNet_prob_g)

        return True

def MuonTree_2016(): return MuonSampleProducerScouting(year=2016)
def MuonTree_2017(): return MuonSampleProducerScouting(year=2017)
def MuonTree_2018(): return MuonSampleProducerScouting(year=2018)
def MuonTree_2021(): return MuonSampleProducerScouting(year=2021)
def MuonTree_2022(): return MuonSampleProducerScouting(year=2022)
def MuonTree_2023(): return MuonSampleProducerScouting(year=2023)
def MuonTree_2024(): return MuonSampleProducerScouting(year=2024)
