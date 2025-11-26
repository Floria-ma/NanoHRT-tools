from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ..helpers.utils import deltaPhi, polarP4, configLogger
from ..helpers.muonCorrector import MuonScaleResCorrector
from ..helpers.triggerHelper import passTrigger
from .HeavyFlavBaseProducer import HeavyFlavBaseProducer

import logging
logger = logging.getLogger('nano')
configLogger('nano', loglevel=logging.INFO)

def convert_prob(jet, sigs, bkgs=None, prefix=''):
    if not jet:
        return -1

    def get(name):
        if isinstance(jet, dict):
            return jet[name]
        else:
            return getattr(jet, name)

    if bkgs is None:
        bkgs = [prefix + n for n in ['QCDbb', 'QCDb', 'QCDcc', 'QCDc', 'QCDothers']]
    else:
        if not isinstance(bkgs, (list, tuple)):
            bkgs = [bkgs]
        bkgs = [prefix + n for n in bkgs]
    bkgsum = sum([get(name) for name in bkgs])

    if sigs is None:
        return bkgsum
    else:
        if not isinstance(sigs, (list, tuple)):
            sigs = [sigs]
        sigs = [prefix + n for n in sigs]
    sigsum = sum([get(name) for name in sigs])

    try:
        return sigsum / (sigsum + bkgsum)
    except ZeroDivisionError:
        return -1



class MuonSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(MuonSampleProducer, self).__init__(channel='muon', **kwargs)

        # ParticleNetAK4 -- exclusive b- and c-tagging categories
        # 5x: b-tagged; 4x: c-tagged; 0: light
        if self.year in (2017, 2018):
            self.jetTagWPs = {
                54: '(pn_b_plus_c>0.5) & (pn_b_vs_c>0.99)',
                53: '(pn_b_plus_c>0.5) & (0.96<pn_b_vs_c<=0.99)',
                52: '(pn_b_plus_c>0.5) & (0.88<pn_b_vs_c<=0.96)',
                51: '(pn_b_plus_c>0.5) & (0.70<pn_b_vs_c<=0.88)',
                50: '(pn_b_plus_c>0.5) & (0.40<pn_b_vs_c<=0.70)',

                44: '(pn_b_plus_c>0.5) & (pn_b_vs_c<=0.05)',
                43: '(pn_b_plus_c>0.5) & (0.05<pn_b_vs_c<=0.15)',
                42: '(pn_b_plus_c>0.5) & (0.15<pn_b_vs_c<=0.40)',
                41: '(0.2<pn_b_plus_c<=0.5)',
                40: '(0.1<pn_b_plus_c<=0.2)',

                0: '(pn_b_plus_c<=0.1)',
            }
        elif self.year in (2015, 2016):
            self.jetTagWPs = {
                54: '(pn_b_plus_c>0.35) & (pn_b_vs_c>0.99)',
                53: '(pn_b_plus_c>0.35) & (0.96<pn_b_vs_c<=0.99)',
                52: '(pn_b_plus_c>0.35) & (0.88<pn_b_vs_c<=0.96)',
                51: '(pn_b_plus_c>0.35) & (0.70<pn_b_vs_c<=0.88)',
                50: '(pn_b_plus_c>0.35) & (0.40<pn_b_vs_c<=0.70)',

                44: '(pn_b_plus_c>0.35) & (pn_b_vs_c<=0.05)',
                43: '(pn_b_plus_c>0.35) & (0.05<pn_b_vs_c<=0.15)',
                42: '(pn_b_plus_c>0.35) & (0.15<pn_b_vs_c<=0.40)',
                41: '(0.17<pn_b_plus_c<=0.35)',
                40: '(0.1<pn_b_plus_c<=0.17)',

                0: '(pn_b_plus_c<=0.1)',
            }
        else:
            self.jetTagWPs = {}


    def evalJetTag(self, j, default=0):
        for wp, expr in self.jetTagWPs.items():
            if eval(expr, j.__dict__):
                return wp
        return default

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(MuonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # trigger variables
        self.out.branch("passMuTrig", "O")

        logger.info('Running year %s ',
            str(self.year))

        # event variables
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_miniIso", "F")
        self.out.branch("leptonicW_pt", "F")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # muon selection
        event._allMuons = Collection(event, "ScoutingMuonVtx")
        event.muons = [mu for mu in event._allMuons if mu.pt > 55 and abs(mu.eta) < 2.4 and abs(
            mu.dxy) < 0.2 and abs(mu.dz) < 0.5 and mu.tightId and mu.miniPFRelIso_all < 0.10]

        if len(event.muons) != 1:
            return False

        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        # met selection
        if event.met.pt < 50.0:
            return False

        # leptonic W pt cut
        event.mu = event.muons[0]
        event.leptonicW = polarP4(event.mu) + event.met.p4()
        if event.leptonicW.Pt() < 100.0:
            return False

        # at least one b-jet, in the same hemisphere of the muon
        if (self.year < 2019):
            for j in event.ak4jets:
                # attach ParticleNet scores
                j.pn_b = convert_prob(j, ['b', 'bb'], ['c', 'cc', 'uds', 'g'], 'ParticleNetAK4_prob')
                j.pn_c = convert_prob(j, ['c', 'cc'], ['b', 'bb', 'uds', 'g'], 'ParticleNetAK4_prob')
                j.pn_uds = convert_prob(j, 'uds', ['b', 'bb', 'c', 'cc', 'g'], 'ParticleNetAK4_prob')
                j.pn_g = convert_prob(j, 'g', ['b', 'bb', 'c', 'cc', 'uds'], 'ParticleNetAK4_prob')
                j.pn_b_plus_c = j.pn_b + j.pn_c
                j.pn_b_vs_c = j.pn_b / j.pn_b_plus_c
                j.tag = self.evalJetTag(j)

            event.bjets = [j for j in event.ak4jets if j.tag >= 51 and
                            abs(deltaPhi(j, event.mu)) < 2]

        else:
            event.bjets = [j for j in event.ak4jets if j.btagPNetB > self.PNet_WP_M and
                        abs(deltaPhi(j, event.mu)) < 2]        
        
        if len(event.bjets) == 0:
            return False

        # require fatjet away from the muon
        probe_jets = [fj for fj in event.fatjets if abs(deltaPhi(fj, event.mu)) > 2]
        if len(probe_jets) == 0:
            return False

        probe_jets = probe_jets[:1]
        self.loadGenHistory(event, probe_jets)
        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)

        # fill
        self.out.fillBranch("passMuTrig", passTrigger(event, ['HLT_Mu50']))
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_miniIso", event.mu.miniPFRelIso_all)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def MuonTree_2016(): return MuonSampleProducer(year=2016)
def MuonTree_2017(): return MuonSampleProducer(year=2017)
def MuonTree_2018(): return MuonSampleProducer(year=2018)
def MuonTree_2021(): return MuonSampleProducer(year=2021)
def MuonTree_2022(): return MuonSampleProducer(year=2022)
def MuonTree_2023(): return MuonSampleProducer(year=2023)
def MuonTree_2024(): return MuonSampleProducer(year=2024)
