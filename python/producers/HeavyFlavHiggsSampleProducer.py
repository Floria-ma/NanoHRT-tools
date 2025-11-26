from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from .HeavyFlavBaseProducer import HeavyFlavBaseProducer


class HiggsSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(HiggsSampleProducer, self).__init__(channel='higgs', **kwargs)
        # self._fill_sv = False # for QCD sample, do not fill SV info

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(HiggsSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        self.out.branch("nfj", "I")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""


        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        # require at least one fatjet
        if len(event.fatjets) < 1:
            return False

        self.loadGenHistory(event, event.fatjets)
        self.evalTagger(event, event.fatjets)
        self.evalMassRegression(event, event.fatjets)

        sorted_fatjet = sorted(event.fatjets, key=lambda x: x.dr_H)
        sorted_fatjet = sorted_fatjet[:1] # only select the fatjet nearest to higgs

        # match to SV
        self.selectSV(event)
        self.matchSVToFatJets(event, sorted_fatjet)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, sorted_fatjet)
        self.out.fillBranch("nfj", len(event.fatjets))

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def HiggsTree_2016(): return HiggsSampleProducer(year=2016)
def HiggsTree_2017(): return HiggsSampleProducer(year=2017)
def HiggsTree_2018(): return HiggsSampleProducer(year=2018)
