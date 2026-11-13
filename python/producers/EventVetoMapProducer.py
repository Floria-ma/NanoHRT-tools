import numpy as np
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from ..helpers.utils import closest
from ..helpers.jetmetCorrector import JetMETCorrector, _jmeJson

# def debug(msg, activated=True):
#     if activated:
#         print(' > ', msg)

'''https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Jet_veto_maps'''
'''https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps'''


class EventVetoMapProducer(Module):

    def __init__(self, year, **kwargs):
        self._opts = {}
        self._opts.update(**kwargs)
        self._usePuppiJets = True #self._opts['usePuppiJets']

        name = {2015: 'Summer19UL16_V1', 2016: 'Summer19UL16_V1',
                2017: 'Summer19UL17_V1', 2018: 'Summer19UL18_V1',
                20220: 'Summer22_23Sep2023_RunCD_V1', 20221: 'Summer22EE_23Sep2023_RunEFG_V1',
                20230: 'Summer23Prompt23_RunC_V1', 20231: 'Summer23BPixPrompt23_RunD_V1',
                2024: 'Summer24Prompt24_RunBCDEFGHI_V1'}[year]

        # era json directory from the shared JME pin table (see helpers/jetmetCorrector.py)
        filename = _jmeJson(JetMETCorrector._ERAS[year][0], 'jetvetomaps.json.gz')
        self.corr = correctionlib.CorrectionSet.from_file(filename)[name]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch('jetVetoMapEventVeto', "O")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        allJets = event._allJets
        muons = Collection(event, "Muon")

        eventVeto = False
        for j in allJets:
            # https://gitlab.cern.ch/cms-analysis/jme/jerc-application-tutorial/-/blob/60281d1527be8b1b4b7d392bffc7d23a5e4d5c03/JvmApplication.cpp#L24
            if not (j.pt > 15 and (j.jetId & 4)):
                # pt, tightId
                continue
            #if not self._usePuppiJets and not (j.pt > 50 or j.puId >= self.puID_WP):
            #    # apply jet puId only for CHS jets
            #    continue
            #if closest(j, muons)[1] < 0.2:
            #    continue
            if j.neEmEF + j.chEmEF > 0.9:
                continue
            # debug('jet phi:%.1f, eta:%.2f' % (j.phi, j.eta))
            jveto = self.corr.evaluate("jetvetomap", np.clip(j.eta, -5.19, 5.19), np.clip(j.phi, -3.14, 3.14))
            # debug('veto: %s' % str(jveto))
            if jveto:
                eventVeto = True
                break
        # debug('event veto: %s' % str(eventVeto))

        self.out.fillBranch('jetVetoMapEventVeto', eventVeto)

        return True


def eventVeto_2015():
    return EventVetoMapProducer(year=2015)


def eventVeto_2016():
    return EventVetoMapProducer(year=2016)


def eventVeto_2017():
    return EventVetoMapProducer(year=2017)


def eventVeto_2018():
    return EventVetoMapProducer(year=2018)
