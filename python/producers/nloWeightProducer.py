import ROOT
import os
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.utils import clip


class NLOWeightProducer(Module, object):

    def __init__(self, year, **kwargs):
        self._opts = {
            'fillSystWeights': True,
        }
        self._opts.update(**kwargs)

        self._vhewk_histnames = ['SignalWeight_nloEWK_rebin',
                                 'SignalWeight_nloEWK_up_rebin', 'SignalWeight_nloEWK_down_rebin']

        self._fWpH = os.path.expandvars(
            '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/gen/Wp_nloEWK_weight_unnormalized.root')
        self._hWpH = None
        self._normfactor_WpH = (474233.6 / 435142.8, 474233.6 / 436091.1, 474233.6 / 434194.5)

        self._fWmH = os.path.expandvars(
            '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/gen/Wm_nloEWK_weight_unnormalized.root')
        self._hWmH = None
        self._normfactor_WmH = (256350.0 / 236154.8, 256350.0 / 236568.9, 256350.0 / 235740.6)

        self._fZllH = os.path.expandvars(
            '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/gen/Zll_nloEWK_weight_unnormalized.root')
        self._hZllH = None
        # self._normfactor_ZllH = (367873.6 / 334293.4, 335519.0 / 335519.0, 335519.0 / 333067.7)
        self._normfactor_ZllH = (367873.6 / 334293.4, 367873.6 / 335519.0, 367873.6 / 333067.7)

        self._fZnnH = os.path.expandvars(
            '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/gen/Znn_nloEWK_weight_unnormalized.root')
        self._hZnnH = None
        self._normfactor_ZnnH = (225368.7 / 215478.1, 225368.7 / 215612.1, 225368.7 / 215344.1)

        self._fVV = os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/gen/vv_nnlo_wgt.root')
        self._hWW = None
        self._hWZ = None
        self._hZZ = None

    def loadHist(self, filename, hname, normfactor=1.):
        tf = ROOT.TFile.Open(filename)
        hist = tf.Get(hname)
        hist.SetDirectory(0)
        tf.Close()
        hist.Scale(normfactor)
        return hist

    def loadHistGroup(self, filename, histnames, normfactors=None):
        if isinstance(histnames, dict):
            return {n: self.loadHist(filename, histnames[n], 1. if normfactors is None else normfactors[n]) for n in histnames}
        else:
            return [self.loadHist(filename, n, norm) for n, norm in zip(histnames, len(histnames) * [1.] if normfactors is None else normfactors)]

    def getBinContent(self, hist, xvalue, useUnderFlow=False, useOverFlow=False):
        lower = 0 if useUnderFlow else 1
        upper = hist.GetNbinsX() + 1 if useOverFlow else hist.GetNbinsX()
        bin = int(clip(hist.GetXaxis().FindFixBin(xvalue), lower, upper))
        return hist.GetBinContent(bin)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('genV_pt', "F")
            self.out.branch('v1_pt', "F")
            self.out.branch('v2_pt', "F")
            self.out.branch('topptWeight', "F")
            self.out.branch('vptWeightEWK', "F")
            self.out.branch('wptWeightEWK', "F")
            self.out.branch('zptWeightEWK', "F")
            self.out.branch('vhWeightEWK', "F")
            self.out.branch('vhWeightEWK_UP', "F")
            self.out.branch('vhWeightEWK_DOWN', "F")
            self.out.branch('vvWeightNNLO', "F")
            self.out.branch('vvWeightNNLO_UP', "F")
            self.out.branch('vvWeightNNLO_DOWN', "F")

            # accumulation hists
            self.h_genwgt = ROOT.TH1D('genweight', 'genweight', 1, 0, 1)
            self.h_ttbar_topptwgt = ROOT.TH1D('TT_toppt_weight', 'TT_toppt_weight', 1, 0, 1)
            self.h_vh_ewkwgt = ROOT.TH1D('VH_EWK_weight', 'VH_EWK_Weight', 1, 0, 1)
            self.h_vh_ewkwgt_UP = ROOT.TH1D('VH_EWK_weight_UP', 'VH_EWK_Weight_UP', 1, 0, 1)
            self.h_vh_ewkwgt_DOWN = ROOT.TH1D('VH_EWK_weight_DOWN', 'VH_EWK_Weight_DOWN', 1, 0, 1)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.isMC:
            cwd = ROOT.gDirectory
            outputFile.cd()
            self.h_genwgt.Write()
            self.h_ttbar_topptwgt.Write()
            self.h_vh_ewkwgt.Write()
            self.h_vh_ewkwgt_UP.Write()
            self.h_vh_ewkwgt_DOWN.Write()
            cwd.cd()

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

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

        def isInvisible(gp):
            if len(gp.dauIdx) == 0:
                raise ValueError('Particle has no daughters!')
            isInv = True
            for idx in gp.dauIdx:
                if abs(genparts[idx].pdgId) not in (12, 14, 16):
                    isInv = False
                    break
            return isInv

        genTops = []
        genWs = []
        genZs = []
        genHs = []

        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                # 13: isLastCopy
                continue
            if abs(gp.pdgId) == 6:
                genTops.append(gp)
            elif abs(gp.pdgId) == 24:
                genWs.append(gp)
            elif abs(gp.pdgId) == 23:
                genZs.append(gp)
            elif abs(gp.pdgId) == 25:
                genHs.append(gp)

        genV_pt = 0.
        vlead_pt = 0.
        vsub_pt = 0.
        topptWeight = 1.
        vptWeightEWK = 1.
        try:
            wptWeightEWK = (-0.830041 + 7.93714 * (math.pow(event.LHE_Vpt + 877.978, -0.213831))
                            ) if (event.LHE_Vpt >= 100 and event.LHE_Vpt < 3000) else 1.
            zptWeightEWK = (-0.1808051 + 6.04146 * (math.pow(event.LHE_Vpt + 759.098, -0.242556))
                            ) if (event.LHE_Vpt >= 100 and event.LHE_Vpt < 3000) else 1.
        except RuntimeError as e:
            wptWeightEWK = 1.
            zptWeightEWK = 1.
        vhWeightEWK = 1.
        vhWeightEWK_UP = 1.
        vhWeightEWK_DOWN = 1.
        vvWeightNNLO = 1.
        vvWeightNNLO_UP = 1.
        vvWeightNNLO_DOWN = 1.

        if len(genTops) == 2:
            # ttbar (+X ?)
            def wgt(pt):
                return math.exp(0.0615 - 0.0005 * clip(pt, 0, 800))

            def wgt_nnlo(pt):
                x = clip(pt, 0, 2000)
                return 0.103 * math.exp(-0.0118 * x) - 0.000134 * x + 0.973

            # topptWeight = math.sqrt(wgt(genTops[0].pt) * wgt(genTops[1].pt))
            topptWeight = math.sqrt(wgt_nnlo(genTops[0].pt) * wgt_nnlo(genTops[1].pt))

        elif len(genTops) == 0:
            if len(genWs) == 1 and len(genZs) == 0 and len(genHs) == 0:
                # W+jets
                wpt = genWs[0].pt
                genV_pt = wpt
                vptWeightEWK = (-0.830041 + 7.93714 * (math.pow(wpt + 877.978, -0.213831))
                                ) if (wpt >= 100 and wpt < 3000) else 1.
            elif len(genWs) == 0 and len(genZs) == 1 and len(genHs) == 0:
                # Z+jets
                zpt = genZs[0].pt
                genV_pt = zpt
                vptWeightEWK = (-0.1808051 + 6.04146 * (math.pow(zpt + 759.098, -0.242556))
                                ) if (zpt >= 100 and zpt < 3000) else 1.
            elif len(genWs) == 1 and len(genZs) == 0 and len(genHs) == 1:
                # WH
                if genWs[0].pdgId > 0:
                    # W+
                    if self._hWpH is None:
                        self._hWpH = self.loadHistGroup(self._fWpH, self._vhewk_histnames, self._normfactor_WpH)
                    vhWeightEWK, vhWeightEWK_UP, vhWeightEWK_DOWN = [
                        self.getBinContent(h, genWs[0].pt) for h in self._hWpH]
                else:
                    # W-
                    if self._hWmH is None:
                        self._hWmH = self.loadHistGroup(self._fWmH, self._vhewk_histnames, self._normfactor_WmH)
                    vhWeightEWK, vhWeightEWK_UP, vhWeightEWK_DOWN = [
                        self.getBinContent(h, genWs[0].pt) for h in self._hWmH]
            elif len(genWs) == 0 and len(genZs) == 1 and len(genHs) == 1:
                # ZH
                if isInvisible(genZs[0]):
                    if self._hZnnH is None:
                        self._hZnnH = self.loadHistGroup(self._fZnnH, self._vhewk_histnames, self._normfactor_ZnnH)
                    vhWeightEWK, vhWeightEWK_UP, vhWeightEWK_DOWN = [
                        self.getBinContent(h, genZs[0].pt) for h in self._hZnnH]
                else:
                    # this can also be Z(qq)H though...
                    if self._hZllH is None:
                        self._hZllH = self.loadHistGroup(self._fZllH, self._vhewk_histnames, self._normfactor_ZllH)
                    vhWeightEWK, vhWeightEWK_UP, vhWeightEWK_DOWN = [
                        self.getBinContent(h, genZs[0].pt) for h in self._hZllH]
            elif len(genWs) == 2 and len(genZs) == 0 and len(genHs) == 0:
                # WW
                vlead_pt, vsub_pt = sorted([v.pt for v in genWs + genZs], reverse=True)
                if self._hWW is None:
                    self._hWW = self.loadHistGroup(self._fVV, ['WW_nom', 'WW_up', 'WW_down'])
                vvWeightNNLO, vvWeightNNLO_UP, vvWeightNNLO_DOWN = [self.getBinContent(h, vsub_pt) for h in self._hWW]
            elif len(genWs) == 1 and len(genZs) == 1 and len(genHs) == 0:
                # WZ
                vlead_pt, vsub_pt = sorted([v.pt for v in genWs + genZs], reverse=True)
                if self._hWZ is None:
                    self._hWZ = self.loadHistGroup(self._fVV, ['WZ_nom', 'WZ_up', 'WZ_down'])
                vvWeightNNLO, vvWeightNNLO_UP, vvWeightNNLO_DOWN = [self.getBinContent(h, vsub_pt) for h in self._hWZ]
            elif len(genWs) == 0 and len(genZs) == 2 and len(genHs) == 0:
                # ZZ
                vlead_pt, vsub_pt = sorted([v.pt for v in genWs + genZs], reverse=True)
                if self._hZZ is None:
                    self._hZZ = self.loadHistGroup(self._fVV, ['ZZ_nom', 'ZZ_up', 'ZZ_down'])
                vvWeightNNLO, vvWeightNNLO_UP, vvWeightNNLO_DOWN = [self.getBinContent(h, vsub_pt) for h in self._hZZ]

        self.out.fillBranch('genV_pt', genV_pt)
        self.out.fillBranch('v1_pt', vlead_pt)
        self.out.fillBranch('v2_pt', vsub_pt)
        self.out.fillBranch('topptWeight', topptWeight)
        self.out.fillBranch('vptWeightEWK', vptWeightEWK)
        self.out.fillBranch('wptWeightEWK', wptWeightEWK)
        self.out.fillBranch('zptWeightEWK', zptWeightEWK)
        self.out.fillBranch('vhWeightEWK', vhWeightEWK)
        self.out.fillBranch('vhWeightEWK_UP', vhWeightEWK_UP)
        self.out.fillBranch('vhWeightEWK_DOWN', vhWeightEWK_DOWN)
        self.out.fillBranch('vvWeightNNLO', vvWeightNNLO)
        self.out.fillBranch('vvWeightNNLO_UP', vvWeightNNLO_UP)
        self.out.fillBranch('vvWeightNNLO_DOWN', vvWeightNNLO_DOWN)
        # accumulation hists
        self.h_genwgt.Fill(0.5, event.genWeight)
        self.h_ttbar_topptwgt.Fill(0.5, topptWeight * event.genWeight)
        self.h_vh_ewkwgt.Fill(0.5, vhWeightEWK * event.genWeight)
        self.h_vh_ewkwgt_UP.Fill(0.5, vhWeightEWK_UP * event.genWeight)
        self.h_vh_ewkwgt_DOWN.Fill(0.5, vhWeightEWK_DOWN * event.genWeight)

        return True
