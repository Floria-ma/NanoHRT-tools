import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from ..helpers.utils import deltaR

era_dict = {2015: '2016preVFP_UL', 2016: '2016postVFP_UL', 2017: '2017_UL', 2018: '2018_UL',
            20220: '2022_Summer22', 20221: '2022_Summer22EE', 20230: '2023_Summer23', 20231: '2023_Summer23BPix', 2024: '2024_Summer24'}


#keys for egamma SFs
era_label =  {2015: '2016preVFP', 2016: '2016postVFP', 2017: '2017', 2018: '2018',
              20220:"2022Re-recoBCD", 20221:"2022Re-recoE+PromptFG", 20230:"2023PromptC", 20231:"2023PromptD", 2024:"2024Prompt"}


def debug(msg, activated=False):
    if activated:
        print(' > ', msg)

def get_sf_root(hist, x, y, inflate_unc_out_bounds=True, no_unc=False):
    bin_x0 = hist.GetXaxis().FindFixBin(x)
    bin_y0 = hist.GetYaxis().FindFixBin(y)
    binX = int(np.clip(bin_x0, 1, hist.GetNbinsX()))
    binY = int(np.clip(bin_y0, 1, hist.GetNbinsY()))
    val = hist.GetBinContent(binX, binY)
    if no_unc:
        return np.array([val, val, val])
    err = hist.GetBinError(binX, binY)
    if binX != bin_x0 or binY != bin_y0:
        if inflate_unc_out_bounds:
            err *= 2.
    return np.array([val, val + err, val - err])


class SFHelper(object):

    def __init__(self, prefix, file=None, key=None, flipInputOrder=False):
        if (file is None) or (key is None):
            raise IOError("No input file or key provided")

        self.vars = key.split('/')[-1] if '/' in key else None
        filepath = os.path.join(prefix, file)
        if filepath.endswith('.root'):
            self._rtfile = ROOT.TFile(filepath)
            self.obj = self._rtfile.Get(key)
        else:
            raise RuntimeError(
                'Invalid file %s. Only root are supported.' % filepath)
        self.flipInputOrder=flipInputOrder

    def __call__(self, x, y, inflate_unc_out_bounds=True, no_unc=False):
        if self.obj is None:
            return np.ones(3)
        if self.flipInputOrder:
            return get_sf_root(self.obj, y, x, inflate_unc_out_bounds, no_unc)
        else:
            return get_sf_root(self.obj, x, y, inflate_unc_out_bounds, no_unc)

class ElectronSFProducer(Module, object):

    def __init__(self, year, **kwargs):
        self.year = year
        self.era = era_dict[self.year]
        if self.year<=2018:
            correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electron.json.gz'            
            self.corr = correctionlib.CorrectionSet.from_file(correction_file)['UL-Electron-ID-SF']
        elif self.year==2024:
            # the latest central EGM metadata on cvmfs, key '2024Prompt'
            correction_file = '/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/electron.json.gz'
            self.corr = correctionlib.CorrectionSet.from_file(correction_file)['Electron-ID-SF']
            #correction_file2 = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electronID_v1.json.gz'
            #self.corr2 = correctionlib.CorrectionSet.from_file(correction_file2)['Electron-ID-SF']
        else:
            correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electron.json.gz'
            self.corr = correctionlib.CorrectionSet.from_file(correction_file)['Electron-ID-SF']

        self._opts = {'fillSystWeights': True}
        self._opts.update(**kwargs)

    def SetIDWP(self, WPname):
        self._wp_ID = WPname

    def SetIsoWP(self, WPname):#useless for now
        self._wp_Iso = WPname
        
    def get_sf(self, sf_type, lep):
        if abs(lep.pdgId) != 11:
            raise RuntimeError('Input lepton is not a electron')

        wp = None
        if sf_type == 'Reco':
            if lep.pt < 20:
                wp = 'RecoBelow20'
            elif 20 <= lep.pt < 75:
                wp = 'Reco20to75'
            else:
                wp = 'RecoAbove75'
            
        elif sf_type == 'ID':
            if self._wp_ID=="Tight":
                assert lep.mvaIso_WP80
                wp = 'wp80iso'
            elif self._wp_ID=="Loose":
                assert lep.mvaIso_WP90
                wp = 'wp90iso'
            else:
                raise NameError(f'ele id working point "f{self._wp_ID}" not supported')

        label = era_label[self.year] 

        if self.year in [20230,20231] :
            scale_factors = np.array([self.corr.evaluate(label, syst, wp, lep.etaSC, lep.pt, lep.phi)
                                      for syst in ('sf', 'sfup', 'sfdown')])
        elif self.year==2024:
            #try:
            scale_factors = np.array([self.corr.evaluate(label, syst, wp, lep.etaSC, lep.pt)
                                          for syst in ('sf', 'sfup', 'sfdown')])
            #except:
            #    scale_factors = np.array([self.corr2.evaluate(label, syst, self._wp_ID, lep.etaSC, lep.pt)
            #                              for syst in ('sf', 'sfup', 'sfdown')])
                
        else:
            scale_factors = np.array([self.corr.evaluate(label, syst, self._wp_ID, lep.etaSC, lep.pt)
                                      for syst in ('sf', 'sfup', 'sfdown')])
        # debug(f'Electron pt:{lep.pt:.1f}, etaSC:{lep.etaSC:.2f}, {wp}, SF({sf_type}) = {scale_factors} (nom, up, down)')
        return scale_factors

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('elEffWeight', "F", limitedPrecision=10)
            self.out.branch('elIDWeight', "F", limitedPrecision=10)
            self.out.branch('elRecoWeight', "F", limitedPrecision=10)
            if self._opts['fillSystWeights']:
                self.out.branch('elEffWeight_UP', "F", limitedPrecision=10)
                self.out.branch('elEffWeight_DOWN', "F", limitedPrecision=10)
                self.out.branch('elIDWeight_UP', "F", limitedPrecision=10)
                self.out.branch('elIDWeight_DOWN', "F", limitedPrecision=10)
                self.out.branch('elRecoWeight_UP', "F", limitedPrecision=10)
                self.out.branch('elRecoWeight_DOWN', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgtReco = np.ones(3)
        wgtID = np.ones(3)

        for lep in event.selectedLeptons[:2]:
            # consider only up to two leading leptons
            if abs(lep.pdgId) != 11:
                continue
            wgtReco *= self.get_sf('Reco', lep)
            wgtID *= self.get_sf('ID', lep)

        eventWgt = wgtReco * wgtID
        self.out.fillBranch('elEffWeight', eventWgt[0])
        self.out.fillBranch('elIDWeight', wgtID[0])
        self.out.fillBranch('elRecoWeight', wgtReco[0])
        if self._opts['fillSystWeights']:
            self.out.fillBranch('elEffWeight_UP', eventWgt[1])
            self.out.fillBranch('elEffWeight_DOWN', eventWgt[2])
            self.out.fillBranch('elIDWeight_UP', wgtID[1])
            self.out.fillBranch('elIDWeight_DOWN', wgtID[2])
            self.out.fillBranch('elRecoWeight_UP', wgtReco[1])
            self.out.fillBranch('elRecoWeight_DOWN', wgtReco[2])

        return True


class MuonSFProducer(Module, object):

    def __init__(self, year, **kwargs):
        self.year = year
        self.era = era_dict[self.year]
        if self.year == 2024:
            correction_file = '/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/muon_Z.json.gz'
        else:
            correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_Z.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)
        self._opts = {'fillSystWeights': True}
        self._opts.update(**kwargs)

    def SetIDWP(self, WPname):
        self._wp_ID = WPname

    def SetIsoWP(self, WPname):
        self._wp_Iso = WPname

    def get_sf(self, sf_type, lep):
        if abs(lep.pdgId) != 13:
            raise RuntimeError('Input lepton is not a muon')
        if sf_type == 'ID':
            if self._wp_ID=="Tight":
                assert lep.tightId
            elif self._wp_ID=="Loose":
                assert lep.looseId
            else:
                raise NameError(f'mu id working point "f{self._wp_ID}" not supported')

            if self.year<=2018:
                key = f'NUM_{self._wp_ID}ID_DEN_genTracks'
            else:
                key = f'NUM_{self._wp_ID}ID_DEN_TrackerMuons'
        elif sf_type == 'Iso':
            if self._wp_Iso=="Tight": #cut at 0.06
                assert lep.pfIsoId>=5
            elif self._wp_Iso=="Loose": # cut at 0.25
                assert lep.pfIsoId>=2
            else:
                raise NameError(f'mu isolation working point "f{self._wp_Iso}" not supported')

            if self.year<=2018:
                if self._wp_ID=="Tight":
                    key = f'NUM_{self._wp_Iso}RelIso_DEN_TightIDandIPCut'
                else:
                    key = f'NUM_{self._wp_Iso}RelIso_DEN_{self._wp_ID}ID'
            else:
                key = f'NUM_{self._wp_Iso}PFIso_DEN_{self._wp_ID}ID'

        scale_factors = np.array([self.corr[key].evaluate(abs(lep.eta), lep.pt, syst)
                                  for syst in ('nominal', 'systup', 'systdown')])
        # debug(f'Muon pt:{lep.pt:.1f}, eta:{lep.eta:.2f}, wp:({lep._wp_ID}, {lep._wp_Iso}), SF({sf_type}) = {scale_factors} (nom, up, down)')
        return scale_factors

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('muEffWeight', "F", limitedPrecision=10)
            self.out.branch('muIsoWeight', "F", limitedPrecision=10)
            self.out.branch('muIDWeight', "F", limitedPrecision=10)
            if self._opts['fillSystWeights']:
                self.out.branch('muEffWeight_UP', "F", limitedPrecision=10)
                self.out.branch('muEffWeight_DOWN', "F", limitedPrecision=10)
                self.out.branch('muIsoWeight_UP', "F", limitedPrecision=10)
                self.out.branch('muIsoWeight_DOWN', "F", limitedPrecision=10)
                self.out.branch('muIDWeight_UP', "F", limitedPrecision=10)
                self.out.branch('muIDWeight_DOWN', "F", limitedPrecision=10)


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgtID = np.ones(3)
        wgtIso = np.ones(3)

        for lep in event.selectedLeptons[:2]:
            # consider only up to two leading leptons
            if abs(lep.pdgId) != 13:
                continue
            wgtID *= self.get_sf('ID', lep)
            wgtIso *= self.get_sf('Iso', lep)

        eventWgt = wgtID * wgtIso
        self.out.fillBranch('muEffWeight', eventWgt[0])
        self.out.fillBranch('muIsoWeight', wgtIso[0])
        self.out.fillBranch('muIDWeight', wgtID[0])
        if self._opts['fillSystWeights']:
            self.out.fillBranch('muEffWeight_UP', eventWgt[1])
            self.out.fillBranch('muEffWeight_DOWN', eventWgt[2])
            self.out.fillBranch('muIsoWeight_UP', wgtIso[1])
            self.out.fillBranch('muIsoWeight_DOWN', wgtIso[2])
            self.out.fillBranch('muIDWeight_UP', wgtID[1])
            self.out.fillBranch('muIDWeight_DOWN', wgtID[2])

        return True


class TriggerSF(Module, object):

    def __init__(self, year, **kwargs):
        self.year = year
        self.era = era_dict[self.year]
        self._opts = {'fillSystWeights': True}
        self._opts.update(**kwargs)

        eleptthrs={2015:27., 2016:27., 2017:32., 2018:32., 20220:30., 20221:30., 20230:30., 20231:30., 2024:30.}
        muptthrs={2015:24., 2016:24., 2017:27., 2018:24., 20220:24., 20221:24., 20230:24., 20231:24., 2024:24.}
        self.eleptmin=eleptthrs[self.year]
        self.muptmin=muptthrs[self.year]
        
    def SetChannel(self, channel):
        self.channel = channel

    def LoadSFs(self):
        if self.channel == '0L':
            pass #correction added in the histo production step
        if self.channel == '1L':
            # mu
            if self.year<2017:
                mukey='NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight'
            elif self.year==2017:
                mukey='NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight'
            else:
                mukey='NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight'
            
            if self.year == 2024:
                mufile = '/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/muon_Z.json.gz'
            else:
                mufile = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_Z.json.gz'
            print(f"Loading muon trigger SFs from {mufile}")
            self.corr_mu = correctionlib.CorrectionSet.from_file(mufile)[mukey]

            # ele
            self.path=None
            if self.year<2017 or self.year>2018:
                if self.year<2017:
                    # from https://indico.cern.ch/event/1470465/#7-trigger-sfs-for-hig-24-018
                    elepath="$CMSSW_BASE/src/PhysicsTools/NanoTrees/data/scale_factors/trigger/singleEle/1L/"
                    fname=f'{self.era}.json'
                    self.corr_el = correctionlib.CorrectionSet.from_file(os.path.expandvars(f'{elepath}/{fname}'))
                else:
                    # official EGM SFs 
                    if self.year == 2024:
                        elepath='/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/'
                    else:
                        elepath=f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/'
                    fname='electronHlt.json.gz'
                    self.path="HLT_SF_Ele30_MVAiso80ID"
                    self.corr_el = correctionlib.CorrectionSet.from_file(os.path.expandvars(f'{elepath}/{fname}'))["Electron-HLT-SF"]
                print(f"Loading electron trigger SFs from {elepath}/{fname}")

            else:
                # from https://indico.cern.ch/event/1288547/#25-trigger-sfs-for-gamma-gamma
                elepath="$CMSSW_BASE/src/PhysicsTools/NanoTrees/data/scale_factors/trigger/singleEle/1L/"
                fname=f'{self.era}.root'
                print(f"Loading electron trigger SFs from {elepath}/{fname}")
                self.corr_el = SFHelper(os.path.expandvars(elepath), fname, 'SF2D')
            
        elif self.channel == '2L':
            #mu
            if self.year<2017:
                mukey='NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight'
            elif self.year==2017:
                mukey='NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight'
            else:
                mukey='NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight'

            if self.year == 2024:
                mufile = '/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/muon_Z.json.gz'
            else:
                mufile = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_Z.json.gz'
            self.corr_mu = correctionlib.CorrectionSet.from_file(mufile)[mukey]
            
            #ele
            if self.year<=2018:
                #from https://indico.cern.ch/event/1146225/contributions/4836371/attachments/2429668/4197633/TauStar_SingleElSFsApprovalTalk_Updated12May.pdf
                elepath="$CMSSW_BASE/src/PhysicsTools/NanoTrees/data/scale_factors/trigger/singleEle/2L/"
                fname=f'{self.era}.root'
                print(f"Loading electron trigger SFs from {elepath}/{fname}")
                self.corr_el =  SFHelper(os.path.expandvars(elepath), fname, 'EGamma_SF2D', flipInputOrder=True)
                self.MCeff_el = SFHelper(os.path.expandvars(elepath), fname, 'EGamma_EffMC2D', flipInputOrder=True)
            else:
                # official EGM SFs
                if self.year == 2024:
                    elepath='/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/'
                else:
                    elepath=f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/'
                fname='electronHlt.json.gz'
                self.path="HLT_SF_Ele30_MVAiso90ID"
                print(f"Loading electron trigger SFs from {elepath}/{fname}")
                self.corr_el = correctionlib.CorrectionSet.from_file(os.path.expandvars(f'{elepath}/{fname}'))["Electron-HLT-SF"]
                self.MCeff_el = correctionlib.CorrectionSet.from_file(os.path.expandvars(f'{elepath}/{fname}'))["Electron-HLT-McEff"]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        # TrigObj_id/filterBits, without kinematics
        self.hasTrigObjpT = bool(inputTree.GetBranch('TrigObj_pt'))
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('elTrigWeight', "F", limitedPrecision=10)
            self.out.branch('muTrigWeight', "F", limitedPrecision=10)
            if self._opts['fillSystWeights']:
                self.out.branch('elTrigWeight_UP', "F", limitedPrecision=10)
                self.out.branch('elTrigWeight_DOWN', "F", limitedPrecision=10)
                self.out.branch('muTrigWeight_UP', "F", limitedPrecision=10)
                self.out.branch('muTrigWeight_DOWN', "F", limitedPrecision=10)

    def get_muon_trigger_SF(self, eta, pt):
        # 2024 muon trigger SFs are binned in signed eta, other years in abs(eta)
        if self.year != 2024:
            eta = abs(eta)
        if (self.year!=2017 and pt<26):
            print("warning found muon with pt<26 GeV --> I will take the sf for pt=26 and inflate its uncertainty")
            # systup/systdown are the shifted SFs (not the uncertainty)
            nom = self.corr_mu.evaluate(eta, 26.1, "nominal")
            trigWgt = np.array([nom,
                                nom + 2*(self.corr_mu.evaluate(eta, 26.1, "systup") - nom),
                                nom + 2*(self.corr_mu.evaluate(eta, 26.1, "systdown") - nom)])
        elif (self.year==2017 and pt<29):
            print("warning found muon with pt<29 GeV --> I will take the sf for pt=29 and inflate its uncertainty")
            # systup/systdown are the shifted SFs (not the uncertainty)
            nom = self.corr_mu.evaluate(eta, 29.1, "nominal")
            trigWgt = np.array([nom,
                                nom + 2*(self.corr_mu.evaluate(eta, 29.1, "systup") - nom),
                                nom + 2*(self.corr_mu.evaluate(eta, 29.1, "systdown") - nom)])
        else:
            trigWgt = np.array([self.corr_mu.evaluate(eta, pt, syst)
                                for syst in ('nominal', 'systup', 'systdown')])
        return trigWgt

    def get_muon_MCeff(self, eta, pt):
        # 2024 muon json has MC trigger efficiencies
        if self.year == 2024:
            if pt < 26:
                return self.corr_mu.evaluate(eta, 26.1, "nominal_MCeff")
            else:
                return self.corr_mu.evaluate(eta, pt, "nominal_MCeff")
        return 1.0 
    
    def get_ele_trigger_SF(self, etaSC, pt):
        # use different methods depeding on the correction format
        if isinstance(self.corr_el,SFHelper):
            trigWgt = np.array( self.corr_el(pt, etaSC) )
        else:
            try:            
                if (self.year>2018 and pt<25):
                    print("warning: found electron  with pt<25 GeV --> I will take the sf for pt=25 and inflate its uncertainty")
                    trigWgt = np.array([self.corr_el.evaluate(era_label[self.year], 'sf',     self.path, etaSC, 25.1),
                                        self.corr_el.evaluate(era_label[self.year], 'sfup',   self.path, etaSC, 25.1),
                                        self.corr_el.evaluate(era_label[self.year], 'sfdown', self.path, etaSC, 25.1)])
                                        
                else:
                    trigWgt = np.array([self.corr_el.evaluate(era_label[self.year], syst, self.path, etaSC, pt)
                                        for syst in ('sf', 'sfup', 'sfdown')])
            except:
                try:
                    if (abs(etaSC)>=1.4442 and abs(etaSC)<=1.556):
                        print("warning: found electron in barrel-endcap transition region --> I will set trigger SF to 1 and apply a conservative 6% uncertainty")
                        trigWgt = np.array([1.0, 1.06, 0.94])
                    else:
                        trigWgt = np.array([self.corr_el[syst].evaluate(pt, etaSC)
                                            for syst in ('trigger_SF_nom','trigger_SF_up','trigger_SF_down')])
                except:
                    raise RuntimeError("Fail in accessing trigger SF")
        return trigWgt

    def get_ele_MCeff(self, etaSC, pt):
        # use different methods depeding on the correction format
        if isinstance(self.MCeff_el,SFHelper):
            eff = np.array( self.MCeff_el(pt, etaSC, True, True) )
            return eff[0]
        else:
            if (self.year>2018 and pt<25):
                print("warning: found electron  with pt<25 GeV --> I will take the sf for pt=25 and inflate its uncertainty")
                return self.MCeff_el.evaluate(era_label[self.year], "nom", self.path, etaSC, 25.1)
            else:
                return self.MCeff_el.evaluate(era_label[self.year], "nom", self.path, etaSC, pt)
    
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if not self.isMC:
            return True

        trigWgt = np.ones(3)
        fl=0

        if self.channel == '0L':
            return True

        elif self.channel == '1L':
            lep = event.selectedLeptons[0]
            fl=abs(lep.pdgId)
            if fl == 13:
                trigWgt=self.get_muon_trigger_SF(lep.eta, lep.pt)
            elif fl == 11:
                trigWgt=self.get_ele_trigger_SF(lep.etaSC, lep.pt)
            else:
                raise ValueError("Unsupported lepton")

        elif self.channel == '2L':

            #lepton flavour
            assert abs(event.selectedLeptons[0].pdgId)==abs(event.selectedLeptons[1].pdgId)
            fl=abs(event.selectedLeptons[0].pdgId)
            assert(fl in [11,13])

            #find what leptons fired the trigger by requiring deltaR matching to the trigger objects
            SFs=[]
            EFFs=[]
            for lep in event.selectedLeptons[:2]:
                if not self.hasTrigObjpT:
                    # no trigger object pT, no matching possible: SFs stays empty
                    # TODO: discuss how to fix this
                    break
                for trigObj in event.trigobjs:
                    if abs(trigObj.id)!=fl: # trigobj must be an electron or a muon
                        continue
                    if fl==11:
                        if not(trigObj.pt>self.eleptmin and trigObj.filterBits&2): # pt thr and tight id
                            continue
                    elif fl==13:
                        if not( trigObj.pt>self.muptmin and trigObj.filterBits&2): # pt thr and isolation
                             continue
                    if deltaR(lep,trigObj)<0.1: # deltaR matching
                        if fl==13:
                            SFs.append( self.get_muon_trigger_SF(lep.eta, lep.pt) )
                            EFFs.append( self.get_muon_MCeff(lep.eta, lep.pt) )
                        else:
                            SFs.append( self.get_ele_trigger_SF(lep.etaSC, lep.pt) )
                            EFFs.append( self.get_ele_MCeff(lep.etaSC, lep.pt) )
                        break
                            
            # calculate SF depending on the number of leptons firing the trigger
            if len(SFs)==0:
                # very rare case in which probably a third lepton fired the trigger 
                trigWgt = [1.0, 1.05, 0.95]
            elif len(SFs)==1:
                trigWgt = SFs[0]
            else:
                mc1, mc2 = EFFs[0], EFFs[1]
                sf1, sf2 = SFs[0], SFs[1]
                mc_eff = mc1 + mc2 - mc1 * mc2
                data_eff = (sf1 * mc1) + (sf2 * mc2) - (sf1 * mc1) * (sf2 * mc2)
                trigWgt = data_eff / mc_eff

        # fill trig weight associated to affected lepton and set the other to one
        if fl==11:
            eltrigWgt = trigWgt
            mutrigWgt = np.ones(3)
        else:
            eltrigWgt = np.ones(3)
            mutrigWgt = trigWgt
        self.out.fillBranch('elTrigWeight', eltrigWgt[0])
        self.out.fillBranch('muTrigWeight', mutrigWgt[0])
        if self._opts['fillSystWeights']:
            self.out.fillBranch('elTrigWeight_UP', eltrigWgt[1])
            self.out.fillBranch('elTrigWeight_DOWN', eltrigWgt[2])
            self.out.fillBranch('muTrigWeight_UP', mutrigWgt[1])
            self.out.fillBranch('muTrigWeight_DOWN', mutrigWgt[2])

        return True

