import logging
import os
import ROOT
import correctionlib
ROOT.PyConfig.IgnoreCommandLineOptions = True

from JetMET.CMSJMECalculators.CMSJMECalculators import loadJMESystematicsCalculators
from JetMET.CMSJMECalculators.CMSJMECalculators import config as calcConfigs
from JetMET.CMSJMECalculators.CMSJMECalculators.utils import toRVecFloat, toRVecInt, getMETUnclDeltaXY

from .utils import configLogger, deltaR2

logger = logging.getLogger('jme')
configLogger('jme', loglevel=logging.INFO)


def rndSeed(event, jets, extra=0):
    seed = (event.run << 20) + (event.luminosityBlock << 10) + event.event + extra
    if len(jets) > 0:
        seed += int(jets[0].eta / 0.01)
    return seed


def _safeget(obj, name, default):
    '''Object attribute with fallback; NanoAODTools raises RuntimeError for missing branches.'''
    try:
        return getattr(obj, name)
    except RuntimeError:
        return default


# JME correction files from the central analysis metadata repository on cvmfs.
# Each era is pinned to a dated snapshot (the latest one available in the repo)
# rather than the moving 'latest' link: a dated snapshot is frozen, so the exact
# correction version used is fixed by this code and remote updates stay invisible.
# To pick up newer corrections, bump the date to the latest one listed on cvmfs.
_JME_BASE = '/cvmfs/cms-griddata.cern.ch/cat/metadata/JME'
_JME_SNAPSHOT = {
    'Run2-2016preVFP-UL-NanoAODv15': '2026-06-05',
    'Run2-2016postVFP-UL-NanoAODv15': '2026-06-05',
    'Run2-2017-UL-NanoAODv15': '2026-06-05',
    'Run2-2018-UL-NanoAODv15': '2026-06-05',
    'Run3-22CDSep23-Summer22-NanoAODv12': '2026-06-05',
    'Run3-22EFGSep23-Summer22EE-NanoAODv12': '2026-06-05',
    'Run3-23CSep23-Summer23-NanoAODv12': '2026-07-15',
    'Run3-23DSep23-Summer23BPix-NanoAODv12': '2026-07-15',
    'Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15': '2026-07-16',
    'JER-Smearing': '2026-09-18',
}


def _jmeJson(era, filename):
    '''Path of a JME correction file in its pinned dated snapshot of the metadata repository.'''
    return os.path.join(_JME_BASE, era, _JME_SNAPSHOT[era], filename)


class JetMETCorrector(object):
    '''Jet (re-)correction, smearing and Type-1 MET correction, based on CMSJMECalculators
    (JetVariationsCalculator / Type1METVariationsCalculator) and the JME POG json files.

    jec: re-apply jet energy correction (True|False)
    jes: Jet energy scale options
        - None or 'nominal': no variation (update central values only)
        - 'up', 'down': up/down variation, using the total uncertainty
        - combined with jes_source to select an uncertainty source
    jes_source: name of the JES uncertainty source (e.g. 'Absolute_2016'); '' for the total uncertainty
    jes_uncertainty_file_prefix: prefix for grouped sources, '' or None for the full set,
        non-empty (e.g. 'Regrouped_') for the reduced (regrouped) set;
        the prefixes in the JME jsons are 'Regrouped_'
    jer: Jet energy resolution options
        - None: do nothing
        - 'nominal': apply nominal smearing
        - 'up', 'down': up/down variation
    met_unclustered: MET unclustered energy options
        - None: do nothing
        - 'up', 'down': up/down variation of the unclustered energy
    smearMET: propagate the nominal jet smearing to (Type-1 smeared) MET
    applyHEMUnc: apply the 2018 HEM15/16 issue correction (as down-variation of jets and MET)
    '''

    # JEC/JER era configuration: (json directory, MC jecTag, MC jerTag,
    # [(first run, data jecTag), ...] sorted by run range start)
    # The JME releases ship one DATA tag per era (no per-run IOVs).
    _ERAS = {
        2015: ('Run2-2016preVFP-UL-NanoAODv15', 'Summer20UL16APVNanoV15_V1_MC', 'Summer20UL16APV_JRV5_MC', [
            (0, 'Summer20UL16APVNanoV15_V1_DATA'),
        ]),
        2016: ('Run2-2016postVFP-UL-NanoAODv15', 'Summer20UL16NanoV15_V1_MC', 'Summer20UL16_JRV5_MC', [
            (0, 'Summer20UL16NanoV15_V1_DATA'),
        ]),
        2017: ('Run2-2017-UL-NanoAODv15', 'Summer20UL17NanoV15_V1_MC', 'Summer19UL17_JRV4_MC', [
            (0, 'Summer20UL17NanoV15_V1_DATA'),
        ]),
        2018: ('Run2-2018-UL-NanoAODv15', 'Summer20UL18NanoV15_V1_MC', 'Summer19UL18_JRV3_MC', [
            (0, 'Summer20UL18NanoV15_V1_DATA'),
        ]),
        20220: ('Run3-22CDSep23-Summer22-NanoAODv12', 'Summer22_22Sep2023_V4_MC', 'Summer22_22Sep2023_JRV2_MC', [
            (0, 'Summer22_22Sep2023_V4_DATA'),
        ]),
        20221: ('Run3-22EFGSep23-Summer22EE-NanoAODv12', 'Summer22EE_22Sep2023_V4_MC', 'Summer22EE_22Sep2023_JRV2_MC', [
            (0, 'Summer22EE_22Sep2023_V4_DATA'),
        ]),
        20230: ('Run3-23CSep23-Summer23-NanoAODv12', 'Summer23Prompt23_V4_MC', 'Summer23Prompt23_RunCv1234_JRV3_MC', [
            (0, 'Summer23Prompt23_V4_DATA'),
        ]),
        20231: ('Run3-23DSep23-Summer23BPix-NanoAODv12', 'Summer23BPixPrompt23_V4_MC', 'Summer23BPixPrompt23_RunD_JRV3_MC', [
            (0, 'Summer23BPixPrompt23_V4_DATA'),
        ]),
        2024: ('Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15', 'Summer24Prompt24_V5_MC', 'Summer24Prompt24_JRV2_MC', [
            (0, 'Summer24Prompt24_V5_DATA'),
        ]),
    }

    def __init__(
            self, year, jetType="AK4PFchs", jec=False, jes=None, jes_source=None,
            jes_uncertainty_file_prefix=None, jer='nominal', met_unclustered=None,
            smearMET=True, applyHEMUnc=False):

        self.year = year
        self.jetType = jetType
        self.jec = jec
        self.jes = jes
        self.jes_source = '' if jes_source is None else jes_source
        self.jes_uncertainty_file_prefix = '' if jes_uncertainty_file_prefix is None else jes_uncertainty_file_prefix
        self.jer = jer
        self.met_unclustered = met_unclustered
        self.correctMET = jetType.startswith('AK4')
        self.smearMET = smearMET
        self.applyHEMUnc = applyHEMUnc

        # at most one active variation can be propagated (no combined variations)
        active = [x for x in (jes in ['up', 'down'], jer in ['up', 'down'],
                              met_unclustered is not None, applyHEMUnc) if x]
        if len(active) > 1:
            raise RuntimeError('Only one systematic variation can be active at a time: %s' % str(
                dict(jes=jes, jer=jer, met_unclustered=met_unclustered, applyHEMUnc=applyHEMUnc)))

        try:
            self.era, self.globalTag, self.jerTag, self.dataTags = self._ERAS[self.year]
        except KeyError:
            raise RuntimeError('Invalid year: %s' % (str(self.year)))

        # json producer tag and jet algorithm:
        # - AK4: jet_jerc.json.gz; AK8: fatJet_jerc.json.gz
        # the JME releases ship native AK4PFPuppi/AK8PFPuppi payloads,
        # so no cross-algorithm fallback is needed
        self.jetAlgo = jetType
        jsonName = 'fatJet_jerc.json.gz' if 'AK8' in jetType else 'jet_jerc.json.gz'
        self.jsonFile = _jmeJson(self.era, jsonName)
        self.jsonFileSmearingTool = _jmeJson('JER-Smearing', 'jer_smear.json.gz')
        self.genMatchDR = 0.4 if 'AK8' in jetType else 0.2

        # PUPPI jets have a trivial (identity) L1 correction; it is still passed,
        # as Type1METVariationsCalculator requires an L1 JEC to be configured
        self.l1JecTag = 'L1FastJet'

        # active variation and JES uncertainty source name in the json
        self._jetVariation = 'nominal'
        self._metVariation = 'nominal'
        self._jesSourceName = None
        if applyHEMUnc:
            self._jetVariation = self._metVariation = 'jesHEMIssuedown'
        elif jes in ['up', 'down']:
            prefix = 'Regrouped_' if self.jes_uncertainty_file_prefix else ''
            self._jesSourceName = prefix + (self.jes_source if self.jes_source else 'Total')
            self._jetVariation = self._metVariation = 'jes' + self._jesSourceName + jes
        elif jer in ['up', 'down']:
            self._jetVariation = self._metVariation = 'jer' + jer
        if met_unclustered is not None:
            self._metVariation = 'unclustEn' + met_unclustered

        self._jetCalcs = {}      # jecTag -> JetVariationsCalculator
        self._metCalcs = {}      # jecTag -> Type1METVariationsCalculator
        self._seed = 0

    def _config(self, cls, jecTag, jesUncertainties, jerTag, jecLevel):
        cfg = cls(jsonFile=self.jsonFile, jetAlgo=self.jetAlgo)
        cfg.jecTag = jecTag
        cfg.jecLevel = jecLevel
        cfg.jesUncertainties = jesUncertainties
        cfg.jerTag = jerTag
        if jerTag:
            cfg.jsonFileSmearingTool = self.jsonFileSmearingTool
        cfg.genMatchDR = self.genMatchDR
        cfg.addHEM2018Issue = self.applyHEMUnc
        if issubclass(cls, calcConfigs.METVariationsConfigBase):
            cfg.l1JecTag = self.l1JecTag
            cfg.isT1SmearedMET = self.smearMET
        return cfg

    def _variationIndex(self, available, target):
        try:
            return available.index(target)
        except ValueError:
            raise RuntimeError('Variation %s not found among %s' % (target, str(available)))

    def beginJob(self):
        loadJMESystematicsCalculators()

        # MC calculators: reapply JEC if requested, smear if jer is set, JES uncertainties if requested
        jesUncMC = [self._jesSourceName] if self._jesSourceName is not None else []
        jerTagMC = self.jerTag if self.jer is not None else ''
        cfgMC = self._config(calcConfigs.JetVariations, self.globalTag, jesUncMC, jerTagMC,
                             'L1L2L3Res' if self.jec else '')
        self._jetCalcs['MC'] = cfgMC.create()
        self._jetIdxMC = self._variationIndex(list(self._jetCalcs['MC'].available()), self._jetVariation)
        if self.correctMET:
            cfgMETMC = self._config(calcConfigs.METVariations, self.globalTag, jesUncMC, jerTagMC, 'L1L2L3Res')
            self._metCalcs['MC'] = cfgMETMC.create()
            self._metIdxMC = self._variationIndex(list(self._metCalcs['MC'].available()), self._metVariation)

        # data calculators (one per run IOV): reapply JEC; no smearing nor uncertainties
        for _, tag in self.dataTags:
            cfg = self._config(calcConfigs.JetVariations, tag, [], '', 'L1L2L3Res' if self.jec else '')
            self._jetCalcs[tag] = cfg.create()
            if self.correctMET:
                cfgMET = self._config(calcConfigs.METVariations, tag, [], '', 'L1L2L3Res')
                self._metCalcs[tag] = cfgMET.create()
        self._jetIdxDATA = self._metIdxDATA = 0  # nominal only
        logger.info('Init JetMETCorrector: year %s, jetType %s (json %s), jecTag MC=%s DATA=%s, jerTag %s',
                    str(self.year), self.jetType, self.jetAlgo, self.globalTag,
                    ','.join(t for _, t in self.dataTags), jerTagMC or '<none>')

    def endJob(self):
        pass

    def setSeed(self, seed):
        self._seed = int(seed) % (2**31)

    def _dataJecTag(self, runNumber):
        return next(tag for iov, tag in reversed(self.dataTags) if iov <= runNumber)

    def _jetArgs(self, jets, genjets, isMC, rho, forMET=False):
        args = [
            toRVecFloat([j.pt for j in jets]),
            toRVecFloat([j.eta for j in jets]),
            toRVecFloat([j.phi for j in jets]),
            toRVecFloat([j.mass for j in jets]),
            toRVecFloat([_safeget(j, 'rawFactor', 0.) for j in jets]),
            toRVecFloat([_safeget(j, 'area', 0.) for j in jets]),
        ]
        if forMET:
            args += [
                toRVecFloat([_safeget(j, 'muonSubtrFactor', 0.) for j in jets]),
                toRVecFloat([_safeget(j, 'neEmEF', 0.) for j in jets]),
                toRVecFloat([_safeget(j, 'chEmEF', 0.) for j in jets]),
            ]
        args.append(toRVecInt([_safeget(j, 'jetId', 0) for j in jets] if self.applyHEMUnc else []))
        args.append(float(rho))
        if isMC:
            args += [
                toRVecInt([_safeget(j, 'genJetIdx', -1) for j in jets]),
                toRVecInt([_safeget(j, 'partonFlavour', 0) for j in jets]),
                self._seed,
                toRVecFloat([g.pt for g in genjets]),
                toRVecFloat([g.eta for g in genjets]),
                toRVecFloat([g.phi for g in genjets]),
                toRVecFloat([g.mass for g in genjets]),
            ]
        else:
            args += [toRVecInt([]), toRVecInt([]), 0, toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([])]
        return args

    def correctJetAndMET(self, jets, lowPtJets=None, met=None, rawMET=None,
                         rho=None, genjets=None, isMC=True, runNumber=None):
        # Type-1 MET correction: https://github.com/cms-sw/cmssw/blob/master/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h
        if genjets is None:
            genjets = []
        if met is None:
            lowPtJets = []

        metCalc = iVarMET = None
        if isMC:
            jetCalc, iVarJet = self._jetCalcs['MC'], self._jetIdxMC
            if met is not None:
                metCalc, iVarMET = self._metCalcs['MC'], self._metIdxMC
        else:
            tag = self._dataJecTag(runNumber)
            jetCalc, iVarJet = self._jetCalcs[tag], self._jetIdxDATA
            if met is not None:
                metCalc, iVarMET = self._metCalcs[tag], self._metIdxDATA

        # correct (and smear) the jets -- the arguments are built from the *original*
        # (as stored) jet momenta, which the calculators correct from scratch
        # NOTE: the run number is passed to the calculator; only DATA correction
        # payloads with a dedicated 'run' input use it (ignored otherwise)
        runNumber = 0 if runNumber is None else int(runNumber)
        jetArgs = self._jetArgs(jets, genjets, isMC, rho)
        res = jetCalc.produce(*(jetArgs + [runNumber]))
        newPt, newMass = res.pt(iVarJet), res.mass(iVarJet)

        # correct MET
        if met is not None:
            args = jetArgs[:6] + [
                toRVecFloat([_safeget(j, 'muonSubtrFactor', 0.) for j in jets]),
                toRVecFloat([_safeget(j, 'neEmEF', 0.) for j in jets]),
                toRVecFloat([_safeget(j, 'chEmEF', 0.) for j in jets]),
            ] + jetArgs[6:]
            args += [float(rawMET.phi), float(rawMET.pt)]
            if lowPtJets:
                args += [toRVecFloat([_safeget(j, var, 0.) for j in lowPtJets])
                         for var in ('rawPt', 'eta', 'phi', 'area', 'muonSubtrFactor')]
            else:
                args += [toRVecFloat([]) for _ in range(5)]
            args += [toRVecFloat([]), toRVecFloat([])]
            if isMC and self.met_unclustered:
                dx, dy = getMETUnclDeltaXY(float(met.pt), float(met.phi),
                                           _safeget(met, 'ptUnclusteredUp', 0.),
                                           _safeget(met, 'phiUnclusteredUp', 0.))
            else:
                dx, dy = 0., 0.
            args += [dx, dy, runNumber]
            resMET = metCalc.produce(*args)

        # all results are in hand: update the objects
        for i, j in enumerate(jets):
            j.pt, j.mass = float(newPt[i]), float(newMass[i])
        if met is not None:
            met.pt, met.phi = float(resMET.pt(iVarMET)), float(resMET.phi(iVarMET))


class AK15patcher(object):
    '''Fill in the jet composition variables of AK15 (AK15Puppi) jets when they are
    missing from the input NanoAOD.

    The custom NanoAOD inputs store the AK15Puppi collection without the jet
    composition variables (needed e.g. to compute the jetID working points):
        chHEF           charged Hadron Energy Fraction
        neHEF           neutral Hadron Energy Fraction
        chEmEF          charged Electromagnetic Energy Fraction
        neEmEF          neutral Electromagnetic Energy Fraction
        muEF            muon Energy Fraction
        chMultiplicity  (Puppi-weighted) Number of charged particles in the jet
        neMultiplicity  (Puppi-weighted) Number of neutral particles in the jet
    (meanings as documented in the branch titles of the FatJet and Jet collections).

    For each AK15 jet, the missing variables are estimated from all AK8 (FatJet) jets
    geometrically matched within deltaR < drAK8 (default 0.8) of the AK15 axis, plus all
    AK4 (Jet) jets matched within deltaR < drAK15 (default 1.5) of the AK15 axis whose
    axis does not lie within deltaR < drAK8 of an already-matched AK8 jet (those AK4
    jets are already included in the AK8 energy and must not be double-counted).
    The energy fractions are combined as the pt-weighted average of the matched jets,
    the multiplicities as the plain sum -- the aggregation validated against stored
    FatJet values in run/compare_ak8_ak4_poc.py.

    Branches that exist in the input are never modified ("only if missing"); a variable
    is considered missing when reading it raises RuntimeError (NanoAODTools behavior
    for unknown/inactive branches) or AttributeError (plain python objects). Jets with
    no matched AK8/AK4 jet get zero fractions and multiplicities.

    NOTE: this runs before any JEC re-application, so the pt weights of the matched
    AK8/AK4 jets and the AK15 jet are all the stored (uncorrected) values and are
    mutually consistent.
    '''

    _FRACTION_VARS = ('chHEF', 'neHEF', 'chEmEF', 'neEmEF', 'muEF')
    _MULTIPLICITY_VARS = ('chMultiplicity', 'neMultiplicity')

    _MISSING = object()

    def __init__(self, year):
        self._year = year  # currently unused, kept for possible per-year tuning
        self._drAK8 = 0.8    # AK4-AK8 overlap radius
        self._drAK15 = 1.5  # AK4 -> AK15 matching radius

    @staticmethod
    def _get(jet, var, default):
        try:
            return getattr(jet, var)
        except (RuntimeError, AttributeError):
            return default

    def CalculateMissingVars(self, ak4jets, ak8jets, ak15jets):
        '''Compute the missing composition variables of the AK15 jets, in place.

        ak4jets, ak8jets, ak15jets: iterables of jet objects exposing at least
        pt/eta/phi (NanoAODTools Collections; ak4jets is the Jet, ak8jets the
        FatJet and ak15jets the AK15Puppi collection).
        '''
        if len(ak15jets) == 0:
            return
        # branch presence is a property of the collection, so test the first jet only
        missing = [v for v in self._FRACTION_VARS + self._MULTIPLICITY_VARS
                   if self._get(ak15jets[0], v, self._MISSING) is self._MISSING]
        if not missing:
            return
        fracVars = [v for v in missing if v in self._FRACTION_VARS]
        multVars = [v for v in missing if v in self._MULTIPLICITY_VARS]

        match2AK8 = self._drAK8**2
        match2AK15 = self._drAK15**2
        for ak15 in ak15jets:
            eta, phi = ak15.eta, ak15.phi
            matched8 = [j for j in ak8jets if deltaR2(eta, phi, j.eta, j.phi) < match2AK15]
            matched4 = [j for j in ak4jets
                        if deltaR2(eta, phi, j.eta, j.phi) < match2AK15 and
                        all(deltaR2(j8.eta, j8.phi, j.eta, j.phi) >= match2AK8 for j8 in matched8)]
            matched = matched8 + matched4
            ptw = sum(float(j.pt) for j in matched)
            for v in fracVars:
                value = (sum(float(j.pt) * float(self._get(j, v, 0.)) for j in matched) / ptw) \
                    if ptw > 0. else 0.
                setattr(ak15, v, value)
            for v in multVars:
                setattr(ak15, v, int(round(sum(float(self._get(j, v, 0.)) for j in matched))))


class JetIdCalculator(object):
    def __init__(self, year, jetType, nanov15):
        self._year = year
        self._jetType = jetType
        self.era = {2024: 'Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15'}
        self._nanov15 = nanov15

    def beginJob(self):
        if self._nanov15:
            filename = _jmeJson(self.era[self._year], 'jetid.json.gz')
            self.corr = correctionlib.CorrectionSet.from_file(filename)
        else:
            self.corr = None
        
    def CalculateJetId(self, jets):
        for Jet in jets:
            
            if self.corr is None: # jet id is available in nanoaod v12 but it has to be fixed

                if self._year<=2018:
                    setattr(Jet, "TightId", (Jet.jetId & 2) != 0)
                    setattr(Jet, "TightLeptonVetoId", (Jet.jetId & 4) != 0)
                    
                else:
                    tight_bit = (Jet.jetId & 2) != 0
                    pass_tight = tight_bit
                    if abs(Jet.eta) > 2.7 and abs(Jet.eta) <= 3.0:
                        pass_tight = tight_bit and Jet.neHEF < 0.99
                    elif abs(Jet.eta) > 3.0:
                        pass_tight = tight_bit and Jet.neEmEF < 0.4

                    pass_tight_lep_veto = pass_tight
                    if abs(Jet.eta) <= 2.7:
                        pass_tight_lep_veto = (
                            pass_tight
                            and Jet.muEF < 0.8
                            and Jet.chEmEF < 0.8
                        )
                    setattr(event.Jet, "TightId", pass_tight)
                    setattr(event.Jet, "TightLeptonVetoId", pass_tight_lep_veto)
                    Jet.jetId = int(3*(pass_tight) + 4*(pass_tight_lep_veto))
                    
            else: # need to compute jet id from scratch
                variables = {
                    "eta": Jet.eta,
                    "chHEF": Jet.chHEF,
                    "neHEF": Jet.neHEF,
                    "chEmEF": Jet.chEmEF,
                    "neEmEF": Jet.neEmEF,
                    "muEF": Jet.muEF,
                    "chMultiplicity": Jet.chMultiplicity,
                    "neMultiplicity": Jet.neMultiplicity,
                    "multiplicity": Jet.chMultiplicity + Jet.neMultiplicity,
                }
            for WP in ["Tight","TightLeptonVeto"]:
                setattr(Jet, WP+"Id", self.corr[self._jetType+"_"+WP].evaluate(*(variables[input_.name] for input_ in self.corr[self._jetType+"_"+WP].inputs)))
            setattr(Jet, "jetId", int(3*Jet.TightId + 4*Jet.TightLeptonVetoId))

