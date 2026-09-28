import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from FWCore.ParameterSet.VarParsing import VarParsing

from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackG4e_cfi import (
    residualGlobalCorrectionMakerTwoTrackG4e,
)


_ONIA_SETTINGS = {
    "jpsi": {
        "source": "ALCARECOTkAlJpsiMuMu",
        "mass": 3.0969,
        "mc_input": "/store/mc/RunIISummer20UL16RECO/JPsiToMuMu_Pt8toInf-pythia8/ALCARECO/TkAlJpsiMuMu-106X_mcRun2_asymptotic_v13-v2/2820000/02536347-8E42-C74C-9ED0-7EAF4290973A.root",
        "data_input": "/store/data/Run2016H/Charmonium/ALCARECO/TkAlJpsiMuMu-21Feb2020_UL2016-v1/240000/326FE464-5AF9-4C48-B0B7-1172297E9B71.root",
        "triggers": (
            "HLT_Dimuon0_Jpsi_Muon",
            "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
            "HLT_Dimuon0er16_Jpsi_NoVertexing",
            "HLT_Dimuon10_Jpsi_Barrel",
            "HLT_Dimuon13_PsiPrime",
            "HLT_Dimuon16_Jpsi",
            "HLT_Dimuon20_Jpsi",
            "HLT_Dimuon8_PsiPrime_Barrel",
            "HLT_DoubleMu4_3_Bs",
            "HLT_DoubleMu4_3_Jpsi_Displaced",
            "HLT_DoubleMu4_JpsiTrk_Displaced",
            "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
            "HLT_Mu7p5_Track2_Jpsi",
            "HLT_Mu7p5_Track3p5_Jpsi",
            "HLT_Mu7p5_Track7_Jpsi",
        ),
    },
    "upsilon": {
        "source": "ALCARECOTkAlUpsilonMuMu",
        "mass": 9.4603,
        "mc_input": "/store/mc/RunIISummer20UL16RECO/Upsilon_MultiResonanceFilter_Pt8to13-pythia8/ALCARECO/TkAlUpsilonMuMu-106X_mcRun2_asymptotic_v13-v2/30000/0387C9FE-F2A1-D14A-9064-72B99FF3B170.root",
        "data_input": "/store/data/Run2016G/MuOnia/ALCARECO/TkAlUpsilonMuMu-31Aug2023_UL2016_WMass-v1/2540000/088954B9-4392-3D46-8D78-6C9B4CFCB8A2.root",
        "triggers": (
            "HLT_Dimuon0_Phi_Barrel",
            "HLT_Dimuon0_Upsilon_Muon",
            "HLT_Dimuon13_Upsilon",
            "HLT_Dimuon8_Upsilon_Barrel",
            "HLT_Mu16_TkMu0_dEta18_Onia",
            "HLT_Mu16_TkMu0_dEta18_Phi",
            "HLT_Mu25_TkMu0_dEta18_Onia",
            "HLT_Mu7p5_L2Mu2_Upsilon",
            "HLT_Mu7p5_Track2_Upsilon",
            "HLT_Mu7p5_Track3p5_Upsilon",
            "HLT_Mu7p5_Track7_Upsilon",
            "HLT_QuadMuon0_Dimuon0_Upsilon",
        ),
    },
}


def onia_settings(mode):
    key = mode.lower()
    if key not in _ONIA_SETTINGS:
        raise RuntimeError("Unsupported onia mode '{}'; use 'jpsi' or 'upsilon'.".format(mode))
    return _ONIA_SETTINGS[key]


def make_onia_producer(src, mode, is_mc, output_file="globalcor", use_ideal_geometry=None):
    settings = onia_settings(mode)
    if use_ideal_geometry is None:
        use_ideal_geometry = is_mc

    output_prefix = output_file[:-5] if output_file.endswith(".root") else output_file
    return residualGlobalCorrectionMakerTwoTrackG4e.clone(
        src=cms.InputTag(src),
        massConstraint=cms.double(settings["mass"]),
        triggers=cms.vstring(*settings["triggers"]),
        useIdealGeometry=cms.bool(use_ideal_geometry),
        outprefix=cms.untracked.string(output_prefix),
    )


def make_onia_process(default_mode="jpsi", default_is_mc=True, default_do_3d_field_map=False):
    options = VarParsing("analysis")
    options.outputFile = "globalcor"
    options.maxEvents = -1
    options.register(
        "mode",
        default_mode,
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "Resonance to process: jpsi or upsilon.",
    )
    options.register(
        "isMC",
        int(default_is_mc),
        VarParsing.multiplicity.singleton,
        VarParsing.varType.int,
        "Use MC conditions and the MC default input.",
    )
    options.register(
        "do3DFieldMap",
        int(default_do_3d_field_map),
        VarParsing.multiplicity.singleton,
        VarParsing.varType.int,
        "Use the 3D polynomial magnetic-field map.",
    )
    options.parseArguments()

    mode = options.mode.lower()
    is_mc = bool(options.isMC)
    settings = onia_settings(mode)
    input_files = list(options.inputFiles)
    if not input_files:
        input_files = [settings["mc_input"] if is_mc else settings["data_input"]]

    process = cms.Process("RECO2", Run2_2016)
    process.load("Configuration.StandardSequences.Services_cff")
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.load("Configuration.EventContent.EventContent_cff")
    process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.load("Configuration.StandardSequences.EndOfProcess_cff")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.load("Configuration.StandardSequences.GeometrySimDB_cff")
    process.load("TrackPropagation.Geant4e.geantRefit_cff")

    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
    process.source = cms.Source(
        "PoolSource",
        fileNames=cms.untracked.vstring(*input_files),
        secondaryFileNames=cms.untracked.vstring(),
    )
    process.options = cms.untracked.PSet(
        numberOfThreads=cms.untracked.uint32(1),
        numberOfStreams=cms.untracked.uint32(1),
        numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
    )
    process.configurationMetadata = cms.untracked.PSet(
        annotation=cms.untracked.string("HitAnalyzer onia ntuplizer"),
        name=cms.untracked.string("Applications"),
        version=cms.untracked.string("$Revision: 1.19 $"),
    )

    process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc" if is_mc else "auto:run2_data", "")
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(
            record=cms.string("GeometryFileRcd"),
            tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
            label=cms.untracked.string("Extended"),
        )
    )
    process.XMLFromDBSource.label = cms.string("Extended")
    process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")
    process.globalCor = make_onia_producer(
        settings["source"],
        mode,
        is_mc,
        output_file=options.outputFile,
    )

    if bool(options.do3DFieldMap):
        from MagneticField.ParametrizedEngine.parametrizedMagneticField_PolyFit3D_cfi import (
            ParametrizedMagneticFieldProducer,
        )

        process.PolyFit3DMagneticFieldProducer = ParametrizedMagneticFieldProducer
        field_label = "PolyFit3DMf"
        process.PolyFit3DMagneticFieldProducer.label = field_label
        process.geopro.MagneticFieldLabel = field_label
        process.Geant4ePropagator.MagneticFieldLabel = field_label
        process.stripCPEESProducer.MagneticFieldLabel = field_label
        process.StripCPEfromTrackAngleESProducer.MagneticFieldLabel = field_label
        process.siPixelTemplateDBObjectESProducer.MagneticFieldLabel = field_label
        process.templates.MagneticFieldLabel = field_label
        process.TransientTrackBuilderESProducer.MagneticFieldLabel = field_label
        process.globalCor.MagneticFieldLabel = field_label

    process.reconstruction_step = cms.Path(process.geopro * process.offlineBeamSpot * process.globalCor)
    process.schedule = cms.Schedule(process.reconstruction_step)

    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

    associatePatAlgosToolsTask(process)

    from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands

    process = customiseLogErrorHarvesterUsingOutputCommands(process)

    from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

    process = customiseEarlyDelete(process)
    return process
