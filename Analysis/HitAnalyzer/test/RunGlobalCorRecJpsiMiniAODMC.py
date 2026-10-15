import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.AlCa.GlobalTag import GlobalTag
from FWCore.ParameterSet.VarParsing import VarParsing
from Analysis.HitAnalyzer.OniaNtuplizer_cff import onia_settings


options = VarParsing("analysis")
options.outputFile = "globalcor"
options.maxEvents = -1
options.register(
    "ptMin",
    3.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum slimmedMuon pT used to build the refitted track collection.",
)
options.register(
    "year",
    2016,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Data-taking year used to select the CMSSW era and MC conditions (2016 or 2018).",
)
options.register(
    "doGen",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Fill the gen-level branches of the output tree (required for the muon scale "
    "response systematics and MC-truth checks). Default on for MC; pass "
    "doGen=0 only when rerunning legacy productions that must stay identical.",
)
options.parseArguments()

onia = onia_settings("jpsi")
if options.year not in (2016, 2018):
    raise RuntimeError("Unsupported year {}; use year=2016 or year=2018.".format(options.year))

year_configuration = {
    2016: {
        "era": Run2_2016,
        "global_tag": "106X_mcRun2_asymptotic_v17",
        "geometry_tag": "XMLFILE_Geometry_2016_81YV1_Extended2016_mc",
    },
    2018: {
        "era": Run2_2018,
        "global_tag": "102X_upgrade2018_realistic_v12",
        "geometry_tag": None,
    },
}[options.year]

input_files = list(options.inputFiles)
if not input_files:
    input_files = [
        {
            2016: (
                "/store/mc/RunIISummer20UL16MiniAODv2/"
                "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/"
                "MINIAODSIM/BPH_106X_mcRun2_asymptotic_v17-v2/60000/"
                "CA5DA24E-EB62-6948-B15B-CE63005A8A95.root"
            ),
            2018: (
                "/store/mc/RunIIAutumn18MiniAOD/"
                "BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen/"
                "MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/"
                "E4F2E0BE-896D-2843-B56D-808688E567B7.root"
            ),
        }[options.year]
    ]

output_base = options.outputFile[:-5] if options.outputFile.endswith(".root") else options.outputFile


process = cms.Process("RECO2", year_configuration["era"])

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.load("TrackPropagation.Geant4e.geantRefit_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackG4e_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

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
    annotation=cms.untracked.string("globalCor JPsi displaced ntuples from BPH MINIAODSIM"),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)

process.GlobalTag = GlobalTag(process.GlobalTag, year_configuration["global_tag"], "")
if year_configuration["geometry_tag"]:
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(
            record=cms.string("GeometryFileRcd"),
            tag=cms.string(year_configuration["geometry_tag"]),
            label=cms.untracked.string("Extended"),
        ),
    )
    process.XMLFromDBSource.label = cms.string("Extended")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.load("Analysis.HitAnalyzer.TrackProducerFromPatMuons_cfi")
process.tracksFromPatMuons.src = cms.InputTag("slimmedMuons")
process.tracksFromPatMuons.innerTrackOnly = cms.bool(True)
process.tracksFromPatMuons.ptMin = cms.double(options.ptMin)

import RecoTracker.TrackProducer.TrackRefitter_cfi

process.refittedMuonTracks = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone(
    src="tracksFromPatMuons",
    NavigationSchool="",
    Fitter="FlexibleKFFittingSmoother",
)

process.globalCor = process.residualGlobalCorrectionMakerTwoTrackG4e.clone(
    src=cms.InputTag("refittedMuonTracks"),
    doGen=cms.bool(options.doGen),
    genParticles=cms.InputTag("prunedGenParticles"),
    pileupInfo=cms.InputTag("slimmedAddPileupInfo"),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(True),
    muons=cms.InputTag("slimmedMuons"),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(True),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(True),
    doVtxConstraint=cms.bool(False),
    doMassConstraint=cms.bool(False),
    massConstraint=cms.double(onia["mass"]),
    massConstraintWidth=cms.double(1e-5),
    triggers=cms.vstring(*onia["triggers"]),
    MagneticFieldLabel=cms.string(""),
    outprefix=cms.untracked.string(output_base),
)

process.reconstruction_step = cms.Path(
    process.geopro
    * process.offlineBeamSpot
    * process.tracksFromPatMuons
    * process.refittedMuonTracks
    * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)
