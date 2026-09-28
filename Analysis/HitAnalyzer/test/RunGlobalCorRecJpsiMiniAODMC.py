import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag
from FWCore.ParameterSet.VarParsing import VarParsing


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
options.parseArguments()

output_base = options.outputFile[:-5] if options.outputFile.endswith(".root") else options.outputFile


process = cms.Process("RECO2", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.load("TrackPropagation.Geant4e.geantRefit_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*options.inputFiles),
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

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
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

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerTwoTrackG4e",
    src=cms.InputTag("refittedMuonTracks"),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(True),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(False),
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
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    doVtxConstraint=cms.bool(False),
    doMassConstraint=cms.bool(False),
    massConstraint=cms.double(3.0969),
    massConstraintWidth=cms.double(1e-5),
    corFiles=cms.vstring(),
    triggers=cms.vstring(
        "HLT_Dimuon0_Jpsi_Muon",
        "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
        "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
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
