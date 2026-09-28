import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag


options = VarParsing("analysis")
options.parseArguments()

input_files = list(options.inputFiles)
if not input_files:
    input_files = ["file:TkAlDstToD0Pi.root"]

process = cms.Process("D0ALCANTUPLE", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("TrackPropagation.Geant4e.geantRefit_cff")

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

process.TFileService = cms.Service("TFileService", fileName=cms.string("globalcor_d0_alca_mc.root"))

process.GlobalTag = GlobalTag(process.GlobalTag, "106X_mcRun2_asymptotic_v13", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.load("Analysis.HitAnalyzer.DstToD0PiCascadeRefitProducer_cfi")
process.load("Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackKPiG4e_cfi")

process.selectedD0Tracks = process.DstToD0PiCascadeRefitProducer.clone(src="ALCARECOTkAlDstToD0Pi")

process.globalCorD0 = process.globalCorD0.clone(
    src="selectedD0Tracks",
    fitFromGenParms=False,
    fitFromSimParms=False,
    fillTrackTree=True,
    fillGrads=False,
    fillJac=True,
    fillRunTree=True,
    doGen=False,
    doSim=False,
    requireGen=False,
    useIdealGeometry=False,
    respectTrackOrder=True,
    minPairMass=1.65,
    maxPairMass=2.05,
    outprefix="globalcor_d0_alca_mc",
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.p = cms.Path(process.geopro * process.offlineBeamSpot * process.selectedD0Tracks * process.globalCorD0)
process.schedule = cms.Schedule(process.p)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands

process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)
