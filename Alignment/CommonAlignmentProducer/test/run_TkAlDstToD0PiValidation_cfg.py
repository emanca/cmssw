import FWCore.ParameterSet.Config as cms

process = cms.Process('DSTVAL')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(
        'file:TkAlDstToD0Pi.root'
    )
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
)

process.load('Alignment.CommonAlignmentProducer.TkAlDstToD0PiValidation_cff')

process.TFileService = cms.Service(
    'TFileService',
    fileName = cms.string('TkAlDstToD0PiValidation.root')
)

process.p = cms.Path(process.seqTkAlDstToD0PiValidation)
