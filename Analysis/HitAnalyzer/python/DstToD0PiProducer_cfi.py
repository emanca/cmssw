import FWCore.ParameterSet.Config as cms

dstToD0PiProducer = cms.EDProducer(
    'DstToD0PiProducer',
    src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    deferTreeFill=cms.bool(False),
    vertices=cms.InputTag('offlinePrimaryVertices'),
    beamSpot=cms.InputTag('offlineBeamSpot'),
    doTrigger=cms.bool(False),
    triggerBits=cms.InputTag('TriggerResults', '', 'HLT'),
    hltPaths=cms.vstring(),
    verbose=cms.int32(0),
    hadPtMin=cms.double(0.6),
    softPionPtMin=cms.double(-1.0),
    pvalChi2VtxMin=cms.double(0.005),
    dmD0Max=cms.double(0.05),
    dmDstMax=cms.double(0.15),
    dmDstMinusD0Max=cms.double(0.003),
)
