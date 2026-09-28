import FWCore.ParameterSet.Config as cms


dstToD0PiMCTruthProducer = cms.EDProducer(
    'DstToD0PiMCTruthProducer',
    src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    vertices=cms.InputTag('offlinePrimaryVertices'),
    beamSpot=cms.InputTag('offlineBeamSpot'),
    genParticles=cms.InputTag('genParticles'),
    maxDeltaR=cms.double(0.03),
    hadPtMin=cms.double(0.6),
    softPionPtMin=cms.double(-1.0),
    pvalChi2VtxMin=cms.double(0.005),
    dmD0Max=cms.double(0.05),
    dmDstMax=cms.double(0.15),
    dmDstMinusD0Max=cms.double(0.003),
    fiducialEtaMax=cms.double(2.5),
)
