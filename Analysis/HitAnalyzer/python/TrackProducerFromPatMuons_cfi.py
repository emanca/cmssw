import FWCore.ParameterSet.Config as cms

tracksFromPatMuons = cms.EDProducer(
    "TrackProducerFromPatMuons",
    src=cms.InputTag("slimmedMuons"),
    innerTrackOnly=cms.bool(True),
    ptMin=cms.double(0.0),
)
