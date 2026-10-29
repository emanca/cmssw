import FWCore.ParameterSet.Config as cms

# Muons with all track references resolvable inside an ALCARECO file, see
# plugins/AlCaRecoMuonRelinker.cc. selectedTracks must be set per stream.
AlCaRecoMuonRelinker = cms.EDProducer('AlCaRecoMuonRelinker',
    muons          = cms.InputTag('muons'),
    selectedTracks = cms.InputTag(''),
    sourceTracks   = cms.InputTag('generalTracks'),
    onlyLinked     = cms.bool(True),
)

