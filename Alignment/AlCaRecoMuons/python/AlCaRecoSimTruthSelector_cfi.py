import FWCore.ParameterSet.Config as cms

# Generator-linked SimTracks of the given |pdgId| and their tracker PSimHits,
# see plugins/AlCaRecoSimTruthSelector.cc. Meant to run with module label
# g4SimHits in the ALCA process, so that consumers reading g4SimHits (e.g.
# ResidualGlobalCorrectionMaker*) see the selected products once the SIM-process
# originals are dropped from the output.
_trackerHits = ['TrackerHitsPixelBarrelLowTof', 'TrackerHitsPixelBarrelHighTof',
                'TrackerHitsPixelEndcapLowTof', 'TrackerHitsPixelEndcapHighTof',
                'TrackerHitsTIBLowTof', 'TrackerHitsTIBHighTof',
                'TrackerHitsTIDLowTof', 'TrackerHitsTIDHighTof',
                'TrackerHitsTOBLowTof', 'TrackerHitsTOBHighTof',
                'TrackerHitsTECLowTof', 'TrackerHitsTECHighTof']

AlCaRecoSimTruthSelector = cms.EDProducer('AlCaRecoSimTruthSelector',
    simTracks = cms.InputTag('g4SimHits', '', 'SIM'),
    simHits   = cms.VInputTag([cms.InputTag('g4SimHits', h, 'SIM') for h in _trackerHits]),
    pdgIds    = cms.vint32(13),
)
