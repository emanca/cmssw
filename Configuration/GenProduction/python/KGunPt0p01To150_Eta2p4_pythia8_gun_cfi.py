import FWCore.ParameterSet.Config as cms

# KGunPt0p01To150_Eta2p4: single kaon (K-) gun, flat pT 0.01-150.0 GeV/c, |eta| <= 2.4.
# Uses FlatRandomPtGunProducer (the proven central-UL16-gun mechanism,
# per SMP-RunIISummer20UL16GEN-00001), NOT Pythia8.
generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(150.0),
        MinPt = cms.double(0.01),
        PartID = cms.vint32(-321),
        MinEta = cms.double(-2.4),
        MaxEta = cms.double(2.4),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single kaon (K-) pt 0.01 to 150.0'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)

ProductionFilterSequence = cms.Sequence(generator)
