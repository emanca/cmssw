import FWCore.ParameterSet.Config as cms


# Common two-track configuration used by the J/psi and Upsilon workflows.
# Resonance-specific mass and trigger settings are applied by OniaNtuplizer_cff.
residualGlobalCorrectionMakerTwoTrackG4e = cms.EDProducer(
    "ResidualGlobalCorrectionMakerTwoTrackG4e",
    src=cms.InputTag(""),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(True),
    fillRunTree=cms.bool(True),
    appendStreamIdToOutput=cms.bool(True),
    appendToTFileService=cms.bool(False),
    appendTreeName=cms.string("Events"),
    doGen=cms.bool(False),
    genParticles=cms.InputTag("genParticles"),
    pileupInfo=cms.InputTag("addPileupInfo"),
    requireGen=cms.bool(False),
    doSim=cms.bool(False),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    doMuons=cms.bool(False),
    muons=cms.InputTag("muons"),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(True),
    triggers=cms.vstring(),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(True),
    corFiles=cms.vstring(),
    MagneticFieldLabel=cms.string(""),
    doVtxConstraint=cms.bool(False),
    doMassConstraint=cms.bool(False),
    massConstraint=cms.double(3.0969),
    massConstraintWidth=cms.double(1e-5),
    outprefix=cms.untracked.string("globalcor"),
)
