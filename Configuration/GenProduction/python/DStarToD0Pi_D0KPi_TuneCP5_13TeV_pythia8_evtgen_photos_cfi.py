import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *


generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    maxEventsToPrint=cms.untracked.int32(1),
    pythiaPylistVerbosity=cms.untracked.int32(1),
    filterEfficiency=cms.untracked.double(0.0482),
    pythiaHepMCVerbosity=cms.untracked.bool(False),
    comEnergy=cms.double(13000.0),
    ExternalDecays=cms.PSet(
        EvtGen1=cms.untracked.PSet(
            decay_table=cms.string("GeneratorInterface/EvtGenInterface/data/DECAY_2010.DEC"),
            particle_property_file=cms.FileInPath("GeneratorInterface/EvtGenInterface/data/evt.pdl"),
            user_decay_file=cms.vstring("GeneratorInterface/ExternalDecays/data/Dstar_D0_Kpi.dec"),
            list_forced_decays=cms.vstring("MyD*+", "MyD*-"),
            operates_on_particles=cms.vint32(413, -413),
            convertPythiaCodes=cms.untracked.bool(False),
        ),
        parameterSets=cms.vstring("EvtGen1"),
    ),
    PythiaParameters=cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters=cms.vstring(
            "HardQCD:hardccbar = on",
            "PhaseSpace:pTHatMin = 0.0",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10.0",
        ),
        parameterSets=cms.vstring(
            "pythia8CommonSettings",
            "pythia8CP5Settings",
            "processParameters",
        ),
    ),
)


DstarFilter = cms.EDFilter(
    "PythiaMomDauFilter",
    ParticleID=cms.untracked.int32(413),
    DaughterID=cms.untracked.int32(421),
    ChargeConjugation=cms.untracked.bool(True),
    MinEta=cms.untracked.double(-100.0),
    MaxEta=cms.untracked.double(100.0),
    DaughterIDs=cms.untracked.vint32(421, 211),
    NumberDaughters=cms.untracked.int32(2),
    MomMinPt=cms.untracked.double(3.9),
    NumberDescendants=cms.untracked.int32(2),
    DescendantsIDs=cms.untracked.vint32(-321, 211),
)


ProductionFilterSequence = cms.Sequence(generator * DstarFilter)
