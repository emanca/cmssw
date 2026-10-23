#!/usr/bin/env python3
"""Generate the three UL16 gun fragments (mu, pi, K), flat pT 1-150 GeV/c.

Model: /DoubleMuonGun_Pt3To150/RunIISummer20UL16RECO-NoPU_106X_mcRun2_asymptotic_v13-v2
chain conventions (central UL16 gun, SMP-RunIISummer20UL16RECO-00001 family),
with pT range widened to 1-150 as requested and single-particle guns per species.
Anti-particle production is enabled via PGunParameters/AddAntiParticle=True.

Output: three *_cfi.py files next to this script.
"""
import os

TEMPLATE = """import FWCore.ParameterSet.Config as cms

# {label}: single {name} gun, flat pT {ptmin}-{ptmax} GeV/c, |eta| <= {etamax}.
# FlatRandomPtGunProducer (proven central-UL16-gun mechanism,
# SMP-RunIISummer20UL16GEN-00001), not Pythia8.
generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double({PT}),
        MinPt = cms.double({pt}),
        PartID = cms.vint32({pdgid}),
        MinEta = cms.double(-{eta}),
        MaxEta = cms.double({eta}),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single {name} pt {ptmin} to {ptmax}'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)

ProductionFilterSequence = cms.Sequence(generator)
"""

SPECIES = {
    "MuGunPt0p01To150_Eta2p4": dict(name="muon (mu-)", label="MuGun", pdgid=13),
    "PiGunPt0p01To150_Eta2p4": dict(name="pion (pi-)", label="PiGun", pdgid=211),
    "KGunPt0p01To150_Eta2p4": dict(name="kaon (K-)", label="KGun", pdgid=321),
}

PTMIN, PTMAX, ETAMAX = 0.01, 150.0, 2.4

KIN = """MinP = cms.double(0.0),
        MaxP = cms.double(999999.0),
        MinPt = cms.double({pt}),
        MaxPt = cms.double({PT}),
        MinEta = cms.double(-{eta}),
        MaxEta = cms.double({eta}),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359)"""

here = os.path.dirname(os.path.abspath(__file__))
for fname, spec in SPECIES.items():
    kin_block = KIN.format(pt=PTMIN, PT=PTMAX, eta=ETAMAX)
    text = TEMPLATE.format(
        label=fname, name=spec["name"], pdgid=spec["pdgid"],
        ptmin=PTMIN, ptmax=PTMAX, etamax=ETAMAX,
        kin=kin_block,
    )
    out = os.path.join(here, fname + "_pythia8_gun_cfi.py")
    with open(out, "w") as fh:
        fh.write(text)
    print("wrote", out)
