import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
ALCARECOTkAlDstToD0PiHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    andOr = True,
    HLTPaths = cms.vstring(
        'HLT_IsoMu24_v*',
        'HLT_IsoTkMu24_v*',
    ),
    throw = False
)

import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlDstToD0PiDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter = cms.bool(True),
    AndOr = cms.bool(True),
    DebugOn = cms.untracked.bool(False)
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi
ALCARECOTkAlDstToD0Pi = Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi.AlignmentTrackSelector.clone()
ALCARECOTkAlDstToD0Pi.filter = True
ALCARECOTkAlDstToD0Pi.src = 'generalTracks'

ALCARECOTkAlDstToD0Pi.applyBasicCuts = True
ALCARECOTkAlDstToD0Pi.ptMin = 0.5
ALCARECOTkAlDstToD0Pi.etaMin = -3.5
ALCARECOTkAlDstToD0Pi.etaMax = 3.5
ALCARECOTkAlDstToD0Pi.nHitMin = 0

ALCARECOTkAlDstToD0Pi.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlDstToD0Pi.GlobalSelector.applyIsolationtest = False

ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyMassrangeFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.minXMass = 1.860
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.maxXMass = 2.160
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyIntermediateMassrangeFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.minIntermediateMass = 1.81483
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.maxIntermediateMass = 1.91483
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyMassDifferenceFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.minMassDifference = 0.14243
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.maxMassDifference = 0.14843
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.firstDaughterMass = 0.493677
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.secondDaughterMass = 0.139570
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.thirdDaughterMass = 0.139570
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyChargeFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.charge = 1
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.useUnsignedCharge = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.numberOfCandidates = 1

seqALCARECOTkAlDstToD0Pi = cms.Sequence(ALCARECOTkAlDstToD0PiHLT+
                                        ALCARECOTkAlDstToD0PiDCSFilter+
                                        ALCARECOTkAlDstToD0Pi)
