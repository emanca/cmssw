import os

import FWCore.ParameterSet.Config as cms


_cfg_dir = os.getcwd()
if "__file__" in globals():
    _cfg_dir = os.path.dirname(os.path.abspath(__file__))

_base_cfg = os.path.join(_cfg_dir, "DStarToD0Pi_D0KPi_TuneCP5_RECOALCA_cfg.py")
if not os.path.exists(_base_cfg):
    _base_cfg = "DStarToD0Pi_D0KPi_TuneCP5_RECOALCA_cfg.py"

execfile(_base_cfg)

process.maxEvents.input = cms.untracked.int32(-1)

# For CRAB publication keep only the target ALCARECO stream. The validated local
# cfg also writes a GEN-SIM-RECO side output, which is useful for debugging but
# unnecessarily large for the production target.
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.recosim_step,
    process.pathALCARECOTkAlDstToD0Pi,
    process.endjob_step,
    process.ALCARECOStreamTkAlDstToD0PiOutPath,
)
