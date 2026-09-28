# HitAnalyzer ntuplizers

This package contains the track-level ntuplizers used for the onia and D0
studies. The physics producers are intentionally separate; the Python
configuration is shared wherever the input format and reconstruction are the
same.

## Supported workflows

### J/psi and Upsilon

`ResidualGlobalCorrectionMakerTwoTrackG4e` is the common two-track producer
for both resonances. The ALCARECO source, mass, and trigger list are selected
in `python/OniaNtuplizer_cff.py`; common producer defaults are in
`python/ResidualGlobalCorrectionMakerTwoTrackG4e_cfi.py`.

ALCARECO entry points:

```text
test/RunGlobalCorRecJpsiAlca.py
test/RunGlobalCorRecJpsiData.py
test/RunGlobalCorRecUpsilonAlca.py
test/RunGlobalCorRecUpsilonData.py
```

These four wrappers read ALCARECO collections. They share the `mode=jpsi` or
`mode=upsilon`, `isMC`, and `do3DFieldMap` options; the J/psi and Upsilon files
provide the corresponding defaults. Use `isMC=0` for data and
`do3DFieldMap=1` when the 3D field map is needed.

The separate MiniAOD path is J/psi-specific:
`test/RunGlobalCorRecJpsiMiniAODMC.py`. It builds tracks from the inner tracks
of `slimmedMuons`, refits them, and then runs the same two-track producer. Its
configuration exposes `ptMin`; it does not select Upsilon through `mode`.

### D0 -> K pi and D* -> D0 pi

The D entry points provide two ALCARECO selection paths:

1. `runDstMiniAODTransientEmulationFromAlca.py` runs `DstToD0PiProducer` on
   `ALCARECOTkAlDstToD0Pi`, then passes its `dstarTracks` to the soft-pion
   CVH/G4e refit and K/pi ntuplizer. With `isMC=1`, it also adds truth paths.
2. `runGlobalCorRecD0FromSelectedAlca.py` and
   `RunGlobalCorRecD0AlcaMC.py` run `DstToD0PiCascadeRefitProducer` on the
   same ALCARECO collection and pass selected `dstarTracks` downstream. The
   former enables soft-pion CVH inputs; the latter is the MC-tagged cascade
   variant without that CVH module.

`ResidualGlobalCorrectionMakerTwoTrackKPiG4e` supplies the K/pi ntuplizer.
The runnable D configurations enable Jacobians, not gradients; soft-pion CVH
is enabled only by the transient-emulation and selected-ALCARECO paths.

Main entry points:

```text
test/runDstMiniAODTransientEmulationFromAlca.py
test/runGlobalCorRecD0FromSelectedAlca.py
test/RunGlobalCorRecD0AlcaMC.py
```

## Build and configuration checks

From the CMSSW project area containing `src/`:

```bash
cd /path/to/CMSSW_10_6_26
scram b -j8
cd src/Analysis/HitAnalyzer/test
cmsRun RunGlobalCorRecJpsiAlca.py maxEvents=0 isMC=0
cmsRun RunGlobalCorRecUpsilonAlca.py maxEvents=0 isMC=0
cmsRun RunGlobalCorRecJpsiMiniAODMC.py maxEvents=0
cmsRun runDstMiniAODTransientEmulationFromAlca.py maxEvents=0 inputFiles=file:TkAlDstToD0Pi.root
cmsRun runGlobalCorRecD0FromSelectedAlca.py maxEvents=0 inputFiles=file:TkAlDstToD0Pi.root
cmsRun RunGlobalCorRecD0AlcaMC.py maxEvents=0 inputFiles=file:TkAlDstToD0Pi.root
```

The D commands require a suitable `TkAlDstToD0Pi.root` input; the MiniAOD
J/psi command uses its configured default unless `inputFiles` is supplied.

Do not commit ROOT outputs, CRAB project directories, logs, or generated
Python cache files.
