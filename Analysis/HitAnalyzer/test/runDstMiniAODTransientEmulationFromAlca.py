import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag
import os
import subprocess

options = VarParsing('analysis')
options.register(
    'inputDir',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'Directory containing TkAlDstToD0Pi ROOT files. Supports local, /eos, and /store paths.',
)
options.register(
    'isMC',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Use MC conditions when set to 1.',
)
options.register(
    'transientHadPtMin',
    0.6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    'Transient producer hadron pT threshold in GeV.',
)
options.register(
    'transientSoftPionPtMin',
    -1.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    'Transient producer soft-pion pT threshold in GeV; negative keeps the common hadron threshold.',
)
options.register(
    'transientPvalMin',
    0.005,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    'Transient producer vertex-fit p-value threshold.',
)
options.register(
    'transientD0MassMax',
    0.05,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    'Transient fitted D0 mass half-window in GeV.',
)
options.register(
    'transientDstMassMax',
    0.15,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    'Transient fitted D* mass half-window in GeV.',
)
options.register(
    'transientDeltaMMax',
    0.003,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    'Transient fitted DeltaM half-window in GeV.',
)
options.parseArguments()


def file_name_for_pool_source(path):
    if path.startswith('/eos/uscms/'):
        return 'root://cmseos.fnal.gov//' + path[len('/eos/uscms/'):]
    if path.startswith('/eos/'):
        return 'root://eoscms.cern.ch//' + path
    if path.startswith('/store/'):
        return 'root://cms-xrd-global.cern.ch//' + path
    if path.startswith('root://') or path.startswith('file:'):
        return path
    return 'file:' + path


def list_xrdfs_root_files(server, remote_path):
    try:
        listing = subprocess.check_output(['xrdfs', server, 'ls', '-R', remote_path])
    except Exception:
        listing = subprocess.check_output(['xrdfs', server, 'ls', remote_path])
    if not isinstance(listing, str):
        listing = listing.decode('utf-8')
    return [
        server + '/' + name
        for name in sorted(listing.splitlines())
        if name.endswith('.root')
    ]


def list_root_files(input_dir):
    input_dir = input_dir.rstrip('/')
    if os.path.isdir(input_dir):
        return [
            file_name_for_pool_source(os.path.join(input_dir, name))
            for name in sorted(os.listdir(input_dir))
            if name.endswith('.root')
        ]

    if input_dir.startswith('/eos/uscms/'):
        eos_path = input_dir[len('/eos/uscms'):]
        return list_xrdfs_root_files('root://cmseos.fnal.gov', eos_path)

    if input_dir.startswith('/eos/'):
        return list_xrdfs_root_files('root://eoscms.cern.ch', input_dir)

    if input_dir.startswith('/store/'):
        return list_xrdfs_root_files('root://eoscms.cern.ch', '/eos/cms' + input_dir)

    raise RuntimeError('Input directory does not exist locally and no remote listing rule is known: {}'.format(input_dir))


input_files = list(options.inputFiles)
if options.inputDir:
    input_files = list_root_files(options.inputDir)

if not input_files:
    input_files = ['file:TkAlDstToD0Pi.root']

process = cms.Process('DSTTTEMU', eras.Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('TrackPropagation.Geant4e.geantRefit_cff')
process.load('Analysis.HitAnalyzer.DstToD0PiProducer_cfi')
process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackKPiG4e_cfi')
if options.isMC:
    process.load('Analysis.HitAnalyzer.DstToD0PiMCTruthProducer_cfi')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(*input_files)
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)

process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string(options.outputFile),
)

process.GlobalTag = GlobalTag(
    process.GlobalTag,
    '106X_mcRun2_asymptotic_v13' if options.isMC else 'auto:run2_data',
    '',
)

process.dstToD0PiProducer = process.dstToD0PiProducer.clone(
    src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    vertices=cms.InputTag('offlinePrimaryVertices'),
    beamSpot=cms.InputTag('offlineBeamSpot'),
    hadPtMin=cms.double(options.transientHadPtMin),
    softPionPtMin=cms.double(options.transientSoftPionPtMin),
    pvalChi2VtxMin=cms.double(options.transientPvalMin),
    dmD0Max=cms.double(options.transientD0MassMax),
    dmDstMax=cms.double(options.transientDstMassMax),
    dmDstMinusD0Max=cms.double(options.transientDeltaMMax),
    deferTreeFill=cms.bool(True),
)

process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

# The transient producer above remains the authoritative selector and keeps
# writing its original tree.  Its selected best triplet is also passed to the
# existing single-track and double-track CVH/G4e refit modules below.
process.softPionCor = cms.EDProducer(
    'ResidualGlobalCorrectionMakerG4e',
    src=cms.InputTag('dstToD0PiProducer', 'dstarTracks'),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(False),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(False),
    fillRunTree=cms.bool(False),
    doGen=cms.bool(False),
    genParticles=cms.InputTag('genParticles'),
    pileupInfo=cms.InputTag('addPileupInfo'),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(False),
    muons=cms.InputTag('muons'),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    triggers=cms.vstring(),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(False),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    corFiles=cms.vstring(),
    MagneticFieldLabel=cms.string(''),
    outprefix=cms.untracked.string('soft_pion_cvh')
)

process.globalCorD0 = process.globalCorD0.clone(
    src=cms.InputTag('dstToD0PiProducer', 'dstarTracks'),
    useIdealGeometry=cms.bool(False),
    fillRunTree=cms.bool(True),
    fillJac=cms.bool(True),
    respectTrackOrder=cms.bool(True),
    useDstarTriplets=cms.bool(True),
    enforceUniqueTracks=cms.bool(True),
    writeLegacyBranches=cms.bool(False),
    appendStreamIdToOutput=cms.bool(False),
    appendToTFileService=cms.bool(True),
    appendTreeName=cms.string('dstToD0PiProducer/Events'),
    useSoftPionCvh=cms.bool(True),
    saveGenParticles=cms.bool(bool(options.isMC)),
    # Keep the transient-selected entry when the optional CVH/G4e refit fails.
    fillTruthMatchOnCvhFailure=cms.bool(True),
    genParticles=cms.InputTag('genParticles'),
    softPionCvhPt=cms.InputTag('softPionCor', 'corPt'),
    softPionCvhEta=cms.InputTag('softPionCor', 'corEta'),
    softPionCvhPhi=cms.InputTag('softPionCor', 'corPhi'),
    softPionCvhCharge=cms.InputTag('softPionCor', 'corCharge'),
    softPionCvhEdmval=cms.InputTag('softPionCor', 'edmval'),
    softPionCvhMomCov=cms.InputTag('softPionCor', 'momCov'),
    # The triplet has already passed the transient D0 fit.  Do not apply a
    # second independent raw K/pi mass selection in the CVH module.
    minPairMass=cms.double(-1.0),
    maxPairMass=cms.double(-1.0),
    outprefix=cms.untracked.string('globalcor_d0_cvh')
)

if options.isMC:
    process.dstToD0PiMCTruthProducer = process.dstToD0PiMCTruthProducer.clone(
        src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
        vertices=cms.InputTag('offlinePrimaryVertices'),
        beamSpot=cms.InputTag('offlineBeamSpot'),
        genParticles=cms.InputTag('genParticles'),
        hadPtMin=cms.double(options.transientHadPtMin),
        softPionPtMin=cms.double(options.transientSoftPionPtMin),
        pvalChi2VtxMin=cms.double(options.transientPvalMin),
        dmD0Max=cms.double(options.transientD0MassMax),
        dmDstMax=cms.double(options.transientDstMassMax),
        dmDstMinusD0Max=cms.double(options.transientDeltaMMax),
    )

    process.truthSoftPionCor = process.softPionCor.clone(
        src=cms.InputTag('dstToD0PiMCTruthProducer', 'truthMatchedDstarTracks'),
        outprefix=cms.untracked.string('truth_soft_pion_cvh'),
    )

    process.truthGlobalCorD0 = process.globalCorD0.clone(
        src=cms.InputTag('dstToD0PiMCTruthProducer', 'truthMatchedDstarTracks'),
        fillRunTree=cms.bool(False),
        appendTreeName=cms.string('dstToD0PiMCTruthProducer/MatchedCandidates'),
        saveGenParticles=cms.bool(False),
        fillTruthMatchOnCvhFailure=cms.bool(True),
        softPionCvhPt=cms.InputTag('truthSoftPionCor', 'corPt'),
        softPionCvhEta=cms.InputTag('truthSoftPionCor', 'corEta'),
        softPionCvhPhi=cms.InputTag('truthSoftPionCor', 'corPhi'),
        softPionCvhCharge=cms.InputTag('truthSoftPionCor', 'corCharge'),
        softPionCvhEdmval=cms.InputTag('truthSoftPionCor', 'edmval'),
        softPionCvhMomCov=cms.InputTag('truthSoftPionCor', 'momCov'),
        outprefix=cms.untracked.string('truth_globalcor_d0_cvh'),
    )

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string('GeometryFileRcd'),
        tag=cms.string('XMLFILE_Geometry_2016_81YV1_Extended2016_mc'),
        label=cms.untracked.string('Extended'),
    ),
)
process.XMLFromDBSource.label = cms.string('Extended')

process.p = cms.Path(
    process.geopro *
    process.offlineBeamSpot *
    process.dstToD0PiProducer *
    process.softPionCor *
    process.globalCorD0
)

if options.isMC:
    process.p += (
        process.dstToD0PiMCTruthProducer *
        process.truthSoftPionCor *
        process.truthGlobalCorD0
    )
