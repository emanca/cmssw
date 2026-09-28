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
        try:
            return list_xrdfs_root_files('root://cmseos.fnal.gov', eos_path)
        except Exception as exc:
            raise RuntimeError('Could not list {} locally or via xrdfs: {}'.format(input_dir, exc))

    if input_dir.startswith('/eos/'):
        try:
            return list_xrdfs_root_files('root://eoscms.cern.ch', input_dir)
        except Exception as exc:
            raise RuntimeError('Could not list {} locally or via xrdfs: {}'.format(input_dir, exc))

    if input_dir.startswith('/store/'):
        try:
            return list_xrdfs_root_files('root://eoscms.cern.ch', '/eos/cms' + input_dir)
        except Exception as exc:
            raise RuntimeError('Could not list {} locally or via CERN EOS xrdfs: {}'.format(input_dir, exc))

    raise RuntimeError('Input directory does not exist locally and no remote listing rule is known: {}'.format(input_dir))


input_files = list(options.inputFiles)
if options.inputDir:
    input_files = list_root_files(options.inputDir)

if not input_files:
    input_files = ['file:TkAlDstToD0Pi.root']

process = cms.Process('TEST', eras.Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('TrackPropagation.Geant4e.geantRefit_cff')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(*input_files)
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)

process.load('Analysis.HitAnalyzer.DstToD0PiCascadeRefitProducer_cfi')
process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackKPiG4e_cfi')

process.selectedD0Tracks = process.DstToD0PiCascadeRefitProducer.clone(
    src='ALCARECOTkAlDstToD0Pi'
)

process.softPionCor = cms.EDProducer(
    'ResidualGlobalCorrectionMakerG4e',
    src=cms.InputTag('selectedD0Tracks', 'dstarTracks'),
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
    src=cms.InputTag('selectedD0Tracks', 'dstarTracks'),
    useIdealGeometry=cms.bool(False),
    fillRunTree=cms.bool(True),
    fillJac=cms.bool(True),
    respectTrackOrder=cms.bool(True),
    useDstarTriplets=cms.bool(True),
    enforceUniqueTracks=cms.bool(True),
    writeLegacyBranches=cms.bool(False),
    appendStreamIdToOutput=cms.bool(False),
    useSoftPionCvh=cms.bool(True),
    saveGenParticles=cms.bool(bool(options.isMC)),
    genParticles=cms.InputTag('genParticles'),
    softPionCvhPt=cms.InputTag('softPionCor', 'corPt'),
    softPionCvhEta=cms.InputTag('softPionCor', 'corEta'),
    softPionCvhPhi=cms.InputTag('softPionCor', 'corPhi'),
    softPionCvhCharge=cms.InputTag('softPionCor', 'corCharge'),
    softPionCvhEdmval=cms.InputTag('softPionCor', 'edmval'),
    softPionCvhMomCov=cms.InputTag('softPionCor', 'momCov'),
    minPairMass=cms.double(1.65),
    maxPairMass=cms.double(2.05),
    outprefix=cms.untracked.string('globalcor_d0_selected')
)

process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc' if options.isMC else 'auto:run2_data', '')
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
    process.selectedD0Tracks *
    process.softPionCor *
    process.globalCorD0
)
