# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step4 --python_filename DStarToD0Pi_D0KPi_TuneCP5_RECOALCA_cfg.py --eventcontent RECOSIM,ALCARECO --datatier GEN-SIM-RECO,ALCARECO --filein file:DStarToD0Pi_D0KPi_TuneCP5_HLT.root --fileout file:DStarToD0Pi_D0KPi_TuneCP5_RECOALCA.root --conditions 106X_mcRun2_asymptotic_v13 --customise_commands process.ALCARECOStreamTkAlDstToD0Pi.outputCommands.extend(process.GeneratorInterfaceAOD.outputCommands)\nprocess.ALCARECOStreamTkAlDstToD0Pi.outputCommands.extend(process.SimGeneralAOD.outputCommands) --step RAW2DIGI,L1Reco,RECO,RECOSIM,ALCA:TkAlDstToD0Pi --geometry DB:Extended --era Run2_2016 --runUnscheduled --mc -n -1 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

process = cms.Process('RECO',Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreamsMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:DStarToD0Pi_D0KPi_TuneCP5_HLT.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step4 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:DStarToD0Pi_D0KPi_TuneCP5_RECOALCA.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.ALCARECOStreamTkAlDstToD0Pi = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlDstToD0Pi')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('TkAlDstToD0Pi')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('TkAlDstToD0Pi.root'),
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'keep *_ALCARECOTkAlDstToD0Pi_*_*', 
        'keep L1AcceptBunchCrossings_*_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_TriggerResults_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep *_offlinePrimaryVertices_*_*'
    )
)

# Other statements
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOTkAlDstToD0Pi_noDrop.outputCommands)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v13', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.ALCARECOStreamTkAlDstToD0PiOutPath = cms.EndPath(process.ALCARECOStreamTkAlDstToD0Pi)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.pathALCARECOTkAlDstToD0Pi,process.endjob_step,process.RECOSIMoutput_step,process.ALCARECOStreamTkAlDstToD0PiOutPath)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

process.ALCARECOStreamTkAlDstToD0Pi.outputCommands.extend(process.GeneratorInterfaceAOD.outputCommands)
process.ALCARECOStreamTkAlDstToD0Pi.outputCommands.extend(process.SimGeneralAOD.outputCommands)
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
