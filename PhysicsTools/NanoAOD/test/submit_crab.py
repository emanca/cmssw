import os

samples={}

samples['w']=[

	'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
	'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2/MINIAODSIM',
	'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v5/MINIAODSIM'

]

for resonance,slist in samples.iteritems():
    for i,s in enumerate(slist):
        filename="crabSubmit_"+resonance+"_"+str(i)+'.py'
        f=open(filename,'w')
        cfg="""

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = '{name}'
config.General.workArea = 'CRAB'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'myNanoProdMc_NANO.py'

config.Data.inputDataset = '{dataset}'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/NanoWMassV4/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'NanoWMass'

config.Site.storageSite = 'T2_IT_Pisa'

""".format(name=resonance+"_"+str(i),dataset=s)

        f.write(cfg)
        f.close()



        os.system('crab submit --config='+filename)
