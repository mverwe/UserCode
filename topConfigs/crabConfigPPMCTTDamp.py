from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'MCTTDamp_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False


config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runForestAOD_pp_MC_75X.py'

config.section_('Data')
config.Data.inputDataset = '/TTbar_hdamp_NNPDF30_TuneCUETP8M1_5020GeV-powheg/HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1/AODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.splitting = "EventAwareLumiBased"
#config.Data.unitsPerJob = 40000

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
config.Data.outLFNDirBase = '/store/group/cmst3/group/hintt/mverweij/PP5TeV/MC'
config.Data.publication = False #True
config.Data.outputDatasetTag = ''

config.section_('User')
config.section_('Site')
#config.Site.whitelist = ['T2_US_MIT']
#config.Site.blacklist = ['T2_US_Nebraska','T2_US_Florida','T2_US_Wisconsin','T2_US_Caltech']
config.Site.storageSite = 'T2_CH_CERN'

