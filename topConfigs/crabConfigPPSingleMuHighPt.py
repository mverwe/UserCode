from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'FilteredSingleMuHighPt_v3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False


config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runForestAOD_pp_DATA_75X_SingleMuHighPt.py'

config.section_('Data')
config.Data.inputDataset = '/SingleMuHighPt/Run2015E-PromptReco-v1/AOD'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 50000
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262328_5TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.outLFNDirBase = '/store/group/cmst3/group/hintt/mverweij/PP5TeV/data'
config.Data.publication = False #True
config.Data.outputDatasetTag = ''

config.section_('User')
config.section_('Site')
#config.Site.whitelist = ['T2_US_MIT']
#config.Site.blacklist = ['T2_US_Nebraska','T2_US_Florida','T2_US_Wisconsin','T2_US_Caltech']
config.Site.storageSite = 'T2_CH_CERN'

