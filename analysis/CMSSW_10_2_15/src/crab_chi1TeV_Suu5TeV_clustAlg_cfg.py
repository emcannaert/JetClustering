from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'Suu_clusteringAlgorithm_0001' #'diquark_QCD_background_data_0002' #
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'clusteringAnalyzer_cfg.py'
#config.Data.inputDataset = '/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
#config.Data.inputDataset = '/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
#config.Data.inputDataset = '/JetHT/Run2018C-17Sep2018-v1/MINIAOD'
config.JobType.inputFiles=['file:Thrust.h','file:analyzers/analyzers/src/Thrust.cc']

config.Data.inputDataset = '/SuuToWbht/ecannaer-CRAB_diquark_MC_prod_chi1TeV_Suu5TeV_MINIAOD-0bd58594e6ade05f64e0c3a8301c3139/USER'
config.Data.inputDBS = 'phys03'
#config.Data.inputDataset = '/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
#config.Data.inputDataset = '/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.publication = False
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.lumiMask = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.outputDatasetTag = 'diquark_5TeV_chi1TeV_clustAlg'
config.Site.storageSite = 'T3_US_FNALLPC'
