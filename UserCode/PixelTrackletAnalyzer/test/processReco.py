import FWCore.ParameterSet.Config as cms

process = cms.Process("HiForest")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('isPP',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Flag if this is a pp simulation")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
#process.Timing = cms.Service("Timing") 

# MC Globaltag for 2015 dN/deta analysis
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag.globaltag = 'auto:run2_data' #'MCRUN2_74_V6B::All'
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

process.GlobalTag.toGet.extend([
 cms.PSet(record = cms.string("HeavyIonRcd"),
 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
 ## 5.02 TeV Centrality Tables
 #tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v740x01_mc"),
 #label = cms.untracked.string("HFtowersHydjetDrum5")
 ## 2.76 TeV Centrality Tables for data
 tag = cms.string("CentralityTable_HFtowers200_Glauber2010A_eff99_run1v750x01_offline"),
 label = cms.untracked.string("HFtowers")
 ),
])

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(options.inputFiles[0])
)


import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/DCSOnly/json_DCSONLY.txt').getVLuminosityBlockRange()


#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#                                                   moduleSeeds = cms.PSet(simMuonRPCDigis = cms.untracked.uint32(49835),
#                                                                          simEcalUnsuppressedDigis = cms.untracked.uint32(49835),
#                                                                          simSiStripDigis = cms.untracked.uint32(49835),
#                                                                          mix = cms.untracked.uint32(49835),
#                                                                          simHcalUnsuppressedDigis = cms.untracked.uint32(49835),
#                                                                          simMuonCSCDigis = cms.untracked.uint32(49835),
#                                                                          VtxSmeared = cms.untracked.uint32(49835),
#                                                                          g4SimHits = cms.untracked.uint32(49835),
#                                                                          simMuonDTDigis = cms.untracked.uint32(49835),
#                                                                          simSiPixelDigis = cms.untracked.uint32(49835)
#                                                                          ),
#                                                   sourceSeed = cms.untracked.uint32(49835)
#                                                   )


process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('offlinePrimaryVerticesWithBS'),
                             trackSrc = cms.untracked.InputTag('generalTracks'),
                             doTracking = cms.untracked.bool(False),
                             doMC = cms.untracked.bool(False),
                             )

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')

process.TFileService = cms.Service('TFileService',
                                   fileName=cms.string(options.outputFile))


process.analyze = cms.Path(process.hltanalysis *
                           process.siPixelRecHits*
                           process.ana)


#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

process.NoScraping = cms.EDFilter("FilterOutScraping",
 applyfilter = cms.untracked.bool(True),
 debugOn = cms.untracked.bool(False),
 numtrack = cms.untracked.uint32(10),
 thresh = cms.untracked.double(0.25)
)

process.pPAprimaryVertexFilter = cms.Path(process.PAprimaryVertexFilter)
process.pBeamScrapingFilter=cms.Path(process.NoScraping)

process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")

process.pVertexFilterCutG = cms.Path(process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.skimanalysis)
