
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: reco -s RAW2DIGI --filein file:/afs/cern.ch/work/k/katatar/public/PixelReadoutDQM/0CEFE112-8E63-E511-93F7-0025905A60D0.root --conditions 75X_dataRun2_HLT_withOfflineCustomisation_v0 --no_exec --data -n 4
import FWCore.ParameterSet.Config as cms
import re

process = cms.Process('RAW2DIGI')

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('isPP',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Flag if this is a pp simulation")
options.parseArguments()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.RawToDigi_Repacked_cff')  # for /HIHighPt/HIRun2011-v1/RAW
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#
process.load('HLTrigger.Configuration.HLT_GRun_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# # Input source
# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring(
#                                                                             
# ## some files from dataset = /DoubleEG/Run2015D-v1/RAW
# #'file:/afs/cern.ch/user/k/katatar/eos/cms/store/data/Run2015D/DoubleEG/RAW/v1/000/256/673/00000/4C4110AD-E95C-E511-B6F0-02163E01396A.root'
# #'file:/afs/cern.ch/work/k/katatar/public/PixelReadoutDQM/0CEFE112-8E63-E511-93F7-0025905A60D0.root' 
#                                       ),
#     secondaryFileNames = cms.untracked.vstring()
# )

process.source = cms.Source("NewEventStreamFileReader",
                            fileNames = cms.untracked.vstring(options.inputFiles[0])   
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

outname = options.outputFile
hltOut=outname.replace('.root','HLT.root')
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(hltOut))#options.outputFile))


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('reco_RAW2DIGI.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_HLTHI_v4', '')

# Path and EndPath definitions
process.RawToDigi_custom = cms.Sequence(
    process.siPixelDigis
    +process.hcalDigis
)

process.raw2digi_step = cms.Path(process.RawToDigi_custom)

process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi")
process.hfrecoPath = cms.Path(process.hfreco)

outname = options.outputFile
siPixOut=outname.replace('.root','SiPixelAnalyzer.root')
siPixOutTree=siPixOut.replace('histos','tree')
process.siPixAna = cms.EDAnalyzer('SiPixelAnalyzer',
                              src = cms.InputTag("siPixelDigis"),
                              srcHFhits = cms.InputTag("hfreco"),
                              outputFile = cms.untracked.string(siPixOut),#options.outputFile)
                              outTreeName = cms.untracked.string(siPixOutTree)#
)

#HLT sequences
process.HLTBeginSequence = cms.Sequence( process.hltTriggerType + process.HLTL1UnpackerSequence )

process.hltSiPixelDigis.InputLabel = cms.InputTag("rawDataRepacker")
process.hltScalersRawToDigi.scalersInputTag = cms.InputTag("rawDataRepacker")
process.HLTDoLocalPixel = cms.Sequence( process.hltSiPixelDigis + process.hltSiPixelClusters)
process.HLTDoLocalHF    = cms.Sequence( process.hltHcalDigis + process.hltHfreco )

process.hltSiPixelClusters.maxNumberOfClusters = cms.int32(-1)

process.hltPixHFAct = cms.EDFilter('HLTPixelActivityHFSumEnergyFilter',
                                   inputTag = cms.InputTag("hltSiPixelClusters"),
                                   HFHitCollection = cms.InputTag("hltHfreco"),
                                   eCut_HF = cms.double(10.),
                                   eMin_HF = cms.double(10000.),
                                   offset = cms.double(-1000.),
                                   slope = cms.double(0.5)
)
process.filterPix = cms.Path(process.hltPixHFAct)

#analyzer
process.pixClusHFAna = cms.EDAnalyzer('PixerClusterHFAnalyzer',
                                      inputTag = cms.InputTag("hltSiPixelClusters"),
                                      HFHitCollection = cms.InputTag("hltHfreco"),
                                      eCut_HF = cms.double(10.),
                                      SelectEvents = cms.untracked.PSet(
                                      SelectEvents = cms.vstring('filterPix'))
)

process.HLT_PixHFAct = cms.Path( process.HLTBeginSequence + process.HLTDoLocalPixel + process.HLTDoLocalHF + process.hltPixHFAct + process.HLTEndSequence )

process.anapath = cms.Path(process.pixClusHFAna)

process.m_HLTSchedule = cms.Schedule( *(process.HLTriggerFirstPath, process.HLT_PixHFAct, process.HLTriggerFinalPath, process.HLTAnalyzerEndpath, process.anapath ))
#process.m_HLTSchedule = cms.Schedule( *(process.HLTriggerFirstPath, process.HLT_PixHFAct, process.HLTriggerFinalPath, process.HLTAnalyzerEndpath ))

process.p = cms.Path(process.siPixAna)

# Path and EndPath definitions
process.endjob_step  = cms.Path(process.endOfProcess)
#process.out_step     = cms.EndPath( process.hltTimer + process.output)
#process.out_step     = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step, process.hfrecoPath, process.p,*(process.HLTriggerFirstPath, process.HLT_PixHFAct, process.HLTriggerFinalPath, process.HLTAnalyzerEndpath, process.anapath ))
#process.schedule.extend([process.endjob_step,process.out_step])
process.schedule.extend([process.endjob_step])

#customization
oldFEDCol=cms.InputTag("rawDataCollector")
newFEDCol=cms.InputTag("rawDataRepacker")
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
for s in process.paths_().keys():
    massSearchReplaceAnyInputTag(getattr(process,s),oldFEDCol,newFEDCol)
