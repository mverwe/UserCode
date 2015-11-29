import FWCore.ParameterSet.Config as cms

PixelAnalyzer = cms.EDAnalyzer("SiPixelAnalyzer",
                               src = cms.InputTag("siPixelDigis"),
                               srcHFhits = cms.InputTag("hfreco"),
                               outputFile = cms.untracked.string("testhist.root")     
)

