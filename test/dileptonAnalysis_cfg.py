import FWCore.ParameterSet.Config as cms

process = cms.Process("DileptonAN")

#input
from CMGTools.HtoZZ2l2nu.localPatTuples_cff import configureFromCommandLine
castorDir, outFile, inputList = configureFromCommandLine()
process.source = cms.Source("PoolSource",
                            fileNames = inputList
                            #cms.untracked.vstring('file:/afs/cern.ch/user/f/federica/public/CMSSW_4_2_0/src/FastSimulation/ParticleFlow/test/PAT_FastSim_modified_cfg.root')
                            )

print inputList    
#load the analyzer
process.load('LIP.Top.TopDileptonEventProducer_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventFilter_cfi')
process.load('CMGTools.HtoZZ2l2nu.PileupNormalizationProducer_cfi')
process.load('LIP.Top.DileptonEventAnalysis_cfi')
process.evAnalyzer.dtag=cms.string('top')
process.TFileService = cms.Service("TFileService", fileName = cms.string(outFile) )

#process the analysis
#process.p = cms.Path(process.puWeights*process.cleanEvent*process.cleanEventFilter*process.evAnalyzer)
process.p = cms.Path(process.cleanEvent*process.cleanEventFilter*process.evAnalyzer)

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

