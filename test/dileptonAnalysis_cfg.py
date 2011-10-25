import FWCore.ParameterSet.Config as cms

process = cms.Process("DileptonAN")

#input
from CMGTools.HtoZZ2l2nu.localPatTuples_cff import configureFromCommandLine
castorDir, outFile, inputList = configureFromCommandLine()
process.source = cms.Source("PoolSource",
                            fileNames = inputList
                            )

print inputList    
#load the analyzer
process.load('CMGTools.HtoZZ2l2nu.PileupNormalizationProducer_cfi')
process.load('LIP.Top.DileptonEventAnalysis_cfi')
process.TFileService = cms.Service("TFileService", fileName = cms.string(outFile) )

if(outFile.find('DYJets')>=0 and outFile.find('M20to50')>=0 ) :
    process.evAnalyzer.Generator.filterDYmassWindow=cms.bool(True) 

from CMGTools.HtoZZ2l2nu.PreselectionSequences_cff import addLumifilter
if(outFile.find('May10ReReco')>=0):
    addLumifilter(process,'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt')
if(outFile.find('PromptReco')>=0):
    addLumifilter(process,'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-178078_7TeV_PromptReco_Collisions11_JSON.txt')
if(outFile.find('05AugReReco')>=0):
    addLumifilter(process,'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt')


#process the analysis
process.p = cms.Path(process.puWeights*process.evAnalyzer)

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

