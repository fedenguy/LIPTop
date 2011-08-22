import FWCore.ParameterSet.Config as cms

# simple study
pairSelStudy = cms.PSet(  doTrain = cms.bool(True), 
                          studyTag   = cms.string("pairSelStudy"),
                          weightsDir = cms.string("${CMSSW_BASE}/src/LIP/Top/weights"),
                          methodList = cms.vstring('Likelihood'),
                          varsList   = cms.vstring('kIntegral', 'kMPV', 'kMean', 'kRMS', 'kSkewness', 'kKurtosis', 'k10p', 'k25p', 'k75p', 'k90p')
                          )
