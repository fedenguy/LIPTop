import FWCore.ParameterSet.Config as cms

# simple study
pairSelStudy = cms.PSet(  doTrain = cms.bool(False), 
                          studyTag   = cms.string("pairSelStudy"),
                          weightsDir = cms.string("${CMSSW_BASE}/src/LIP/Top/weights"),
                          methodList = cms.vstring('MLP','BDT','Likelihood'),
                          #varsList   = cms.vstring("k10p50p","k25p50p","k75p50p","k90p50p","k90p10p","k75p25p","k75p10p","k90p25p","k90p75p")
                          varsList   = cms.vstring('kMPV','kRMS','k10p50p','k25p50p','kIntegral')
                          #varsList   = cms.vstring("k90p50p","k10p50p","k25p50p","kIntegral")
                          )
