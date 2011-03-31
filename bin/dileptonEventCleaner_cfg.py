import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('rfio:///castor/cern.ch/user/p/psilva/Dileptons/DYToEEM20/evHyp_0_5.root'),
    maxEvents   = cms.int32(-1),                       
    outputEvery = cms.uint32(10)                       
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('Histograms.root')
)

from CMGTools.HtoZZ2l2nu.StandardSelections_cfi import *
process.evAnalyzer = cms.PSet(
    Generator = BaseGeneratorSelection.clone(),
    Vertices = BaseVertexSelection.clone(),
    Muons = BaseMuonsSelection.clone(),
    Electrons = BaseElectronsSelection.clone(),
    Dileptons = BaseDileptonSelection.clone(),
    Jets = BaseJetSelection.clone(),
    MET = BaseMetSelection.clone()
    )
