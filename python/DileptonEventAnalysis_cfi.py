import FWCore.ParameterSet.Config as cms

from CMGTools.HtoZZ2l2nu.StandardSelections_cfi import *
evAnalyzer = cms.EDAnalyzer("DileptonEventCleaner",
                            Trigger = BaseTriggerSelection.clone(),
                            Generator = BaseGeneratorSelection.clone(),
                            Vertices = BaseVertexSelection.clone(),
                            Muons = BaseMuonsSelection.clone(),
                            LooseMuons = BaseLooseMuonsSelection.clone(),
                            Electrons = BaseElectronsSelection.clone(),
                            LooseElectrons = BaseLooseElectronsSelection.clone(),
                            Dileptons = BaseDileptonSelection.clone(),
                            Jets = BaseJetSelection.clone(),
                            MET = BaseMetSelection.clone()
                            )
