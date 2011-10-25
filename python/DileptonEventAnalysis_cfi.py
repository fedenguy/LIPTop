import FWCore.ParameterSet.Config as cms

from LIP.Top.StandardSelections_cfi import *
evAnalyzer = cms.EDAnalyzer("TopDileptonEventAnalyzer",
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
evAnalyzer.Generator.filterSignal = cms.bool(False)
