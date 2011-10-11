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

evAnalyzer.Generator.filterSignal = cms.bool(False)
evAnalyzer.Muons.id=cms.string("TMLastStationAngTight")
evAnalyzer.Muons.usePFIso=cms.bool(True)
evAnalyzer.LooseMuons.usePFIso=cms.bool(True)
evAnalyzer.Electrons.usePFIso=cms.bool(True)
evAnalyzer.LooseElectrons.usePFIso=cms.bool(True)
