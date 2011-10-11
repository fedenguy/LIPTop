import FWCore.ParameterSet.Config as cms

from CMGTools.HtoZZ2l2nu.StandardSelections_cfi import *
cleanEvent = cms.EDProducer("TopDileptonEventProducer",
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
cleanEvent.Generator.filterSignal = cms.bool(False)
cleanEvent.Muons.id=cms.string("TMLastStationAngTight")
cleanEvent.Muons.usePFIso=cms.bool(True)
cleanEvent.LooseMuons.usePFIso=cms.bool(True)
cleanEvent.Electrons.usePFIso=cms.bool(True)
cleanEvent.LooseElectrons.usePFIso=cms.bool(True)
