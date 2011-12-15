import FWCore.ParameterSet.Config as cms

# base values for trigger event
from CMGTools.HtoZZ2l2nu.TriggerSequences_cff import getTriggerPaths
doubleEle, doubleMu, muEG, singleEle, singleMu, mcTrigs = getTriggerPaths()
BaseTriggerSelection = cms.PSet( source = cms.InputTag("TriggerResults::HLT"),
                                 mumuTriggerPaths = cms.vstring(doubleMu),
                                 emuTriggerPaths  = cms.vstring(muEG),
                                 eeTriggerPaths   = cms.vstring(doubleEle)
                                 )

# base values for the vertex selection ------------------------------------------
BaseGeneratorSelection = cms.PSet( source = cms.InputTag("prunedGen"),
                                   filterId = cms.int32(6),
                                   genJets=cms.InputTag("ak5GenJets"),
                                   puReweight=cms.InputTag("puWeights:puWeight"),
                                   normPuReweight=cms.InputTag("puWeights:renPuWeight"),
                                   filterDYmassWindow = cms.bool(False),
                                   mcpileup = cms.string('/afs/cern.ch/user/p/psilva/public/Pileup/Summer11MC.root'),
                                   datapileup = cms.string('/afs/cern.ch/user/p/psilva/public/Pileup/Pileup2011A.root')
                                   )


# base values for the vertex selection ------------------------------------------
BaseVertexSelection = cms.PSet( source = cms.InputTag("offlinePrimaryVertices"),
                                beamSpot = cms.InputTag("offlineBeamSpot"),
                                maxZ = cms.double(24),
                                maxRho = cms.double(2.0),
                                minNDOF = cms.int32(4)
                                )

# base values for muon selection ----------------------------------------------
BaseMuonsSelection = cms.PSet( source = cms.InputTag("selectedPatMuonsPFlow"),
                               minPt = cms.double(20),
                               maxEta = cms.double(2.4),
                               requireGlobal = cms.bool(True),
                               requireTracker = cms.bool(True),
                               maxTrackChi2 = cms.double(10),
                               minValidTrackerHits = cms.int32(11),
                               minValidMuonHits=cms.int32(1),
                               maxDistToBeamSpot=cms.double(0.02),
                               id = cms.string(""),
                               maxRelIso = cms.double(0.2),
                               usePFIso = cms.bool(True)
                               )

# base values for loose muon selection ----------------------------------------------
BaseLooseMuonsSelection = BaseMuonsSelection.clone( minPt = cms.double(5),
                                                    maxTrackChi2 = cms.double(9999.),
                                                    id = cms.string(""),
                                                    maxRelIso = cms.double(0.5),
                                                    minValidMuonHits=cms.int32(0)
                                                    )

# base values for electron selection ----------------------------------------------
BaseElectronsSelection = cms.PSet( source = cms.InputTag("selectedPatElectronsPFlow"),
                                   minPt = cms.double(20),
                                   minSuperClusterEt = cms.double(5),
                                   maxEta = cms.double(2.5),
                                   vetoTransitionElectrons = cms.bool(False),
                                   applyConversionVeto = cms.bool(True),
                                   maxDistToBeamSpot=cms.double(0.04),
                                   maxTrackLostHits = cms.int32(1),
                                   id = cms.string("eidVBTF90"),
                                   maxRelIso = cms.double(0.17),
                                   minDeltaRtoMuons = cms.double(0.1),
                                   usePFIso = cms.bool(True)
                                   )

# base values for electron selection ----------------------------------------------
BaseLooseElectronsSelection = BaseElectronsSelection.clone( minPt = cms.double(5),
                                                            minSuperClusterEt = cms.double(5),
                                                            vetoTransitionElectrons = cms.bool(False),
                                                            maxTrackLostHits = cms.int32(2),
                                                            id = cms.string("eidVBTF95"),
                                                            maxRelIso = cms.double(0.5)
                                                            )

#my base values for jet selection -------------------------------------------------
BaseJetSelection = cms.PSet( source = cms.InputTag("selectedPatJetsPFlow"),
                             rho = cms.InputTag("kt6PFJetsPFlow:rho"),
                             jetId = cms.PSet( version = cms.string("FIRSTDATA"), quality = cms.string("LOOSE") ),
                             minPt = cms.double(15.),
                             maxEta = cms.double(5.0),
                             minDeltaRtoLepton = cms.double(-1)
                             )

# base values for the dilepton selection ------------------------------------------
BaseDileptonSelection = cms.PSet( minDileptonMass = cms.double(12),
                                  maxDileptonMass = cms.double(9999999.)
                                  )

# base values for met selection -----------------------------------------------------
BaseMetSelection = cms.PSet( source = cms.InputTag("patMETsPFlow") ) 

