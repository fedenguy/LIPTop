import FWCore.ParameterSet.Config as cms

kinProcess = cms.PSet(
#    input = cms.string("/afs/cern.ch/user/p/psilva/scratch0/CMSSW_4_1_3_patch2/src/LIP/Top/data/TTJets_madgraph_Spring11.root"),
    input = cms.string("~aalves/public/DoubleMuon-v3.root"),
    output = cms.string("kin.root"),
    evStart = cms.int32(0),
    evEnd = cms.int32(-1),
    dirName = cms.string("evAnalyzer/data"),
    kinScheme = cms.string("std"),
    maxTries = cms.int32(1000),
    maxJetMult = cms.int32(2),
    mw = cms.double(80.398),
    mb = cms.double(4.8),
    ptResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt"),
    etaResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_EtaResolution_AK5PF.txt"),
    phiResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt"),
    jesUncFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt")
    )
