import FWCore.ParameterSet.Config as cms

kinProcess = cms.PSet(
#    input = cms.string("/afs/cern.ch/user/p/psilva/scratch0/CMSSW_4_1_3_patch2/src/LIP/Top/data/TTJets_madgraph_Spring11.root"),
    input = cms.string("~aalves/public/DoubleElectron-v7.root"),

# data sample used
    output = cms.string("/castor/cern.ch/user/a/aalves/Dileptons/"),

    evStart = cms.int32(0),
    evEnd = cms.int32(-1),
#rum over all events???
    dirName = cms.string("evAnalyzer/data"),
    kinScheme = cms.string("std"),
    maxTries = cms.int32(10000),
# by deflaut should be 10000, more of x the algorithm solve the equations
    maxJetMult = cms.int32(2),
    mw = cms.double(80.398),
    mb = cms.double(4.8),
    ptResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt"),
    etaResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_EtaResolution_AK5PF.txt"),
    phiResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt"),
    jesUncFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt")
    )
