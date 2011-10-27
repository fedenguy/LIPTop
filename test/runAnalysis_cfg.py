import FWCore.ParameterSet.Config as cms

from LIP.Top.MVAStudyConfig_cfi import *

runProcess = cms.PSet(
    input = cms.string("/data/psilva/Top/ntuples_2011.10.27/DYJetsToLL.root"),
    outdir = cms.string("/tmp/psilva"),
    kindir = cms.string("std"),
    runSystematics = cms.bool(False),
    isMC = cms.bool(False),
    xsec = cms.double(1),
    mctruthmode = cms.int32(0),
    saveSummaryTree = cms.bool(False),
    evStart = cms.int32(0),
    evEnd = cms.int32(-1),
    dirName = cms.string("evAnalyzer/data"),
    ptResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt"),
    etaResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_EtaResolution_AK5PF.txt"),
    phiResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt"),
    jesUncFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt"),
    useMVA = cms.bool(False),
    tmvaInput =	chargedHiggsStudy
    )
