import FWCore.ParameterSet.Config as cms

# simple mass study
pairSelStudy = cms.PSet(  doTrain = cms.bool(False), 
                          studyTag   = cms.string("pairSelStudy"),
                          weightsDir = cms.string("${CMSSW_BASE}/src/LIP/Top/weights"),
                          methodList = cms.vstring('MLP','BDT','Likelihood'),
                          #varsList   = cms.vstring("k10p50p","k25p50p","k75p50p","k90p50p","k90p10p","k75p25p","k75p10p","k90p25p","k90p75p")
                          varsList   = cms.vstring('kMPV','kRMS','k10p50p','k25p50p','kIntegral')
                          #varsList   = cms.vstring("k90p50p","k10p50p","k25p50p","kIntegral")
                          )


# simple study
chargedHiggsStudy = cms.PSet( evCategories = cms.vint32(),
                              input      = cms.string("/data/psilva/Top/chargedHiggs/EventSummaries.root"),
                              studyTag   = cms.string("simpleChHiggsDiscriminator"),
                              weightsDir = cms.string("${CMSSW_BASE}/src/LIP/Top/weights"),
                              methodList = cms.vstring('Likelihood'),#'BDT','MLP'),
                              varsList   = cms.vstring("mindrlj","mljmostiso","mtleastisolep","mtmostisolep"),#"drmostisol","drleastisol"),
                              #                               varsList   = cms.vstring('mljmostiso','ptleastisollep','ptmostisollep','dphimetleastisolep','mtleastisolep'),
                              procList   = cms.vstring('WH_120_Summer11',
                                                       #                             'SingleT_tW',
                                                       #                             'SingleTbar_tW',
                                                       #                             'TTJets',
                                                       'TTJets_signal',
                                                       #                             'WW','WZ','ZZ',
                                                       #                             'DYJetsToLL'
                                                       ),
                              procType   = cms.vint32 (1,2),#2,2,2,2,2,2,2),
                              procWeight = cms.vdouble(1,1),#1,1,1,1,1,1,1)
                               )




