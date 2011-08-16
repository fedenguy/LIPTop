import FWCore.ParameterSet.Config as cms

# simple study
pairSelStudy = cms.PSet(  evCategories = cms.vint32(0,1,2,3),
                          input      = cms.string("EventSummaries.root"),
                          studyTag   = cms.string("pairSelStudy"),
                          weightsDir = cms.string("${CMSSW_BASE}/src/LIP/Top/weights"),
                          methodList = cms.vstring('BDT'),
                          varsList   = cms.vstring('dphill', 'mtsum', 'redMetoverzpt', 'ptsum', 'maxdrlz' ),
                          #                         procList   = cms.vstring( 'VBFtoH200toZZto2L2Nu', 'GGtoH200toZZto2L2Nu', 'GGtoH200toWWto2L2Nu', 'ZZto2l2Nu', 'WWTo2L2Nu', 'WZTo3LNu'),
                          #                         procType   = cms.vint32 ( 1                     , 1                    ,  1                   , 2          , 2          ,  2),
                          #                         procWeight = cms.vdouble( 1                     , 1                    ,  1                   , 1          , 1          ,  1),
                          #                         procList   = cms.vstring( 'VBFtoH200toZZto2L2Nu', 'GGtoH200toZZto2L2Nu',  'ZZto2l2Nu', 'WWTo2L2Nu', 'WZTo3LNu'),
                          #                         procType   = cms.vint32 ( 1                     , 1                    ,  2          , 2          ,  2),
                          #                         procWeight = cms.vdouble( 1                     , 1                    ,  1          , 1          ,  1),
                          procList   = cms.vstring( 'VBFtoH200toZZto2L2Nu',  'GGtoH200toWWto2L2Nu', 'ZZto2l2Nu', 'WWTo2L2Nu', 'WZTo3LNu'),
                          procType   = cms.vint32 ( 1                     ,   1                   , 2          , 2          ,  2),
                          procWeight = cms.vdouble( 1                     ,   1                   , 1          , 1          ,  1),
                          )
