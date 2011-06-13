from PhysicsTools.PatAlgos.patTemplate_cfg import *

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if ( not runOnMC ): process.GlobalTag.globaltag = 'GR_R_42_V2::All'
else:               process.GlobalTag.globaltag = 'START42_V12::All'

# global options
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                      SkipEvent = cms.untracked.vstring('ProductNotFound')
                                      )
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# event source  
from CMGTools.HtoZZ2l2nu.localPatTuples_cff import fillFromCastor
process.source.fileNames=fillFromCastor('/castor/cern.ch/cms/store/cmst3/user/cbern/CMG/TT_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# event output
from CMGTools.HtoZZ2l2nu.OutputConfiguration_cff import configureOutput
configureOutput(process)

# preselection filters
from CMGTools.HtoZZ2l2nu.PreselectionSequences_cff import addPreselectionSequences,addLumifilter
addPreselectionSequences(process)
if(not runOnMC ): addLumifilter(process, '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-166502_7TeV_PromptReco_Collisions11_JSON.txt')

# trigger filter
from CMGTools.HtoZZ2l2nu.TriggerSequences_cff import addTriggerSequence
if(not runOnMC): addTriggerSequence(process,trigFilter)

# pat sequences
from CMGTools.HtoZZ2l2nu.PatSequences_cff import addPatSequence
addPatSequence(process,runOnMC)

# dilepton filters
from CMGTools.HtoZZ2l2nu.DileptonFilterSequences_cff import addDileptonFilters
addDileptonFilters(process)

# event counters
process.startCounter = cms.EDProducer("EventCountProducer")
process.endCounter = process.startCounter.clone()

# define the paths
if(runOnMC):
    process.eePath = cms.Path(process.startCounter * process.patSequence * process.eeCandidateSequence )
    process.mumuPath  = cms.Path(process.startCounter * process.patSequence * process.mumuCandidateSequence )
    process.emuPath  = cms.Path(process.startCounter * process.patSequence * process.emuCandidateSequence )
else:
    process.eePath = cms.Path(process.startCounter * process.preselection * process.trigSequence * process.patSequence * process.eeCandidateSequence )
    process.mumuPath  = cms.Path(process.startCounter * process.preselection * process.trigSequence * process.patSequence * process.mumuCandidateSequence )
    process.emuPath  = cms.Path(process.startCounter * process.preselection * process.trigSequence * process.patSequence * process.emuCandidateSequence )
process.e = cms.EndPath( process.endCounter*process.out )

# all done, schedule the execution
if(runOnMC):
    from CMGTools.HtoZZ2l2nu.GeneratorLevelSequences_cff import addGeneratorLevelSequence
    addGeneratorLevelSequence(process)
    process.schedule = cms.Schedule( process.genLevelPath, process.eePath, process.mumuPath, process.emuPath, process.e )
else :
    process.schedule = cms.Schedule( process.eePath, process.mumuPath, process.emuPath, process.e )

print "Scheduling the following modules"
print process.schedule
print "Events will be selected for the following paths:"
print process.out.SelectEvents.SelectEvents

