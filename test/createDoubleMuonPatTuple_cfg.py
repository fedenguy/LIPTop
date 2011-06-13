import os,sys
runOnMC=False
trigFilter='mumu'
cfgFile=os.path.expandvars('${CMSSW_BASE}/src/LIP/Top/test/createPatTuple_cfg.py')
castorDir='/castor/cern.ch/cms/store/cmst3/user/cbern/CMG/TT_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'
if(len(sys.argv)>2 ): castorDir=sys.argv[2]
execfile(cfgFile)
