import os
runOnMC=False
trigFilter='mumu'
cfgFile=os.path.expandvars('${CMSSW_BASE}/src/LIP/Top/test/createPatTuple_cfg.py')
execfile(cfgFile)
