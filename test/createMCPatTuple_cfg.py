import os
runOnMC=True
trigFilter=''
cfgFile=os.path.expandvars('${CMSSW_BASE}/src/LIP/Top/test/createPatTuple_cfg.py')
execfile(cfgFile)
