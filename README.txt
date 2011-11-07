#
# submit ntuple production
#
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -s True -d patdir -t QCD
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -s True -d patdir -t WJets
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t ZZ
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 True -d patdir -t WZ
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 True -d patdir -t WW
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t TT
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t DY
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t SingleT
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t May10
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t v4
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t Aug05
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t v6
