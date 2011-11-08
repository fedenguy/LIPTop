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

#
# create control distributions and event summaries
#
# submit the local analysis to batch
runLocalAnalysisOverSamples.py -e showControlDistributions -o ${HOME}/scratch0/top -d /castor/cern.ch/cms/store/cmst3/user/psilva/Top/ntuples_2011.10.31 -j data/samples-signal.json -p "@saveSummaryTree=True @runSystematics=True"  -c test/runAnalysis_cfg.py.templ -s True

# plot the results
runPlotter --iLumi 1947 --inDir ~/scratch0/top/ --outDir /tmp/psilva/ --json data/samples-signal.json

# create a single event summary file
hadd EventSummaries.root ${HOME}/scratch0/*_summary.root 
rm ${HOME}/scratch0/*_summary.root


