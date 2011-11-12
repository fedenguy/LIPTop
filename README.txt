#
# PAT Tuple production
#

crab -report -c dir
lumiCalc2.py -c frontier://LumiCalc/CMS_LUMI_PROD -i dir/res/summary.json overview -b stable -norm pp7TeV > dir.lumi

#
# NTUPLE PRODUCTION
# use the lxbatch to run over the pat-tuples and create the ntuples
#
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -s True -d patdir -t QCD
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5-s True -d patdir -t WJets
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t ZZ
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t WZ
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t WW
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t TT
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t DY
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t SingleT
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t May10
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t v4
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t Aug05
runOverSamples.py -j data/samples-signal.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s True -d patdir -t v6

# merge the outputs (hadd)
./test/mergeOutputs.sh

#
# CREATE CONTROL DISTRIBUTIONS AND EVENT SUMMARIES
# it will run an executable and store the outputs in root files named after each sample
#
runLocalAnalysisOverSamples.py -e showControlDistributions -o ${HOME}/scratch0/top -d /castor/cern.ch/cms/store/cmst3/user/psilva/Top/ntuples_2011.11.11 -j data/samples-signal.json -p "@saveSummaryTree=True @runSystematics=True"  -c test/runAnalysis_cfg.py.templ -s True


#
# PLOTTING THE RESULTS
# besides creating the plots it will collect all in a file called plotter.root
#
runPlotter --iLumi 1947 --inDir ~/scratch0/top/ --outDir /tmp/psilva/ --json data/samples-signal.json


#
# DY control
#
# count events under the Z mass peak ee/mumu
python bin/macros/getDileptonCutflow.py -c dilmassctr -b 2 -i plotter.root 

# fit the sum MT distribution for the emu channel
root -l bin/macros/fitDYdistribution.C

#
# EVENT YIELDS
# 
python bin/macros/getDileptonCutflow.py -c evtflow -b 4 -i plotter.root 



# create a single event summary file
hadd EventSummaries.root ${HOME}/scratch0/*_summary.root 
rm ${HOME}/scratch0/*_summary.root


#
# HFC studies
# the scripts are in bin/HFC
# the fit parameters are configured in hfcFitter_cfg.json
# fit b-tag efficiency
python runHFCMeasurement.py -l 1957 -j ../../data/samples-signal.json -i ~/scratch0/top/EventSummary.root -f 1 -a TCHEL
python runHFCMeasurement.py -l 1957 -j ../../data/samples-signal.json -i ~/scratch0/top/EventSummary.root -f 5 -a TCHEL

# fit R


#
# OTHERS
#
# mount user area in pclip11
./test/mount_store.sh

#when you're done unmount it
fusermount -u store
rm -rf store
