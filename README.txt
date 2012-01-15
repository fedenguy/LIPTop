#
# PAT Tuple production
#

crab -report -c dir
lumiCalc2.py -c frontier://LumiCalc/CMS_LUMI_PROD -i dir/res/summary.json overview -b stable -norm pp7TeV > dir.lumi

#
# NTUPLE PRODUCTION
# use the lxbatch to run over the pat-tuples and create the ntuples
#
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 100 -s 8nh -d patdir -t QCD
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 50 -s 8nh -d patdir -t WJets
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t ZZ
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t WZ
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t WW
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t TT
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 5 -s 8nh -d patdir -t DY
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t SingleT
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t May10
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t v4
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t Aug05
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t v6
runOverSamples.py -j data/samples.json -p "-cfg=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/test/dileptonAnalysis_cfg.py -castor=${HOME}/scratch0/CMSSW_4_2_4/src/LIP/Top/data" -n 1 -s 8nh -d patdir -t 2011B

# merge the outputs (hadd)
./test/mergeOutputs.sh

#
# CREATE CONTROL DISTRIBUTIONS AND EVENT SUMMARIES
# it will run an executable and store the outputs in root files named after each sample
#
runLocalAnalysisOverSamples.py -e showControlDistributions -o ${HOME}/scratch0/top -d /store/cmst3/user/psilva/Top_ntuples -j data/samples-signal.json -p "@saveSummaryTree=True @runSystematics=True"  -c test/runAnalysis_cfg.py.templ -s 2nd
runLocalAnalysisOverSamples.py -e showControlDistributions -o ${HOME}/scratch0/top-nosyst -d /store/cmst3/user/psilva/Top_ntuples -j data/samples-signal.json  -c test/runAnalysis_cfg.py.templ  -p "@saveSummaryTree=True" -s 8nh
runLocalAnalysisOverSamples.py -e showControlDistributions -o ${HOME}/scratch0/top-newjec -d /store/cmst3/user/psilva/Top_ntuples_14Jan -j data/samples-signal.json  -c test/runAnalysis_cfg.py.templ -s 8nh


#
# PLOTTING THE RESULTS
# besides creating the plots it will collect all in a file called plotter.root
#
runPlotter --iLumi 2165 --inDir ~/scratch0/top/ --outDir /tmp/psilva/ --json data/samples-signal.json

#
# RUN KIN RECONSTRUCTION
#
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=std" -s 2nd
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=jesup" -s 2nd -t TT
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=jesdown" -s 2nd -t TT
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=jer" -s 2nd -t TT

runKinOverSamples.py -j data/mutau-samples.json -e 100 -d /store/cmst3/user/psilva/Top_mutau_ntuples -p "-out=/castor/cern.ch/user/p/psilva/Top -run=std" -s 2nd
runKinOverSamples.py -j data/mutau-mass-samples.json -e 100 -d /store/cmst3/user/psilva/Top_mutau_ntuples -p "-out=/castor/cern.ch/user/p/psilva/Top -run=std" -s 2nd


# generate plots for signal templates
runLocalAnalysisOverSamples.py -e showMassDistribution -j data/samples-signal.json -d store/Top/kin -o /tmp/psilva -c test/runAnalysis_cfg.py.templ -l 2165 -p "@kindir=std"
runLocalAnalysisOverSamples.py -e showMassDistribution -j data/mass-samples.json -d store/Top/kin -o /tmp/psilva -c test/runAnalysis_cfg.py.templ -l 2165 -p "@kindir=std"
runLocalAnalysisOverSamples.py -e showMassDistribution -j data/syst-samples.json -d store/Top/kin -o /tmp/psilva -c test/runAnalysis_cfg.py.templ -l 2165 -p "@kindir=std"

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
python runHFCMeasurement.py -l 2165 -j ../../data/samples-signal.json -i ~/scratch0/top/EventSummary.root -f 1 -a TCHEL
python runHFCMeasurement.py -l 2165 -j ../../data/samples-signal.json -i ~/scratch0/top/EventSummary.root -f 5 -a TCHEL

# fit R


#
# OTHERS
#
# mount user area in pclip11/maccms316
mount_store.sh

#when you're done unmount it
fusermount -u store
rm -rf store
