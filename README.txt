#
# PAT TUPLE/ NTUPLE PRODUCTION: check README.txt under the CMGTools/HtoZZ2l2nu 
#

#
# CREATE BASE DISTRIBUTION AND EVENT SUMMARIES
#
runLocalAnalysisOverSamples.py -e showControlDistributions -j $CMSSW_BASE/src/LIP/Top/data/samples_2012.json -d /store/cmst3/user/psilva/28May2012_CMSSW444_HZZ2l2v_ntuples -o $CMSSW_BASE/src/LIP/Top/test/results -c test/runAnalysis_cfg.py.templ -p "@sfMetCut=40 @ofMetCut=0 @jetPtCut=30 @runSystematics=True" -s 8nh
runPlotter --iLumi 5041 --inDir test/results/ --outDir test/results/plots --json data/samples_2012.json


#
# RUN KIN RECONSTRUCTION
#
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=std" -s 2nd
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=jesup" -s 2nd -t TT
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=jesdown" -s 2nd -t TT
runKinOverSamples.py -j data/samples.json -e 100 -d /store/cmst3/user/psilva/Top_ntuples_14Jan -p "-out=/castor/cern.ch/user/p/psilva/Top -run=jer" -s 2nd -t TT

runKinOverSamples.py -j data/mutau-samples.json -e 100 -d /store/cmst3/user/psilva/Top_mutau_ntuples -p "-out=/castor/cern.ch/user/p/psilva/Top -run=std" -s 2nd
runKinOverSamples.py -j data/tau-mass-samples.json -e 100 -d /store/cmst3/user/psilva/Top_mutau_ntuples -p "-out=/castor/cern.ch/user/p/psilva/Top -run=std" -s 2nd

runKinOverSamples.py -j data/etau-samples.json -e 10 -d /store/cmst3/user/psilva/Top_etau_ntuples -p "-out=/castor/cern.ch/user/p/psilva/Top_etau -run=std" -s 8nh
runKinOverSamples.py -j data/tau-mass-samples.json -e 10 -d /store/cmst3/user/psilva/Top_etau_ntuples -p "-out=/castor/cern.ch/user/p/psilva/Top_etau -run=std" -s 8nh


# generate plots for signal templates
runLocalAnalysisOverSamples.py -e showMassDistribution -j data/samples-signal.json -d store/Top/ntuples_14Jan/merged/ -o /tmp/psilva -c test/runAnalysis_cfg.py.templ -l 4565 -p "@sfMetCut=40 @ofMetCut=30 @jetPtCut=35 @kindir=std"
runLocalAnalysisOverSamples.py -e showMassDistribution -j data/mass-samples.json -d store/Top/ntuples_14Jan/merged/ -o /tmp/psilva -c test/runAnalysis_cfg.py.templ -l 4565 -p "@sfMetCut=40 @ofMetCut=30 @jetPtCut=35 @kindir=std"
runLocalAnalysisOverSamples.py -e showMassDistribution -j data/syst-samples.json -d store/Top/ntuples_14Jan/merged/ -o /tmp/psilva -c test/runAnalysis_cfg.py.templ -l 4565 -p "@sfMetCut=40 @ofMetCut=30 @jetPtCut=35 @kindir=std"

# to run PDF weights...
#for i in `seq 0 99`; do \
#for i in `seq 100 199`; do \
#for i in `seq 250 299`; do \
for i in `seq 300 499`; do \
let a=50*$i;\
let b=50*$i+50;\
runLocalAnalysisOverSamples.py -e showMassDistribution -j ../data/samples-signal.json -d ../store/Top/ntuples_14Jan/merged -o /tmp/psilva -c ../test/runAnalysis_cfg.py.templ -l 2165 -p "@kindir=std @runSystematics=True @saveSummaryTree=True @evStart=$a @evEnd=$b" -t signal &\
sleep 50
done

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
