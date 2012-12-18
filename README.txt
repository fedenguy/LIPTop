#
# PAT TUPLE/ NTUPLE PRODUCTION: check README.txt under the CMGTools/HtoZZ2l2nu 
#

#
# CREATE BASE DISTRIBUTIONS AND ESTIMATE BACKGROUNDS
#
runLocalAnalysisOverSamples.py -e showControlDistributions -j $CMSSW_BASE/src/LIP/Top/data/samples_2012.json -d /store/cmst3/user/psilva/Moriond2013_ntuples/ -o ~/work/top/2012_raw/ -c test/runAnalysis_cfg.py.templ -p "@sfMetCut=40 @ofMetCut=0 @jetPtCut=30 @applyDYweight=False @runSystematics=False @saveSummaryTree=True" -s 8nh
runLocalAnalysisOverSamples.py -e showControlDistributions -j $CMSSW_BASE/src/LIP/Top/data/syst-samples_2012.json -d /store/cmst3/user/psilva/Moriond2013_ntuples/ -o ~/work/top/2012/ -c test/runAnalysis_cfg.py.templ -p "@sfMetCut=40 @ofMetCut=0 @jetPtCut=30 @applyDYweight=False @runSystematics=False @saveSummaryTree=True" -s 8nh
runPlotter --iLumi 16689 --inDir ~/work/top/2012_raw/ --outDir ~/work/top/2012_raw/plots --json data/samples_2012.json --outFile ~/work/top/plotter_raw.root --noLog  --plotExt .pdf --showUnc
runPlotter --iLumi 16689 --inDir ~/work/top/2012/ --outDir ~/work/top/2012/plots --json data/syst-samples_2012.json --outFile ~/work/top/plotter_syst.root --noPlot

#
# FIT DY
#
fitDYcontribution --in ~/work/top/plotter_raw.root --ttrep ~/work/top/plotter_syst.root  --out Img/ --smooth


runLocalAnalysisOverSamples.py -e showControlDistributions -j $CMSSW_BASE/src/LIP/Top/data/samples_2012.json -d /store/cmst3/user/psilva/Moriond2013_ntuples/ -o ~/work/top/2012/ -c test/runAnalysis_cfg.py.templ -p "@sfMetCut=40 @ofMetCut=0 @jetPtCut=30 @applyDYweight=True @runSystematics=True @saveSummaryTree=True" -s 8nh
runPlotter --iLumi 16689 --inDir ~/work/top/2012/ --outDir ~/work/top/2012/plots --json data/samples_2012.json --outFile ~/work/top/plotter.root --noLog --plotExt .pdf --showUnc

#compute DY fit systematics
fitDYcontribution --in ~/work/top/plotter.root --ttrep ~/work/top/plotter_syst.root --out Img/ --smooth --syst

#
# FIT THE CROSS SECTION
#
fitCrossSection --in ~/work/top/plotter.root --json data/samples_2012.json --syst ~/work/top/plotter_syst.root --bins 2,3,4 --out ~/www/top/xsec > ~/www/top/xsec/result.txt
fitCrossSection --in ~/work/top/plotter.root --json data/samples_2012.json --syst ~/work/top/plotter_syst.root --bins 1 --out ~/www/top/xsec/1jet > ~/www/top/xsec/1jet/result.txt
fitCrossSection --in ~/work/top/plotter.root --json data/samples_2012.json --syst ~/work/top/plotter_syst.root --bins 2 --out ~/www/top/xsec/2jet > ~/www/top/xsec/2jet/result.txt
fitCrossSection --in ~/work/top/plotter.root --json data/samples_2012.json --syst ~/work/top/plotter_syst.root --bins 3 --out ~/www/top/xsec/3jet > ~/www/top/xsec/3jet/result.txt
fitCrossSection --in ~/work/top/plotter.root --json data/samples_2012.json --syst ~/work/top/plotter_syst.root --bins 4 --out ~/www/top/xsec/4jet > ~/www/top/xsec/4jet/result.txt


#
# Lxy MEASUREMENT
#
runPlotter --iLumi 14063 --inDir ~/work/top/2012/ --outDir ~/www/top/secvtx/madgraph --json data/samples_2012.json --noRoot --noLog --showUnc --chi2 --no2D --only jetsv --only jetmass --only lxy
runPlotter --iLumi 14063 --inDir ~/work/top/2012/ --outDir ~/www/top/secvtx/mcatnlo --json data/samples_mcatnlo_2012.json --noRoot --noLog --showUnc --chi2  --no2D --only jetsv --only jetmass --only lxy
runPlotter --iLumi 14063 --inDir ~/work/top/2012/ --outDir ~/www/top/secvtx/powheg --json data/samples_powheg_2012.json --noRoot --noLog --showUnc --chi2  --no2D --only jetsv --only jetmass --only lxy

prepareLxyDistributions --in ~/work/top/2012/plotter.root,~/work/top/2012/plotter_syst.root
FitSecVtxDistributions templ=lxydil.root

#
# CORRECT ASSIGNMENTS MEASUREMENT (INCLUDING CALIBRATION)
#
MljAnalysisCalibration --in test/results/ --json data/samples_2012.json

#
# QCD DIJET ANALYSIS
#
runLocalAnalysisOverSamples.py -e runQCDAnalysis -j $CMSSW_BASE/src/LIP/Top/data/samples_qcd_2012.json -d /store/group/phys_btag/performance/CMSSW_5_3_2_patch4/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A -o $CMSSW_BASE/src/LIP/Top/test/results -c test/runAnalysis_cfg.py.templ -p "@jetPtCut=30 @sfMetCut=40 @ofMetCut=0 @applyDYweight=False" -s 8nh
runPlotter --inDir test/results/ --outDir test/results/plots --json data/samples_qcd_2012.json --noPlot --outFile ~/work/top/plotter_qcd.root
root -b -q "bin/HFC/getQCDWeights.C+(\"~/work/top/plotter.root\",\"~/work/top/plotter_qcd.root\")"
mv QCDweights.root data/
runLocalAnalysisOverSamples.py -e runQCDAnalysis -j $CMSSW_BASE/src/LIP/Top/data/samples_qcd_2012.json -d /store/group/phys_btag/performance/CMSSW_5_3_2_patch4/MC/QCD_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A -o $CMSSW_BASE/src/LIP/Top/test/results -c test/runAnalysis_cfg.py.templ -p "@jetPtCut=30 @sfMetCut=40 @ofMetCut=0 @applyDYweight=False @weightsFile='data/QCDweights.root'" -s 8nh
runPlotter --inDir test/results/ --outDir test/results/plots --json data/samples_qcd_2012.json





#instructions below are not up to date...

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
fusermount -u store && rm -rf store
