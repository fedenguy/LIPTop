#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
CMSSW=${SCRIPTPATH}/../../src
ME=`whoami`
MYLETTER=${ME:0:1}

#configure environment
cd $CMSSW
export SCRAM_ARCH=$ARCH
eval `scram r -sh`
cd $CMSSW_BASE/src/LIP/Top

#parse the command line
input="~aalves/public/DoubleMuon-v3.root"
output="/castor/cern.ch/user/${MYLETTER}/${ME}/Dileptons/"
scheme="std"
evStart=0
evEnd=-1
for i in $*
do
  case $i in
      -src=*)
      input=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
      ;;
      -out=*)
      output=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
      ;;
      -f=*)
      evStart=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
      ;;
      -e=*)
      evEnd=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
      ;;
      -run=*)
      scheme=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
      ;;
      -help*)
      echo "wrapRunKinAnalysis [-src=input_file] [-out=out_dir] [-f=evStart] [-e=evEnd] [-run=std]"
      exit -1
      ;;
  esac
done

#generate cfg file
cfg="/tmp/${RANDOM}_cfg.py"
cat > $cfg <<EOF
import FWCore.ParameterSet.Config as cms

kinProcess = cms.PSet(
    input = cms.string("${input}"),
    output = cms.string("${output}"),
    evStart = cms.int32(${evStart}),
    evEnd = cms.int32(${evEnd}),
    dirName = cms.string("evAnalyzer/data"),
    kinScheme = cms.string("${scheme}"),
    maxTries = cms.int32(10000),
    maxJetMult = cms.int32(2),
    mw = cms.double(80.398),
    mb = cms.double(4.8),
    ptResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt"),
    etaResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_EtaResolution_AK5PF.txt"),
    phiResolFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt"),
    jesUncFileName = cms.string("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt")
    )
EOF

#run with the arguments passed
runKinAnalysis $cfg
rm $cfg
