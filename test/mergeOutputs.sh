#!/bin/bash

tags=(WWtoAnything_Spring11 WZtoAnything_Spring11 ZZtoAnything_Spring11 TToBLNu_tW-channel_Spring11 TToBLNu_t-channel_Spring11 TToBLNu_s-channel_Spring11 TTJets_madgraph_Spring11 WJetsToLNu_Spring11 DYToEE_M-20_Spring11 DYToMuMu_M-20_Spring11 DYToTauTau_M-20_Spring11 QCD_Pt-30to50_Spring11 QCD_Pt-50to80_Spring11 QCD_Pt-80to120_Spring11 QCD_Pt-120to170_Spring11 QCD_Pt-170to300_Spring11 QCD_Pt-300to470_Spring11 QCD_Pt-470to600_Spring11 QCD_Pt-600to800_Spring11 QCD_Pt-800to1000_Spring11 DoubleMuon-v3 DoubleMuon-v5 DoubleMuon-v6 DoubleMuon-v7 DoubleMuon-v8 DoubleElectron-v3 DoubleElectron-v5 DoubleElectron-v6 DoubleElectron-v7 DoubleElectron-v8  MuEG-v3 MuEG-v5 MuEG-v6 MuEG-v7 MuEG-v8)

for i in ${tags[@]}; do
    echo "****** $i *******"
    hadd -f ${i}.root ${i}_*.root
    rm ${i}_*.root
done

