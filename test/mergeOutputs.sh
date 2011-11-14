#!/bin/bash

inputdir=$1
outputdir=${inputdir}/merged
mkdir -p ${outputdir}

tags=(WW ZZ WZ 
    WJetsToLNu DYJetsToLL DYJetsToEE DYJetsToMuMu
    TTJets TTJets_matchingup TTJets_matchingdown TTJets_scaleup TTJets_scaledown
    SingleTbar_tW SingleTbar_t SingleTbar_s SingleT_tW SingleT_t SingleT_s SingleTbar_tW_DS SingleT_tW_DS
    DoubleElectronMay10ReReco DoubleMuMay10ReReco MuEGMay10ReReco
    DoubleElectronPromptRecov4 DoubleMuPromptRecov4 MuEGPromptRecov4
    DoubleElectron05AugReReco DoubleMu05AugReReco MuEG05AugReReco
    DoubleElectronPromptRecov6
    DoubleMuPromptRecov6_172620_173244 MuEGPromptRecov6_172620_173244
    DoubleMuPromptRecov6_173380-173692 MuEGPromptRecov6_173380-173692
)

for i in ${tags[@]}; do
    echo "****** $i *******"
    hadd -f ${outputdir}/${i}.root ${inputdir}/${i}_*.root
#    rm ${i}_*.root
done


cd ${outputdir}
ln -s TTJets.root TTJets_signal.root

