#!/bin/bash
m=(161.5 163.5 166.5 169.5 172.5 175.5 178.5 181.5 184.5)
#m=(172.5)
s=(std)
#s=(lesup lesdown)
#s=(baresignal jer jesup jesdown matchingup matchingdown scaleup scaledown puup pudown)
#s=(std dyup dydown effbup effbdown effqup effqdown)
for i in ${m[@]}; do
    for j in ${s[@]}; do
	TopMassCalibration  --in store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par bin/Mass/MassParFile_exclusive.txt &
#	outdir="${i}_${j}"
#	mkdir -p ${outdir}
#	cd ${outdir}
#	TopMassCalibration ../EventSummaries.root ../bin/Mass/MassParFile_exclusive_signalonly.txt ${j} ${i} ./ &
#	cd -
    done
done