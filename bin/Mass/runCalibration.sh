#!/bin/bash

#m=(161.5 163.5 166.5 169.5 172.5 175.5 178.5 181.5 184.5)
#s=(std)

#m=(172.5)
#s=(baresignal jesup jesdown matchingup matchingdown scaleup scaledown puup pudown mcgen)

m=(172.5)
s=(pushift-5 pushift-4 pushift-3 pushift-2 pushift-1 pushift+1 pushift+2 pushift+3 pushift+4 pushift+5)

#m=(172.5)
#s=(dyup dydown effbup effbdown effqup effqdown)

#m=(172.5)
#s=(lesup lesdown jer)

for i in ${m[@]}; do
    for j in ${s[@]}; do
	#TopMassCalibration  --in ../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par MassParFile_inclusive.txt &
	#TopMassCalibration  --in ../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par MassParFile_exclusive.txt &
	outdir="${i}_${j}"
	mkdir -p ${outdir}
	cd ${outdir}
	TopMassCalibration --in ../../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par ../MassParFile_exclusive_signalonly.txt --syst ${j} &
	cd -
    done
done