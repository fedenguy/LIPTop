#!/bin/bash

m=(161.5 163.5 166.5 169.5 172.5 175.5 178.5 181.5 184.5)
s=(std)

#m=(172.5)
#s=(baresignal jesup jesdown matchingup matchingdown scaleup scaledown puup pudown mcgen)

#m=(172.5)
#s=(pushift-5 pushift-4 pushift-3 pushift-2 pushift-1 pushift+1 pushift+2 pushift+3 pushift+4 pushift+5)

#m=(172.5)
#s=(pdf0 pdf1 pdf2 pdf3 pdf4 pdf5 pdf6 pdf7 pdf8 pdf9 pdf10 pdf11 pdf12 pdf13 pdf14 pdf15 pdf16 pdf17 pdf18 pdf19 pdf20 pdf21 pdf22 pdf23 pdf24 pdf25 pdf26 pdf27 pdf28 pdf29 pdf30 pdf31 pdf32 pdf33 pdf34 pdf35 pdf36 pdf37 pdf38 pdf39 pdf40 pdf41 pdf42 pdf43 pdf44)

#m=(172.5)
#s=(dyup dydown effbup effbdown effqup effqdown)

#m=(172.5)
#s=(lesup lesdown jer)

for i in ${m[@]}; do
    for j in ${s[@]}; do
	#TopMassCalibration  --in ../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par MassParFile_inclusive.txt &
	#TopMassCalibration  --in ../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par MassParFile_exclusive.txt &
	TopMassCalibration  --in ../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 500 --par MassParFile_dilepton_exclusive.txt &
        #outdir="${i}_${j}"
	#mkdir -p ${outdir}
	#cd ${outdir}
	#TopMassCalibration --in ../../../store/Top/ntuples_14Jan/merged/summary/EventSummaries.root --calib ${i} --out ./ --npe 100 --par ../MassParFile_exclusive_signalonly.txt --syst ${j} &
	#TopMassCalibration --in ../../../pdfs/EventSummaries.root --par ../MassParFile_exclusive_signalonly.txt --out ./ --npe 50 --syst ${j} &
	#cd -
    done
done