#!/bin/bash

mycastordir=/castor/cern.ch/user/p/psilva/Dileptons
results=(`rfdir ${mycastordir} | awk '{print $9}'`)
echo "${#results[@]} results are available: start copying"
for r in ${results[@]}; do
    mkdir -p ${r}/std
    cd ${r}/std
    files=(`rfdir ${mycastordir}/${r} | awk '{print $9}'`)
    for f in ${files[@]}; do
	rfcp ${mycastordir}/${r}/$f ./
    done
    cd -
    echo "   - ${#files[@]} files for $r"
done
echo "All done"