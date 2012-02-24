#!/bin/bash

# This is for quick redo of a whole directory of histograms. Just edit it for the correct regexp
# (or use one of the ones I did already develop) and run it from the plots directory
# Pietro

for i in `ls *.C`; do

    # Lumi replacement
    sed -i -e 's/L=2.2/L=2.3/g' $i

    # Resaving pdfs and pngs
    # WARNING: this substitution you should run only once, because it will match everytime the SetSelected,
    # thus accumulating a cascade of SaveAs. It is not destructive (just rewrites with the same histo), but takes more time.
    j=`echo $i | sed 's/\.C//'`
    sed -i -e 's/\(.*[a-zA-Z0-9]\)\->SetSelected(\1);/\1\->SetSelected(\1);\1\->SaveAs("'$j'\.png");\1\->SaveAs("'$j'\.pdf");/' $i
    root -b -q "$i"
done
