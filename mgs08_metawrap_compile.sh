################################################
## 2020 11 08.
## For the meta-genome analysis
##
## mgs_04 binning by using Maxbin 
##  on super computer system
##
## https://sourceforge.net/projects/maxbin/files/
##
################################################

#!/bin/bash

## =================================== ##
## Bash que list
quick=( QUICK 4 72 )
small=( SMALL 12 30 )
apcq=( APC 30 40 )
apcmax=( APC 40 360 )

## =================================== ##
que=${quick[0]}
ncpu=${quick[1]}
mem=${quick[2]}

bins="/aptmp/hiroakif/mgs10_binning_metawrap"
out="/aptmp/hiroakif/bin_fasta_metawrap"

## --------------------------------------------------------  ##
mkdir $out

cd $bins

for i in `ls | grep M`
do
	
	for l in `ls ${i}/reassemble/reassembled_bins`
	do
		
		renamed="${i}_${l}"
		cp ${i}/reassemble/reassembled_bins/${l} $out/${renamed}

	done
		
done
## -------------------------------------------------------- ##
