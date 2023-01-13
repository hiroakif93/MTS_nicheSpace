################################################
## 2020 11 08.
## For the meta-genome analysis
##
## mgs_05 checking binnig quality by using checkm  
##  on super computer system
##
## https://ecogenomics.github.io/CheckM/
##
################################################

#!/bin/bash

## =================================== ##
## Bash que list
quick=( QUICK 4 72 )
small=( SMALL 12 108 )
apcq=( APC 20 50 )
apcmax=( APC 40 360 )

## =================================== ##
que=${apcq[0]}
ncpu=${apcq[1]}
mem=${apcq[2]}

bins="/aptmp/hiroakif/mgs10_binning_metawrap/mgs10_16_representive_bin_coverage"
qsub=mgs10_16_binning_metawrap
## --------------------------------------------------------  ##

echo -e "sample\tbin_id\tmean\trelative_abndunce\tRPKM" > ${bins}/bin_abundacne.txt

for i in `ls ${bins} | grep MTS`
do

	sed -e '1d' ${bins}/${i} | awk -F "\t" -v foo=`echo ${i} | cut -b 1-7` 'BEGIN{ OFS = "\t" } {print foo,$1,$2,$3,$4}' >> ${bins}/bin_abundacne.txt
	
 
done