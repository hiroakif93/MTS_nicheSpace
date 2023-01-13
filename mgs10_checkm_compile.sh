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
apcq=( APC 30 10 )
apcmax=( APC 40 360 )

## =================================== ##
que=${apcq[0]}
ncpu=${apcq[1]}
mem=${apcq[2]}

bins="/aptmp/hiroakif/mgs10_binning_metawrap"
out="/aptmp/hiroakif/bin_fasta_metawrap"
qsub=mgs10_binning_metawrap
## --------------------------------------------------------  ##

echo -e "sample\tbin_id\tcompleteness\tcontamination\tstrain_heterogeneity\tgenome_size\tN50_contig\tGC.content\tcoding_density\tpredicted_genes" > ${bins}/bin_summary.txt

for i in `ls ${bins} | grep MTS`
do

	sed -e '1d' ${bins}/${i}/checkm_qa.txt | awk -F "\t" -v foo=${i} 'BEGIN{ OFS = "\t" } {print foo,$1,$6,$7,$8,$9,$14,$19,$21,$23}' >> ${bins}/bin_summary.txt 
	
 
done

