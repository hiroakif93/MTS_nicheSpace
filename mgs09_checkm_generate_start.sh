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

for i in `ls ${bins} | grep MTS`
do

cat <<EOF > metagenome/${qsub}/${i}_${qsub}.qsub
#!/bin/bash
#PBS -q ${que}
#PBS -N ${i}_${qsub}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${qsub}
#PBS -j eo

## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load CheckM

cd /aptmp/hiroakif/

## -------------------------------------------------------- ##


checkm qa --tab_table -o 2 -f ${bins}/${i}/checkm_qa.txt ${bins}/${i}/reassemble/reassembled_bins.checkm/lineage.ms ${bins}/${i}/reassemble/reassembled_bins.checkm

 
EOF

	chmod +x metagenome/${qsub}/${i}_${qsub}.qsub
	qsub metagenome/${qsub}/${i}_${qsub}.qsub

	echo " "
	echo "********* INFO: ${i}_${qsub}.qsub was generated and started! *********"


## -------------------------------------------------------- ##

done