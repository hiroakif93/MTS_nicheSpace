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

bins="/aptmp/hiroakif/representive_bin_fasta_metawrap"
out="/aptmp/hiroakif/mgs10_binning_metawrap/mgs10_16_representive_bin_coverage"
qsub=mgs10_16_mapping_metawrap
## --------------------------------------------------------  ##

mkdir ${out}
mkdir metagenome/${qsub}

for i in `ls ${bins}`
	do
	
cat <<EOF > metagenome/${qsub}/${i}.qsub
#!/bin/bash
#PBS -q ${que}
#PBS -N ${i}_${qsub}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${qsub}
#PBS -j eo
## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load coverm R/3.6.1

cd /aptmp/hiroakif/

bins=${bins}
out=${out}
## -------------------------------------------------------- ##

EOF

cat <<'EOF' >> metagenome/${qsub}/${i}.qsub

for l in `ls mgs10_binning_metawrap/fasqsf | cut -b 1-7 | uniq`
do

coverm genome --genome-fasta-directory ${bins} -x fa -m mean relative_abundance rpkm -1 mgs10_binning_metawrap/fasqsf/${l}_1.fastq -2 mgs10_binning_metawrap/fasqsr/${l}_2.fastq -t 40 > ${out}/${l}_mapping.tsv

done		
EOF

	chmod +x metagenome/${qsub}/${i}.qsub
	qsub metagenome/${qsub}/${i}.qsub

	echo " "
	echo "********* INFO: ${i}_${out}.qsub was generated and started! *********"

done
## -------------------------------------------------------- ##
