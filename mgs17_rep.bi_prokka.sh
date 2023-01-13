################################################
## 2020 11 08.
## For the meta-genome analysis
##
## mgs_10_21 Finding ortholog
##  on super computer system
##
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

quesys=${small}

que=${quesys[0]}
ncpu=${quesys[1]}
mem=${quesys[2]}


bins="/aptmp/hiroakif/mgs10_binning_metawrap"
fasta="/aptmp/hiroakif/representive_bin_fasta_metawrap"
qsub=mgs10_binning_metawrap
out=prokka_rep.bins

## =================================== ##
mkdir ${bins}/${out}

for i in `ls ${fasta} | grep fa`
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
module load prokka

cd /aptmp/hiroakif/

## -------------------------------------------------------- ##

mkdir ${bins}/${out}/${i}

prokka --force ${fasta}/${i} -o ${bins}/${out}/${i} --metagenome --cpus 40

 
EOF
	chmod +x metagenome/${qsub}/${i}_${qsub}.qsub
	qsub metagenome/${qsub}/${i}_${qsub}.qsub

	echo " "
	echo "********* INFO: ${i}_${qsub}.qsub was generated and started! *********"

done
