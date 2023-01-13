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
small=( SMALL 12 108 )
apcq=( APC 40 300 )
apcmax=( APC 40 360 )

## =================================== ##
que=${small[0]}
ncpu=${small[1]}
mem=${small[2]}


qsub=mgs10_14_binning_metawrap
out="/aptmp/hiroakif/mgs10_binning_metawrap"
fasta="/aptmp/hiroakif/bin_fasta_metawrap"


## --------------------------------------------------------  ##
mkdir metagenome/${qsub}
mkdir ${out}/${qsub}

ls -d ${fasta}/* > ${out}/${qsub}/fastANIlist.txt

cat <<EOF > metagenome/${qsub}/${qsub}.qsub
#!/bin/bash
#PBS -q ${que}
#PBS -N ${qsub}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${out}/${qsub}
#PBS -j eo

## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load fastani

cd /aptmp/hiroakif/

## -------------------------------------------------------- ##

fastANI -t ${ncpu} --ql ${out}/${qsub}/fastANIlist.txt --rl ${out}/${qsub}/fastANIlist.txt -o ${out}/${qsub}/fastani_out.tsv --matrix fastANI.matrix

EOF

chmod +x metagenome/${qsub}/${qsub}.qsub
qsub metagenome/${qsub}/${qsub}.qsub

echo " "
echo "********* INFO: ${qsub}.qsub was generated and started! *********"

## -------------------------------------------------------- ##
