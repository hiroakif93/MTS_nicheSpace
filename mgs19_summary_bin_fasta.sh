

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


#!/bin/bash
#PBS -q ${que}
#PBS -N $${qsub}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${qsub}
#PBS -j eo


# Run on local pc

cd /aptmp/hiroakif/representive_bin_fasta_metawrap

for i in `ls | grep .fa$`
do

	id=`echo ${i} |  awk -F ".fa$" '{print $1}' `

 	sed -e "s/>/>${id}_/g" ${i} > rename_${id}.fa
 		
 	 	
done
 
cat rename_* > ../representive_bin_summary.fasta