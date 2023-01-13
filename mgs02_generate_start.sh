
################################################
## 2020 10 21.
## For the meta-genome analysis
##
## mgs_02 genome assembly by using SPAdes 
## 			and take more than 2500 bp scaffolds by using seqkit
##  on super computer system
##
## https://github.com/ablab/spades
## https://bioinf.shenwei.me/seqkit/usage/
##
################################################

#!/bin/bash
que="APC"
ncpu="40"
mem="300"

load="/aptmp/hiroakif/mgs01_fastp"
out=mgs02_genome_assembly

mkdir metagenome/${out}

for i in `ls ${load} | grep fastq.gz | cut -b 1-7 | uniq`
	do

	fwd=`ls ${load} | grep $i | grep R1`
	rvs=`ls ${load} | grep $i | grep R2`

	cat <<EOF > metagenome/${out}/${i}_${out}.qsub
	
#!/bin/csh
#PBS -q ${que}
#PBS -N ${i}_${out}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${out}
#PBS -j eo

source /etc/profile.d/modules.sh
module load SPAdes seqkit

cd /aptmp/hiroakif
		
mkdir ${out}/${i}


printf '## ----- Start genome assembly of %s and %s. ----- ##\n' $fwd $rvs
start_time=`date +%s`

spades.py --meta -1 ${load}/$fwd -2 ${load}/$rvs -o ${out}/${i} -k 21,33,55,77,99,121 -m 350 -t 40

end_time=`date +%s`
run_time=$((end_time - start_time))

seqkit seq -g -m 2500 ${out}/${i}/scaffolds.fasta > ${out}/${i}_scaffolds_m2500.fasta

printf '## ---------    Done. It takes %s sec.    --------- ##\n' $run_time
	
EOF


	chmod +x metagenome/${out}/${i}_${out}.qsub
	qsub metagenome/${out}/${i}_${out}.qsub

	echo " "
	echo "********* INFO: ${i}_${out}.qsub was generated and started! *********"

done

