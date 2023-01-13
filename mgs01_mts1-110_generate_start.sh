
################################################
## 2020 11 02.
## For the meta-genome analysis
##
## For MTS-1 ~110 excepting MTS24 and MTS90
##
## mgs_01 Filtering low quality and checking fastaq by using fastp
##	on super computer system in kyoto univ. 
## https://cutadapt.readthedocs.io/en/stable/installation.html
## https://github.com/OpenGene/fastp
##
################################################

#!/bin/bash
que="APC"
ncpu="40"
mem="200"


load="/aptmp/hiroakif/No_90082717_1-110"
out=mgs01_cut_and_filt
		
out1="mgs01_cutadapt"
out2="mgs01_fastp"

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
module load cutadapt

cd /aptmp/hiroakif

mkdir ${out1}
mkdir ${out2}

printf '## ----- Start genome assembly of %s and %s. ----- ##\n' $fwd $rvs
start_time=`date +%s`

cutadapt -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b CTGTCTCTTATACACATCT TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CTGTCTCTTATACACATCTGACGCTGCCGACGA AGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B AGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG CTGTCTCTTATACACATCTCCGAGCCCACGAGAC CTGTCTCTTATACACATCT TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CTGTCTCTTATACACATCTGACGCTGCCGACGA -o ${out1}/${i}_R1.fastq.gz -p ${out1}/${i}_R2.fastq.gz ${load}/$fwd ${load}/$rvs -j ${ncpu} -m 50 --max-n 0
/user4/kyoto1/hiroakif/fastp -q 20 -w ${ncpu} -3 -W 6 -M 20 -i ${out1}/${i}_R1.fastq.gz -I ${out1}/${i}_R2.fastq.gz -o ${out2}/${i}_R1.fastq.gz -O ${out2}/${i}_R2.fastq.gz -h ${out2}/report_${i}.html

end_time=`date +%s`
run_time=$((end_time - start_time))
	
printf '## ---------    Done. It takes %s sec.    --------- ##\n' $run_time
	
EOF


	chmod +x metagenome/${out}/${i}_${out}.qsub
	qsub metagenome/${out}/${i}_${out}.qsub

	echo " "
	echo "********* INFO: ${i}_${out}.qsub was generated and started! *********"

done

