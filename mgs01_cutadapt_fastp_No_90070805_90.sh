
################################################
## 2020 10 12.
## For the meta-genome analysis
##
## For MTS-001 ~ MTS-020, MTS-030 ~ MTS-080 and MTS-100 ~ MTS-110
## (Cutting adaptor 5') 
##
## mgs_01 Cutting adaptor by using cutadapt
## 	  Filtering low quality and checking fastaq by using fastp
## 	  on super computer system in kyoto univ. 
## 
## https://cutadapt.readthedocs.io/en/stable/installation.html
## https://github.com/OpenGene/fastp
##
################################################
## -- Install
# python3 -m pip install --user --upgrade cutadapt
# export PATH=$PATH:$HOME/.local/bin
# source ~/.bash_profile
# echo $PATH

#!/bin/bash

out1="mgs01_cutadapt"
out2="mgs01_fastp"

fa="No_90070805_90"

mkdir ${out1} ; mkdir ${out2}

for i in `ls $fa | grep fastq | cut -b 1-7 | uniq`
do
	fwd=`ls $fa | grep $i | grep R1`
	rvs=`ls $fa | grep $i | grep R2`
	
	printf '\n\n##	==============================================================  ##\n'
	printf '##	==========		Start filtering of %s and %s.		==========	##\n' $fwd $rvs
	start_time=`date +%s`
	
	cutadapt -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -A AAGTCGGATCGTAGCCATGTCGTTCTGTGAGC -o ${out1}/${i}_R1.fastq.gz -p ${out1}/${i}_R2.fastq.gz $fa/$fwd $fa/$rvs -j 16 -m 50 --max-n 0
	/user4/kyoto1/hiroakif/fastp -A -q 20 -w 16 -3 -W 6 -M 20 -i ${out1}/${i}_R1.fastq.gz -I ${out1}/${i}_R2.fastq.gz -o ${out2}/${i}_R1.fastq.gz -O ${out2}/${i}_R2.fastq.gz -h ${out2}/report_${i}.html

	end_time=`date +%s`
	run_time=$((end_time - start_time))
	
	printf '\n##	==========		%s Done. It takes %s sec.		==========	##\n' $i $run_time
	printf '##	======================================================	##\n\n'
	
done



#out1="cutadapt"
#out2="fastp"

#mkdir $out1
#mkdir $out2
	
#	cutadapt -a CTGTCTCTTATACACATCT -A GACAGAGAATATGTGTAGA -o $out1/outF_MTS-01.fastq.gz -p $out1/outR_MTS-01.fastq.gz MTS-01_S62_R1_001.fastq.gz MTS-01_S62_R2_001.fastq.gz -j 8
#	fastp -A -q 20 -i $out1/outF_MTS-01.fastq.gz -I $out1/outR_MTS-01.fastq.gz -o $out2/#F_MTS-01.fq.gz -O $out2/R_MTS-01.fq.gz -h $out2/report_MTS-01.html -w 8

#out1="cutadapt2"
#out2="fastp2"

#mkdir $out1
#mkdir $out2

#	cutadapt -a CTGTCTCTTATACACATCT -g GACAGAGAATATGTGTAGA -A CTGTCTCTTATACACATCT -G GACAGAGAATATGTGTAGA -o $out1/outF_MTS-01.fastq.gz -p $out1/outR_MTS-01.fastq.gz MTS-01_S62_R1_001.fastq.gz MTS-01_S62_R2_001.fastq.gz -j 8
#	fastp -A -q 20 -l 75 -i $out1/outF_MTS-01.fastq.gz -I $out1/outR_MTS-01.fastq.gz -o $out2/#F_MTS-01.fq.gz -O $out2/R_MTS-01.fq.gz -h $out2/report_MTS-01.html -w 8


# 	trimmomatic PE -threads 40 -phred33 -trimlog mgs01/trimlog.txt mgs01/fastp/MTS-001_R1.fastq.gz mgs01/fastp/MTS-001_R2.fastq.gz MTS-001_R1_paired.fq.gz MTS-001_R1_unpaired.fq.gz MTS-001_R2_paired.fq.gz MTS-001_R2_unpaired.fq.gz LEADING:20 HEADCROP:10 MINLEN:75
