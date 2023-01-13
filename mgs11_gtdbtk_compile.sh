################################################
## 2020 11 08.
## For the meta-genome analysis
##
##  mgs_06 taxonomy annotation by using GTDBtk
##
## Extract higher quality bin (completeness > 90%, contamination < 5)
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
que=${quick[0]}
ncpu=${quick[1]}
mem=${quick[2]}


bins="/aptmp/hiroakif/mgs10_binning_metawrap"
fasta="/aptmp/hiroakif/bin_fasta_metawrap"
qsub=mgs10_binning_metawrap

## =================================== ##
mkdir ${bins}/prokka

for i in `ls ${bins} | grep MTS`
do
	
	sed -e '1d' ${bins}/${i}/gtdbtk/gtdbtk.bac120.summary.tsv | awk -F "\t" -v foo=${i} 'BEGIN{ OFS = "\t" } {print foo,$1,$2}' >> ${bins}/gtdbtk.bac120.summary.txt

done 

## =================================== ##

cat <<EOF > metagenome/${qsub}/${qsub}_13_compile.qsub
#!/bin/bash
#PBS -q ${que}
#PBS -N ${i}_${qsub}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${qsub}
#PBS -j eo

## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load R/3.6.1

cd /user4/kyoto1/hiroakif

cat <<'ROF' > metagenome/${qsub}/${qsub}_13_compile.R
taxa_sum <- c()
bar=c()

	df <- read.table('gtdbtk.bac120.summary.txt', header=F, sep='\t')
	
	splt <- lapply(strsplit(as.character(df[,3]), ';'), function(x) sapply( strsplit(x, '__'), '[', 2) )
	
	taxa <- cbind(sample=as.character(df[,1]), bin_id=as.character(df[,2]), do.call(rbind, splt))
	taxa[is.na(taxa)] <- 'unidentified'
	
	colnames(taxa) <- c( 'sample', 'bin_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' )
	
	taxa_sum <- rbind(taxa_sum, taxa)


write.table(taxa_sum, 'taxa_list.txt', sep='\t', quote=F, row.names=F)

ROF

cd ${bins}
Rscript --vanilla /user4/kyoto1/hiroakif/metagenome/${qsub}/${qsub}_13_compile.R


EOF

chmod +x metagenome/${qsub}/${qsub}_13_compile.qsub
qsub metagenome/${qsub}/${qsub}_13_compile.qsub

echo " "
echo "********* INFO: metagenome/${qsub}/${qsub}.qsub was generated and started! *********"
