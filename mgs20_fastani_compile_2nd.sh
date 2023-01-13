################################################
## 2020 12 24.
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
small=( SMALL 12 30 )
apcq=( APC 12 40 )
apcmax=( APC 40 360 )

## =================================== ##
que=${apcq[0]}
ncpu=${apcq[1]}
mem=${apcq[2]}


out=mgs10_15
## =================================== ##
mkdir metagenome/${out}

## =================================== ##

cat <<'RRR' > metagenome/${out}/mgs10_15_extract_representive_bins.R
	
	setwd('/aptmp/hiroakif/')

	
	
df <- read.table('mgs10_binning_metawrap/mgs10_14_binning_metawrap/fastani_out.tsv', header=F)
bin.info <- read.table('mgs10_binning_metawrap/bin_summary.txt', header=T)
taxa.info <- read.table('mgs10_binning_metawrap/taxa_list.txt', header=T, sep='\t')

taxa.info$id <- apply(taxa.info[,1:2], 1, paste, collapse='_')
bin.info$id <- apply(bin.info[,1:2], 1, paste, collapse='_')
bin.info$name <- paste(bin.info$id, 'fa', sep='.')


df[,1] <- sapply( strsplit(as.character(df[,1]), '/'), '[', 5 )
df[,2] <- sapply( strsplit(as.character(df[,2]), '/'), '[', 5 )

df2 <- df3 <- df[ df[,3] > 99, ]

dupli <- c()
for(i in unique(df2[,1]) ){ # i=unique(df2[,1])[88] ; i='MTS-060_bin.1.orig.fa'
	
	bin <- df3[ df3[,1]==i, ]
	
	bin.ex <- bin[bin[,1] != bin[,2],]
	
	if( nrow(bin.ex)>0 & sum(df3[, 1]%in%bin.ex[,2])>0 ){
		dupli <- c(dupli, which(df2[, 1]%in%bin.ex[,2]) )
		
		df3 <- df3[-which(df3[, 1]%in%bin.ex[,2]),] ; print(dim(df3)) ;print(i)
	}
}

## ------------------------------------------------------ ##
	
rep.bin.list <- rep.bin.mat <- set.info.summary <- c()
for(i in unique(df3[,1]) ){ # i=unique(df3[,1])[2]
	
	## ------------------------------------------------------ ##
	
	bin.set <- as.vector( df3[df3[,1]==i, 2] )
	
	set.info <- bin.info[bin.info$name %in% bin.set, ]
	
	represent.bin <- set.info[set.info$completeness==max(set.info$completeness), ]
	represent.bin2 <-  represent.bin[represent.bin$N50_contig==max(represent.bin$N50_contig), ]
	
	rep.bin.list <- rbind(rep.bin.list, represent.bin2[, 'name'] )
	
	## ------------------------------------------------------ ##
	
	set.info2 <- as.data.frame(set.info[,-ncol(set.info)])
	representive <-  paste( as.matrix(represent.bin2)[,1:2], collapse='_')
	
	set.info3 <- cbind(set.info2[,1:2], representive=paste( as.matrix(represent.bin2)[,1:2], collapse='_'),
						taxa.info[ which(taxa.info$id %in% set.info2$id ), -c(1:2)],
					    set.info2[,-c(1:2, ncol(set.info2))])
						
	set.info.summary <- rbind( set.info.summary,  set.info3)
	
	## ------------------------------------------------------ ##
	
	rep.taxa <- unique(taxa.info[ which(taxa.info$id %in% represent.bin$id ), -c(1:2, ncol(taxa.info) ) ])
	
    mean <- colMeans(represent.bin[,-c(1:2,-1:0+ncol(represent.bin) )]) 
	sd=apply(represent.bin[,-c(1:2,-1:0+ncol(represent.bin) )], 2, sd)
	
	rep.bin.mat <- rbind(rep.bin.mat, c(as.matrix(represent.bin2[,c(1:2)]), id=represent.bin2$id, as.matrix(rep.taxa), mean=mean, sd=sd ) )
	
	## ------------------------------------------------------ ##
	
}

set.info.summary2 <- set.info.summary[order(set.info.summary$sample), ]

colnames(rep.bin.mat)[1:2] <- c('sample', 'bin_id')
colnames(rep.bin.mat)[4:10] <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

write.table(rep.bin.list, 'mgs10_binning_metawrap/representive_bin.txt', quote=F, row.names=F, col.names=F)
write.table(set.info.summary2, 'mgs10_binning_metawrap/bin_summary2.txt', quote=F, row.names=F, col.names=T, sep='\t')
write.table(rep.bin.mat, 'mgs10_binning_metawrap/representive_bin_summary.txt', quote=F, row.names=F, col.names=T, sep='\t')

## ------------------------------------------------------ ##
dir.create('representive_bin_fasta_metawrap')
dir.create('mgs10_binning_metawrap/representive_bin_fasta_metawrap')

files <- list.files('bin_fasta_metawrap')[list.files('bin_fasta_metawrap') %in% rep.bin.list]

sapply(files, function(x) file.copy( sprintf('bin_fasta_metawrap/%s',x), 'mgs10_binning_metawrap/representive_bin_fasta_metawrap' ) )

## ------------------------------------------------------ ##
	
RRR

## =================================== ##

cat <<EOF > metagenome/${out}/${out}.qsub
#!/bin/bash
#PBS -q ${que}
#PBS -N ${i}_${out}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${out}
#PBS -j eo

## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load rnammer R/3.6.1

cd /user4/kyoto1/hiroakif

## -------------------------------------------------------- ##

## -------------------------------------------------------- ##

Rscript --vanilla metagenome/${out}/mgs10_15_extract_representive_bins.R

echo "Finish R part"
## -------------------------------------------------------- ##
EOF



	chmod +x metagenome/${out}/${out}.qsub
	qsub metagenome/${out}/${out}.qsub

	echo " "
	echo "********* INFO: ${out}.qsub was generated and started! *********"


## -------------------------------------------------------- ##
