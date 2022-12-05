############################################################################
####
#### R script for Fujita (2019)
####
#### Metagenome analysis
#### 2020.11.12 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### setwd('../')
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c("seqinr") )

## -- Load table
gene.fasta <- read.fasta('Table/ShotgunMetagenome/representive_MAG_protain.faa')
mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds')
repmag  <- mag[unique(mag$represent), ] 

## ======================================================================== ##

## -- Get the gene names
gene.list <- sapply( sapply(strsplit(sapply(gene.fasta, getAnnot), ' '), '[', -1), paste, collapse='.' )
unique(paste(sapply( strsplit(names(gene.list), '_'), '[', 1), sapply( strsplit(names(gene.list), '_'), '[', 2))) 
 
## -- Count gene family each representive bin
gene.table <- lapply(rownames(repmag), function(x) table(gene.list[grep(x, names(gene.list))]))
names(gene.table) <- rownames(repmag)
lapply(gene.table, head)

gene.list2 <- unique(do.call(c, sapply(gene.table, names)))
gene.count.mat <- matrix(0, nrow=length(gene.table), ncol=length(gene.list2),
						 dimnames=list(rownames(repmag), gene.list2))

for(i in names(gene.table)){ gene.count.mat[i, names(gene.table[[i]])] <- gene.table[[i]] }	
			 
## ======================================================================== ##
gene.count.mat2 <- gene.count.mat
gene.count.mat2 <- gene.count.mat2[,-which(colnames(gene.count.mat2)=='hypothetical.protein')]

gene.count.mat2[gene.count.mat2>0] <- 1
freq <- colSums(gene.count.mat2)/nrow(gene.count.mat2)
genemat <- gene.count.mat2

## ======================================================================== ##

saveRDS(gene.table, 'Table/04_03_gene_table_represent.rds')
saveRDS(gene.count.mat, 'Table/04_03_gene_count.matrix_represent.rds')
saveRDS(genemat, 'Table/04_03_gene_binary.matrix_represent.rds')

## ======================================================================== ##
