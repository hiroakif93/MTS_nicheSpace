############################################################################
####
#### R script for Fujita (2019)
####
#### Metagenome analysis
#### 2020.11.12 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_nichespace/')
#### setwd('../')
############################################################################

## -- Loading Function and Library
source('functions/functions.R')

load.lib(c('seqinr', 'dada2'))

# -- Load data table
seqtab <- readRDS("Table/00_seqtab_rmchimera.rds")
seq.fna <- read.fasta("Table/01_Fasta_16S.fasta")
asvlist <- readRDS('Table/Taxa_list.rds')

dlist <- readRDS('Table/matrixList.rds')
ts <- dlist[["Water/Medium-A"]]
ts <- ts[rownames(sml), colSums(ts[rownames(sml), ])>0]

sml <- readRDS('Table/sample_info.rds') ; sml <- sml[which(sml $treat1=="Water/Medium-A" & sml $replicate.id==5), ]


############################################################################

fna.mat <- as.data.frame(cbind(names(seq.fna), do.call(rbind, lapply(seq.fna, paste, collapse='')))) ; rownames(fna.mat) <- names(seq.fna)
asvlist$sequence <- fna.mat[rownames(asvlist),2]

## ---------------------------------------------------------- ##
reference <- "Table/bac120_ssu_reps_r95_reforamt_append_std.fasta"
taxa <- assignTaxonomy(seqtab, reference, multithread=TRUE)

taxa.print <- as.data.frame(taxa)

gbtk <- as.matrix(taxa.print[toupper(asvlist$sequence), ])
gbtk[is.na(gbtk)] <- 'unidentified'

sum <- cbind(asvlist[,-ncol(asvlist)], gbtk=gbtk, asvlist[,ncol(asvlist)])
head(sum[,-ncol(sum)])
saveRDS(sum[,-ncol(sum)], 'Table/04_Taxa_list.rds')

