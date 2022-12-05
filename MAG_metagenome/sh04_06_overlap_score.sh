################################################
## 2021 01 23.
#### CCM
##  on super computer system
##
################################################

#!/bin/bash

## =================================== ##
## Bash que list
quick=( QUICK 1 18 )
small=( SMALL 12 108 )
apcq=( APC 12 108 )
apcmax=( APC 40 360 )
sdf=( SDF 60 700 )

## =================================== ##
que=${apcq[0]}
ncpu=${apcq[1]}
mem=${apcq[2]}

## --------------------------------------------------------  ##
jobplace=MTS_metagenome/MAG_metagenome/

## --------------------------------------------------------  ##


log=${jobplace}/scl_logs/${scriptno}
qname=${scriptno}

for i in `seq 1 13`
do
scriptno=0406_${i}
qname=${scriptno}


mkdir -p ${jobplace}/scl_logs
mkdir -p ${jobplace}/scl_logs/${scriptno}

log=${jobplace}/scl_logs/${scriptno}

RscriptName=${log}/${qname}_${i}.R
qsubSript=${log}/${qname}_${i}.qsub

cat <<EOF > ${RscriptName}
############################################################################
####
#### R script for Fujita (2019)
####
#### Metagenome analysis
#### 2020.11.12 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_metagenome/')
#### setwd('../')
############################################################################

## -- For super computer system

# -- Load library and functions
#　lapply(lib, function(x){ install.packages(lib, lib=lib.at)})
lib <-  c( 'vegan', 'bipartite', 'maxnodf')
lib.at <- '/user4/kyoto1/hiroakif/R/x86_64-pc-linux-gnu-library/4.1'

invisible(lapply(lib, function(x){ library(package = x, character.only=TRUE, lib=lib.at)}) )
invisible(sapply(lib, function(x) cat(sprintf('%s %s\n',x, packageVersion(x, lib=lib.at)) ) ))

setwd('/user4/kyoto1/hiroakif/aptmp/MTS_metagenome')
############################################################################

source('functions/functions.R')

# -- make directory to save results
dir <- make.dir('MAG_metagenome/Niche_overlap')

# -- Load data table
gene.count.mat <- readRDS('Table/04_03_gene_count.matrix_represent.rds')
mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds') 
mag[mag ==''] <- 'Unidentified'

repmag <- mag[unique(mag[,'represent']),]
ts <- readRDS('Table/04_MAG_timeSeries.rds')
ts <- ts/rowSums(ts)
cols <- readRDS('Table/04_color_palette.rds')


############################################################################

rand=10

require(doParallel)
cluster = makeCluster(${ncpu}, "FORK")
registerDoParallel(cluster)

day=${i}
EOF

cat <<'EOF' >> ${RscriptName}
############################################################################

freq <- colSums(gene.count.mat)/nrow(gene.count.mat)
bi <- gene.count.mat; bi[bi>0] <- 1

## ==================================================== ##
## -- Niche overlap score

s <- Sys.time()
statniche <- c()
#for(i in as.character(unique(mag$sample))){ #i=as.character(unique(mag $sample))[5]
i=as.character(unique(mag$sample))[day]    
    stmp <- proc.time()[3]
    ## --------------------------------- ##
    present <- colnames(ts)[which(ts[i,]>0.001)]
	sub=bi[present,]
	sub <- sub[,colSums(sub)>0]
    ## =========================== ## 
	
	sim <- mean(1-vegdist(sub, method='chao'))
	nodf <- nestednodf(sub)$statistic[3]
	h2 <- H2fun(sub, H2_integer = FALSE)[1]
	nodfc <- NODFc(sub, quality = 0)			
	 ## =========================== ## 
    statniche <- c(sim=sim, nodf, h2, nodfc)
    gc();   gc()
    
    etmp <- proc.time()[3]
    cat(sprintf('Day %s done: elapsed time %0.2f sec',i,  etmp - stmp))
#}

saveRDS(statniche, sprintf('%s/04_06_randamization_%s.rds', dir$rdsdir, i) )
e <- Sys.time()

print(e-s)
############################################################################

EOF


##=====================================================================================================##

cat <<EOF > ${qsubSript}

#!/bin/bash
#PBS -q ${que}
#PBS -N ${qname}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e ${log}
#PBS -j eo
#PBS -l walltime=12:00:00
## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load R

cd /user4/kyoto1/hiroakif/aptmp/

Rscript --vanilla ${RscriptName}
	

EOF

chmod +x ${qsubSript}
qsub ${qsubSript}

echo " "
echo "********* INFO: ${qname}.qsub was generated and started! *********"


##=====================================================================================================##

done