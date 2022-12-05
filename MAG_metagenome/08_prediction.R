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
load.lib( c('ggplot2', 'RColorBrewer', 'tidyr', 'vegan', "ggrepel", 'extrafont'))

# -- make directory to save results
dir <- make.dir('04_Dynamics_in_metagenomic_niche_space/Prediction_abruptChange')

# -- Load data table
rand <- readRDS('Table/04_06_randamization.rds')
smet <- read.table('Table/ShotgunMetagenome/CarveMe2Smetana/smetana/mip_mro_global.tsv', header=TRUE)
smet$time <- formatC(smet[,1], width=3, flag='0')
smet <- smet[order(smet$time), ]

ts <- readRDS('Table/04_MAG_timeSeries.rds')
ts <- ts/rowSums(ts)

############################################################################

distmat <- as.matrix(vegdist(ts))
abrupt <- c()
for(i in 2:nrow(ts)){ abrupt <- c(abrupt, distmat[i, i-1])}
abrupt <- c(abrupt, NA)

df <- cbind(abrupt, niche=rand[,1], smet[,c('mip', 'mro')])
cor( df, use='pairwise.complete.obs', method='spearman')
pairs(cbind(abrupt, -rand[,1], smet[,c('mip', 'mro')]))

exp <- c('mip', 'mro', 'mip+mro')

for(i in exp){
	
	fml <- formula(sprintf('abrupt~%s', i))
	model <- lm(fml, data=df)
	
	print(AIC(model))
}

modelfull <- lm(abrupt~mip+mro, data=df)
step(modelfull)

