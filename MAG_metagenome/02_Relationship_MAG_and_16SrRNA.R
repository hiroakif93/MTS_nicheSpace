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

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('ggplot2', 'RColorBrewer', 'ggtext', 'tidyr', 'extrafont', 'cowplot'))

## -- make directory to save results
dir <- make.dir('MAG_metagenome/MAG_and_16SrRNA_relation')

# -- Load data table
mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds') ; dim(mag)
mag <- mag[unique(mag$represent), ]

mag_ts <- read.table('Table/ShotgunMetagenome/MAG_abundacne.txt', sep='\t', header=T, stringsAsFactors=F); dim(mag_ts)
split <- do.call(rbind, strsplit(mag_ts[,2], '\\.'))
mag_ts[,2] <- apply(split[,1:2], 1, paste, collapse='.')

dlist <- readRDS('Table/matrixList.rds')

sml <- readRDS('Table/sample_info.rds') 
sml <- sml[which(sml $treat1=="Water/medium A" & sml $replicate.id==5), ]
sml <- sml[as.numeric(sml$time)%in%c(1,10,20,24,30,40,50,60,70,80,90,100,110), ]

amp_ts <- dlist[["Water/medium A"]]
amp_ts <- ts[rownames(sml), colSums(ts[rownames(sml), ])>0]

amp_taxa <- as.matrix(readRDS('Table/04_Taxa_list.rds'))
amp_taxa <- amp_taxa[,grep('gbtk', colnames(amp_taxa))]

## -- Parameter
check_taxa <- 'Genus'
abundance='RPKM'

#non.prok <- c('Other than prokaryotes', 'Unidentified')
non.prok <- c('Unidentified')
############################################################################
## -- Making sample-bin matrix

mat <- matrix(0, ncol=nrow(mag), nrow=length( unique(mag_ts[,1])),
			  dimnames=list(unique(mag_ts[,1]), rownames(mag)) )

for(i in rownames(mat)){#i='MTS-001'
	
	tssub <- mag_ts[mag_ts[,1]==i,]
	tssub <- tssub[tssub[,2] %in% colnames(mat), ]
	
	mat[i, tssub[,2] ] <- tssub[,abundance]
	
}

saveRDS(mat, 'Table/04_MAG_timeSeries.rds')

## =========================================== ##
## -- Merging taxonomy information
tstmp <- na.omit(mag_ts)

taxarep <- mag[tstmp[,2], 5:11]
#taxarep[is.na(taxarep)] <- 'Other than prokaryotes'
taxarep[taxarep==''] <- 'Unidentified'

ts.info <- cbind(time=as.numeric(gsub('MTS-', "", tstmp[,1])),
				 tstmp, taxarep)

## =========================================== ##
## -- Colors

sep <- 5
col1 <- c(  'firebrick4', 'firebrick3', 'indianred1' ) ; plot(c(1:length(col1))~1, col=col1, pch=19)
col2 <- c( 'royalblue3', 'dodgerblue2','skyblue2') ; plot(c(1:length(col2))~1, col=col2, pch=19)
col3 <- c( 'olivedrab4','springgreen4', 'green3', 'darkolivegreen3', 'darkolivegreen1') ; plot(c(1:length(col3))~1, col=col3, pch=19) #darkolivegreen4
col4 <- c( 'mediumorchid4', 'mediumpurple') ; plot(c(1:length(col4))~1, col=col4, pch=19) #'mediumorchid','plum2',  
col5 <- c('goldenrod4', 'darkgoldenrod1') ; plot(c(1:length(col5))~1, col=col5, pch=19)
col6 <- rev(c('cornsilk1','cornsilk3')); plot(c(1:length(col6))~1, col=col6, pch=19)
col7 <- rev(c('sienna3')); plot(c(1:length(col7))~1, col=col7, pch=19)

Acidobacteriae <- 'plum2'#'navy'
Alphaproteobacteria <- c(col2, col4); plot(c(1:length(Alphaproteobacteria))~1, col=Alphaproteobacteria, pch=19)
Bacilli <- col6
Bacteroidia <- col5
Gammaproteobacteria <- c(col1, col7, col3); plot(c(1:length(Gammaproteobacteria))~1, col=Gammaproteobacteria, pch=19)

tes <- unique(ts.info[,c('Phylum','Class','Genus')])
tes[order(tes$Class),]
col <- data.frame(row.names=unique(tes[order(tes$Class), 'Genus']))

col[-which(rownames(col) %in% c(non.prok, 'Hydrotalea')), 'color'] <- c(Acidobacteriae, Alphaproteobacteria, Bacilli, Bacteroidia, Gammaproteobacteria)
    
col[rownames(col)=='Unidentified','color'] <- 'grey30'
col[rownames(col)=='Other than prokaryotes','color'] <- 'grey70'
col[rownames(col)=='Hydrotalea','color'] <- 'yellow'
plot(1:21, pch=19, col=col[,1]) #c(col1, col2, col3, col4, col5, col6))

cols <- col[,1]; names(cols) <- rownames(col)
names(cols)[c(17, 19)] <- names(cols)[c(19, 17)]

cols <- c(cols[13:21], cols[10:12], cols[2:7], cols[1], cols[c(8,9)])

saveRDS(cols, 'Table/04_color_palette.rds')

## =========================================== ##
## -- Visualize

names(cols)[ -which(names(cols)%in% non.prok) ] <- paste('*', names(cols)[ -which(names(cols)%in% non.prok) ], '*', sep='')

sum <- aggregate(ts.info[, abundance], by=list(taxa=ts.info[,'Genus'], time=ts.info[,'time']), sum)
sum[-which(sum[,'taxa']%in% non.prok),'taxa'] <- paste('*', sum[-which(sum[,'taxa']%in% non.prok),'taxa'], '*', sep='')

order <- aggregate(sum[,3], by=list(taxa= sum[,1]), sum)
order <- order[ -which(order[,1] %in% non.prok), ]

sum$taxa <- factor(sum$taxa, levels= c(non.prok , order[order(order[,2]),1]))

level <- as.character(sum[order(sum[,2]), 1])

gg <- ggplot(sum) +
	  geom_bar(aes(x= as.factor(time), y=x, fill = as.factor(taxa)), 
	  		   width=1,stat='identity', position='fill', size=0.4, color='black', show.legend=FALSE)+
	  theme_bw(base_size=8)+
	  labs(x='Day', y='Relative abundance based on<br>reads per kilobase of genome<br>per million reads mapped', fill="Genus")+
	  theme(text=element_text(family='Arial'), 
          legend.text= element_markdown(family='Arial'),
          axis.title.y=element_markdown(family='Arial'))+
	  scale_fill_manual(values= cols[c(names(cols)[-which(names(cols)%in%'Unidentified')], non.prok)] )+
	  scale_y_continuous(expand=c(0,0))+scale_x_discrete(expand=c(0,0))

ggsave(plot=gg, w=8/2.5,h=5/2.5,
		sprintf('%s/MAG_assembly.pdf', dir$figdir))	
svglite::svglite(sprintf('%s/MAG_assembly.svgz',dir$figdir), width=8/2.5, height=5/2.5)
plot(gg); dev.off()		

gg <- ggplot(sum) +
	  geom_point(aes(x= 2, y=2, fill = as.factor(taxa)), shape=21, 
	  		   size=3, color='black', show.legend=T)+
	  theme_bw(base_size=8)+
	  theme(text=element_text(family='Arial'), 
          legend.text= element_markdown(family='Arial'),
          axis.title.y=element_markdown(family='Arial'),
        	        legend.background = element_blank() ,
          			legend.spacing.x =unit(0, 'cm'))+
	  scale_fill_manual(values= cols[c(names(cols)[-which(names(cols)%in%'Unidentified')], non.prok)] )+
	  guides(fill=guide_legend(nrow=6))
svglite::svglite(sprintf('%s/MAG_assembly_legend.svgz',dir$figdir), width=15/2.5, height=5/2.5)
plot(gg); dev.off()		  
	  
## =========================================== ##
## -- Relationship between MAG relative abundances 
##    and ASV relative abundance

## -- Sum up ASV abunadnce in Family level abundance
ampsub <- amp_ts[rownames(sml),]
ampGenus <- t(Taxa.mat(ampsub, amp_taxa, 'gbtk.Family'))

## -- Sum up MAG abunadnce in Family level abundance
rel <- mat
binF <- t(Taxa.mat(na.omit(rel), mag[, 5:11], 'Family'))
binF <- binF/rowSums(binF)

## =========================================== ##

## -- ASV and MAG abundance table
df <- unique(rbind( expand.grid(1:nrow(binF), colnames(binF)), expand.grid(1:nrow(binF), colnames(ampGenus)) ))

df$amplicon <- 0
for(i in 1:13){
	
	tmp <- ampGenus[i, ]/sum(ampGenus[i, ])
	tmp2 <- df[which(df[,1]==i & df[,2] %in% names(tmp)),]
	
	df[which(df[,1]==i & df[,2] %in% names(tmp)),3] <- tmp[as.character(tmp2[,2])]
}

df$bin <- 0
for(i in 1:13){
	
	tmp <- binF[i, ]
	tmp2 <- df[which(df[,1]==i & df[,2] %in% names(tmp)),]
	
	df[which(df[,1]==i & df[,2] %in% names(tmp)),4] <- tmp[as.character(tmp2[,2])]
}

## =========================================== ##
## -- Visualize 

total <- aggregate(df[,3:4], by=list(df$Var2), sum, na.rm=TRUE)
total$Group.1 <- gsub('unidentified', 'Unidentified', total$Group.1)

total$sum <- rowSums(total[,2:3])
total <- total[order(total$sum, decreasing=TRUE), ]
total$color <- 'Others'
total[1:11,'color'] <- as.character(total[1:11, 1])

col <- c(brewer.pal(6, 'Set1'), brewer.pal(5, 'Set2')[-2], 'grey80', 'Grey30')
names(col) <- c(setdiff(as.character(total[1:11,'color'] ), 'Unidentified'),  'Unidentified', 'Others')

a <- cor.test(na.omit(df)$amplicon, na.omit(df)$bin, method='spearman')

g <-ggplot(df[rowSums(df[,3:4])>0,], aes(x=amplicon, y=bin, color= Var2))+
	geom_abline(a=1, linetype=2, color='grey40', size=1.5)+
	geom_point(size=5)+	
    theme_bw(base_size=25)+
    theme(text=element_text(family='Arial'), 
    	  plot.margin=unit(c(1,1.5,1,1), 'lines'))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    scale_color_manual(values=col)
    
ggsave(filename=sprintf('%s/family.tiff', dir$figdir),
	   plot=g, w=10,h=6)    

## =========================================== ##

## -- For KEGG-Decoder and SMETANA
dir.create('Table/KO_number_represent')
ko <- read.csv('Table/ShotgunMetagenome/ghost_koala_represent_rename.csv', header=TRUE)
ko[,1] <- gsub('.bin', '_bin', ko[,1])

mat <- mat/rowSums(mat)
smetana.comm <- c()
for(i in 1:nrow(mat)){ #i=1
	
	bin <- colnames(mat)[which(mat[i,]>0.001)]
	
	smetana.comm <- rbind(smetana.comm, cbind(community=rownames(mat)[i], bin= paste(bin, '_rename', sep='')))
	
	rows <- c()
	for(l in bin){
		rows <- c(rows, grep( gsub('_bin', '.bin', l), ko[,1]))
	}
	
	kotmp <- ko[rows,]
	kotmp[,1] <- rownames(mat)[i]

	write.table(kotmp, sprintf('Table/ShotgunMetagenome/KO_number_represent/KO_%s.txt', rownames(mat)[i]),
				sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}

write.table(smetana.comm, 'Table/ShotgunMetagenome/smetana_commmunity.txt',
				sep='\t', quote=FALSE, row.names=FALSE)
				
				
############################################################################
