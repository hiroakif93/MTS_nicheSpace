############################################################################
####
#### R script for Fujita (2019)
####
#### Finding expected early warnings signal by visualizing
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_metagenome/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')

load.lib( c('ggplot2', 'tidyr', 'cowplot', 'ggh4x', 'car', 'RColorBrewer', 'extrafont', 'ggtext'))

# -- Create directory to save
dir <- make.dir('MAG_metagenome/Function_module')

# -- Load data table
path_comp <- read.table('Table/ShotgunMetagenome/keggdecorder_represent_bin.tsv', sep='\t', header=T, row.names=1, check.names = FALSE)
rownames(path_comp) <- gsub('.bin.', '_bin.', rownames(path_comp))
colnames(path_comp) <- paste(toupper(substring(colnames(path_comp), 1, 1)), substring(colnames(path_comp), 2), sep='')

mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds') ; dim(mag)
mag[mag==''] <- 'Unidentified'

ts <- readRDS('Table/04_MAG_timeSeries.rds')
ts <- ts/rowSums(ts)

cols <- readRDS('Table/04_color_palette.rds')
############################################################################
ini.upper <- sapply(strsplit(colnames(path_comp), split=''), function(x){ x[1] <- toupper(x[1]); return(x) })
names <- sapply(ini.upper, paste, collapse='')

names <- gsub('Nitrite', 'NO<sub>2</sub><sup>-</sup>', names)
names <- gsub('CO2', 'CO<sub>2</sub>', names)
names <- gsub('H2', 'H<sub>2</sub>', names)
names <- gsub('B12', 'B<sub>12</sub>', names)
names <- gsub('Nitric oxide', 'NO', names)
names <- gsub('Nitrous-oxide', 'N<sub>2</sub>O', names)
names <- gsub('phosphate', 'PO<sub>4</sub><sup>3-</sup>', names)
names <- gsub('Alpha', '<i>&alpha;</i>', names)

colnames(path_comp) <- names
## ======================================================================== ##
## -- MAG-level pathway
pathorder <- colSums(path_comp ==1)/2+colSums(path_comp)+colSums(path_comp>0)
binorder <- rowSums(path_comp ==1)+ rowSums(path_comp)+ rowSums(path_comp>0)

path_comp2 <- cbind(bin=rownames(path_comp), 
					Family=as.character(mag[ rownames(path_comp), 'Family']),
					Genus=as.character(mag[ rownames(path_comp), 'Genus']), 
					path_comp[,colSums(path_comp)>0] )
path_comp2$id <- path_comp2[,c(1)]#apply(path_comp2[,c(1)], 1, paste, collapse='_')

lf.path <- gather(path_comp2, key, value, -c(1,2,3, ncol(path_comp2)))
head(lf.path)

lf.path$key <- factor(lf.path$key, level=names(sort(pathorder, decreasing=T)))
lf.path$id <- factor(lf.path$id, level= unique(as.character(path_comp2[names(sort(binorder, decreasing=T)),'id']) ))

## ======================================================================== ##
## -- Present / Abesent matrix
tstmp <- ts
tstmp[tstmp>0.001] <- 1
tstmp[tstmp<=0.001] <- 0
df <- cbind(bin=colnames(ts), as.data.frame(t(tstmp)))
lf <- gather(as.data.frame(df), key, value, -1)
lf$key <- gsub('MTS-', '', lf$key)
lfbind <- cbind(lf, mag[as.character(lf[,'bin']), c(7:10)])

lfbind$key <- as.numeric(as.character(lfbind$key))
lfbind$key <- factor(lfbind$key, levels=rev(unique(lfbind$key)))

set <- unique(lfbind[,c(1,ncol(lfbind))])

order <- order(apply(mag[rownames(df),6:10], 1, paste, 1, collapse='_'))
levels <- c('Acidobacteriae', "Alphaproteobacteria","Gammaproteobacteria", "Bacilli")

lf.path$id <- factor(lf.path$id, levels=rownames(df)[order] )
lfbind$id <- factor(lfbind$bin, levels=rownames(df)[order] )

ggpath <- ggplot(lf.path)+
        geom_tile(aes(y=key, x=id, fill=value), color='grey90')+
        #geom_point(aes(x=id, y=-1))+
        scale_fill_gradientn(colors=brewer.pal(9, 'Greens'))+
        theme_bw()+
        scale_x_discrete(expand=c(0,0))+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_markdown(family='Arial', size=10,hjust=1),
              strip.background=element_blank())+
        labs(x='', y='')

shape=c(21:25, 4)
names(shape)=sort(unique(mag$Class))

ggts <- ggplot(lfbind)+
        geom_tile(aes(y=key, x=id, fill=ifelse(value==1, 'Present', 'Absent')), color='grey90')+
        geom_point(aes(x=id, y=0, shape=Class, fill=Genus), show.legend = FALSE, size=3)+
        scale_fill_manual(values=c(cols, Present='darkslategrey', Absent='white'))+
        scale_color_manual(values=cols)+
        theme_bw()+
        scale_x_discrete(expand=c(0,0))+
        scale_y_discrete(expand=c(0,1.8,0,0))+
        theme(axis.text.x=element_text(angle=40, hjust=1, vjust=1),
              axis.text.y=element_text(size=7),
              strip.background=element_blank())+
        scale_shape_manual(values=shape)+
        labs(x='', y='',fill='')

pdf(sprintf('%s/pathway_complete_v2.pdf',dir$figdir), width=13, height=21)
plot_grid(ggpath,NULL, ggts, ncol=1, rel_heights = c(0.85,-0.01, 0.1), align='v')
dev.off()

svglite::svglite(sprintf('%s/pathway_complete_v.svgz',dir$figdir), width=13, height=21)
plot(plot_grid(ggpath,NULL, ggts, ncol=1, rel_heights = c(0.85,-0.01, 0.1), align='v'))
dev.off()

## ======================================================================== ##
## -- ecosystem-level pathway

loaddir <- 'Table/ShotgunMetagenome/KO_number_represent'

filetmp <- list.files(loaddir)
file <- filetmp[grep('KD', filetmp)]

kd <- c()
for(i in 1:length(file)){
	
	df <- read.table( sprintf('%s/%s', loaddir, file[i]), sep='\t', header=T, 
				row.names=1, check.names = FALSE )
	
	mts <- gsub('.txt.tsv', '', strsplit(file[i], '_')[[1]][3]	)
	kd <- rbind(kd, cbind(sample=mts, df))
	
}

## ------------------------------ ##

colnames(kd) <- paste(toupper(substring(colnames(kd), 1, 1)), substring(colnames(kd), 2), sep='')

kdsub <- kd[,-1]
kdsub <- kdsub[,colSums(kdsub>0.9)<12]
kdsub <- kdsub[,colSums(kdsub)>0.9]

order=colSums(kdsub==1)+colSums(kdsub)+colSums(kdsub>0)

ecosys.path <- cbind(time=gsub('MTS-','Day ',rownames(kd)), kdsub )
ecosys.path[,1:10]
lf.path <- gather(ecosys.path, key, value, -c(1))
head(lf.path)

#lf.path$key <- gsub('OAA,', 'OAA,\n', lf.path$key)
lf.path$key <- factor(lf.path$key, level=names(sort(order, decreasing=T)))
lf.path$time2 <- as.character(as.numeric(gsub('Day ','', lf.path$time)))
lf.path$time2 <- factor(lf.path$time2, level=unique(lf.path$time2))

path.gg <- ggplot(lf.path)+
		   geom_tile(aes(x=key, y=time2 , fill=value), color='grey90')+
		   #geom_hline(yintercept=11.5, color='red')+
		   scale_fill_gradientn(breaks=seq(0,1,0.2), colors=brewer.pal(8,'Greens'))+
		   theme_md(bsize=6)+
		   theme(text=element_text(family='Arial'),
		   		 axis.text.x=element_text(size=12, angle=40, hjust=0.98, vjust=1),
		   		 axis.text.y=element_text(size=10))+
		   coord_flip()+
		   scale_x_discrete(expand=c(0,0))+
		   scale_y_discrete(expand=c(0,0))+
		   labs(x='', y='')		   
	
ggsave(plot= path.gg, h=5, w=11,
       filename=sprintf('%s/ecosystem_pathway_complete.pdf',dir$figdir))		

## ======================================================================== ##
path_comp3 <- path_comp2[rownames(path_comp2)[c(26,27,16,29,25)],]
#path_compdomi <- path_compdomi[ , c(1:3,ncol(path_compdomi), which( colSums(path_compdomi[,-c(1:3,ncol(path_compdomi))]>0.8)/nrow(path_compdomi) < 1 & 
#                                                                        colSums(path_compdomi[,-c(1:3,ncol(path_compdomi))]>0.8)/nrow(path_compdomi) > 0.1)+3)]

path_compdomi_magCol <- path_comp3[,c(1:3,ncol(path_comp3))]
path_compdomi_pathCol <- path_comp3[,-c(1:3,ncol(path_comp3))]
path_compdomi_less5 <- path_compdomi_pathCol[,which((colSums(path_compdomi_pathCol>0.1)/nrow(path_compdomi_pathCol)) < 1)]
path_compdomi_less5 <- path_compdomi_less5[,colSums(path_compdomi_less5>0.9)>0]

path_compdomi <- cbind(path_compdomi_magCol, path_compdomi_less5)
path_compdomi$id <- apply(path_compdomi[,c('Genus', 'id')], 1, paste, collapse=' ')

lf.path <- gather(path_compdomi, key, value, -c(1:4))
head(lf.path)

rowSums(path_compdomi_less5)
lf.path$key <- factor(lf.path$key, level=names(sort(colSums(path_compdomi_less5), decreasing=T)))
lf.path$id <- factor(lf.path$id, level= path_compdomi$id[order(rowSums(path_compdomi[,-c(1:(ncol(path_compdomi_magCol)+1))]))] )
head(lf.path)
path.gg <- ggplot(lf.path)+
            geom_tile(aes(y=key, x= id, fill=value), color='grey60')+
            scale_fill_gradientn(breaks=seq(0,1,0.2), colors=brewer.pal(8,'Greens'))+
            theme_bw(base_size=6)+
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3),
                  axis.text.y=element_markdown(family='Arial', size=5,hjust=1),
                  text=element_text(family='Arial'))+	
            scale_y_discrete(expand=c(0,0))+
            scale_x_discrete(expand=c(0,0))+
            labs(x='', y='')

pdf(sprintf('%s/pathway_dominant.pdf',dir$figdir), width=10/2.5, height=13/2.5)
plot(path.gg)
dev.off()
svglite::svglite(sprintf('%s/pathway_dominant.svgz',dir$figdir), width=8/2.5, height=13/2.5)
plot(path.gg)
dev.off()