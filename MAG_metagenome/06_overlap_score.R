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
library(AnalysisHelper)
load.lib( c('ggplot2', 'RColorBrewer', 'tidyr', 'vegan', "ggrepel", 'extrafont','bipartite','maxnodf', 'ggtext', 'FactoClass', 'scatterplot3d'))

# -- make directory to save results
dir <- make.dir('MAG_metagenome/Niche_overlap')

# -- Load data table
gene.count.mat <- readRDS('Table/04_03_gene_count.matrix_represent.rds')
mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds') 
mag[mag ==''] <- 'Unidentified'

repmag <- mag[unique(mag$represent),]
write.csv(repmag, 'Table/ShotgunMetagenome/repMAG_infomation.csv', row.names=TRUE)

ts <- readRDS('Table/04_MAG_timeSeries.rds')
write.csv(ts, 'Table/ShotgunMetagenome/MAG_timeSeries.csv', row.names=TRUE)

ts <- ts/rowSums(ts)
cols <- readRDS('Table/04_color_palette.rds')

for.parallel(8)
############################################################################

shape=c(21:25, 4)
#names(shape) <- c("Gammaproteobacteria", "Alphaproteobacteria", "Bacilli", 
#                  "Bacteroidia", "Acidobacteriae", "Negativicutes"  )


freq <- colSums(gene.count.mat)/nrow(gene.count.mat)
bi <- gene.count.mat#[, freq < 0.95]; 
bi[bi>0] <- 1

## ---------------------------------------- ##
## -- PCoA and clustering
if(F)
{
	dist <- vegdist(bi, 'jaccard')
	pca <- cmdscale(dist, k=3)
	saveRDS(pca, 'Table/04_06_coord.rds')
}

pca <- readRDS('Table/04_06_coord.rds')

pcadf <- cbind(repmag, PC= pca[rownames(repmag),])
shape=c(21:25, 4)
names(shape) <- sort(unique(pcadf$Class) )

g1 <- ggplot(pcadf)+
      geom_point(aes(x=PC.1, y= PC.2, fill= Genus, shape= Class), size=3,
                 show.legend = FALSE)+
      theme_bw()+
      scale_fill_manual(values=cols[(pcadf$Genus)] ) +
      scale_color_manual(values=cols[(pcadf$Genus)] ) +
      theme_bw(base_size=8)+
      theme(text=element_text(family='Arial'),
            strip.background = element_blank(),
            strip.text=element_text(size=16) )+
      scale_shape_manual(values=c(21:25, 4))+
      labs(x='PCoA 1', y='PCoA 2')
      
ggsave(filename=sprintf('%s/Niche_spapce.tiff', dir$figdir),
       plot= g1, w=9, h=5)
svglite::svglite(filename = sprintf('%s/Niche_spapce.svgz', dir$figdir), 
                 width=6/2.5, height=5/2.5)
plot(g1); dev.off()    

g1 <- ggplot(pcadf)+
    geom_point(aes(x=PC.1, y= PC.2,  fill= Genus), size=3,
               show.legend = T,shape=21)+
    guides(fill=guide_legend(nrow=6))+
    scale_fill_manual(values=cols[(pcadf$Genus)] ) +
    theme_bw()+
    theme_bw(base_size=8)+
    theme(text=element_text(family='Arial'),
          strip.background = element_blank(),
          strip.text=element_text(size=16) )+
    labs(x='PCoA 1', y='PCoA 2')

svglite::svglite(filename = sprintf('%s/Niche_spapce_color.svgz', dir$figdir), 
                 width=14/2.5, height=5/2.5)
plot(g1); dev.off()   
## ---------------------------------------- ##
## -- Contribution to PCoA

cormat <- cor(pca, bi, method='spearman')
a <- apply(abs(cormat), 1,function(x) order(x, decreasing=TRUE)[1:5] )
pca_cont <- apply(a, 2, function(x){colnames(cormat)[x]})
sink(sprintf('%s/PCA_contribution_top5.txt', dir$tabledir)); print(pca_cont); sink()
if(F){
	write.csv(cbind(c(rep('PCoA1',5),rep('PCoA2',5),rep('PCoA3',5)), c(pca_cont)),
		  sprintf('%s/PCA_contribution_top5.csv', dir$tabledir), row.names=F, col.names=F)
}

lf.merge <- c()
for(i in 1:3){
	
	tmp <- cbind(mag[rownames(pca), c('Genus', 'Class')], score=pca[,i], bi[,pca_cont[,i]])
	
	lf <- gather(as.data.frame(tmp), key, value, -c(1:3))
	lf.merge <- rbind(lf.merge, cbind(axis=i, lf))
	
	g <- ggplot(lf, aes(x=score, y=value))+
		 geom_point(aes(shape= Class, fill= Genus))+
		 facet_wrap(~key)+
      	 theme_bw(base_size=10)+
      	 scale_fill_manual(values=cols[(pcadf$Genus)] ) +
      	 scale_color_manual(values=cols[(pcadf$Genus)] ) +
      	 scale_shape_manual(values=shape)+
      	 theme(strip.background = element_blank(),
            	strip.text=element_text(size=8) )
    ggsave(filename=sprintf('%s/Niche_spapce_PCoA%s.pdf', dir$figdir, i),
       plot= g, w=9, h=5)
              	
}

## ---------------------------------------- ##
## -- PCoA
pcoaTS <- smetanaComm <- c()
for(i in rownames(ts)){ #i=rownames(ts)[9]
	
	present <- colnames(ts)[which(ts[i,]>0.001)]
	time=paste('Day', as.numeric(strsplit(i, '-')[[1]][[2]]))
	
	pcoaTS <- rbind(pcoaTS, cbind(time=time,pcadf[present, ]))
	
	com <- cbind(time=as.numeric(strsplit(i, '-')[[1]][[2]]), paste(mag[present,'represent'], 'rename', sep='_'))
	smetanaComm <- rbind(smetanaComm, com)
}

write.table(smetanaComm, 'Table/smetanaCommunity.tsv', sep='\t', row.names=FALSE, col.names=FALSE, quote=F)
unique(pcoaTS$Class)


pcoaTS[pcoaTS$time=='Day 70',]
pcoaTS[grep('Terra', pcoaTS$Genus),]
g2 <- ggplot(pcoaTS)+
      geom_point(aes(x=PC.1, y=PC.2, fill= Genus, shape=Class), size=5, show.legend=T)+
      theme_bw()+
      scale_fill_manual(values=cols[(pcoaTS$Genus)] ) +
      theme_bw(base_size=13)+
      theme(text=element_text(family='Arial'),
            strip.background = element_blank(),
            strip.text=element_text(size=13) )+
      scale_shape_manual(values=shape)+
      labs(x='PC 1', y='PC 2')+
      facet_wrap(~time, nrow=2)+
      xlim(c(min(pcoaTS$PC.1)-0.1, max(pcoaTS$PC.1)+0.1))+
      ylim(c(min(pcoaTS$PC.2)-0.1, max(pcoaTS$PC.1)+0.1))
      
ggsave(filename=sprintf('%s/Niche_spapce_day.tiff', dir$figdir),
       plot= g2, w=11, h=4)    

g4 <- ggplot(pcadf)+
	  geom_point(aes(x=genome_size/10e+5, y=GC.content, fill= Genus, shape=Class),
	             size=3, show.legend=FALSE)+
	  #geom_text(aes(x=genome_size/10e+5, y=GC.content, label= Genus),
	  #           size=2, show.legend=FALSE)+
	  scale_fill_manual(values=cols[as.character(pcadf$Genus)] ) +
	  theme_bw(base_size=8)+
	  theme(text=element_text(family='Arial'),
	        strip.background = element_blank(),
        	        legend.background = element_blank()   )+
	  scale_shape_manual(values=shape)+
	  labs(y='GC content(%)', x='Genome size(mb)')
ggsave(filename=sprintf('%s/gc_vs_genomesize.tiff', dir$figdir),
       plot=g4,w=10, h=7)    
svglite::svglite(filename = sprintf('%s/gc_vs_genomesize.svgz', dir$figdir), width=6/2.5, height=5/2.5)
plot(g4); dev.off()



col=cols[(pcadf$Genus)]
svglite::svglite(sprintf('%s/Niche_spapce.svgz', dir$figdir), width=9/2.5, height=9/2.5)
s3d<-scatterplot3d::scatterplot3d(x= pca[,1], z=pca[,2], y=pca[,3],
		  pch="", grid=FALSE, box=FALSE, xlab='PCoA 1', ylab='PCoA 2', zlab='PCoA 3', cex.lab=0.5, cex.axis=0.5)
addgrids3d(x= pca[,1], z=pca[,2], y=pca[,3], grid = c("xy", "xz", "yz"))
s3d$points3d(x= pca[,1], z=pca[,2], y=pca[,3],type = "h",
		  pch=shape[pcadf$Class], bg=col, cex=1.5)       
dev.off()	

pdf(sprintf('%s/Niche_spapce_day.pdf', dir$figdir), width=18, height=11)
par(mfrow=c(3,5))
for(i in rownames(ts)){ #i=rownames(ts)[9]
	
	present <- colnames(ts)[which(ts[i,]>0.001)]
	time=paste('Day', as.numeric(strsplit(i, '-')[[1]][[2]]))
	
	pcaday <- pcadf[present, ]
	
	
	s3d<-scatterplot3d::scatterplot3d(x= pcadf[,'PC.1'], z= pcadf[,'PC.2'], y= pcadf[,'PC.3'],
		  pch="", grid=FALSE, box=FALSE, xlab='PCoA 1', ylab='PCoA 2', zlab='PCoA 3', 
		  cex.lab=0.5, cex.axis=0.5, 
		  main= time)
	addgrids3d(x= pcadf[,'PC.1'], z= pcadf[,'PC.2'], y= pcadf[,'PC.3'], grid = c("xy", "xz", "yz"))
	s3d$points3d(x= pcaday[,'PC.1'], z= pcaday[,'PC.2'], y= pcaday[,'PC.3'],type = "h",
		  pch=shape[pcaday $Class], bg=col[pcaday $Genus], cex=1.5)   
}
dev.off()

## ==================================================== ##
## -- Niche overlap score
#if(F){
		
	simlist <- statniche <- c()	
	for(i in as.character(unique(mag$sample))){ #i=as.character(unique(mag $sample))[5]
		    stmp <- proc.time()[3]
		    ## --------------------------------- ##
		    present <- colnames(ts)[which(ts[i,]>0.001)]
			sub=bi[present,]
			sub= sub[,colSums(sub)>0]
			
		    ## =========================== ## 
			simlist[[i]] <- as.vector(1-vegdist(sub, method='jaccard', binary=TRUE))

			sim <- mean(1-vegdist(sub, method='jaccard', binary=TRUE))
			#nodf <- nestednodf(sub)$statistic[3]
			#h2 <- H2fun(sub, H2_integer = FALSE)[1]
			#nodfc <- NODFc(sub, quality=2)			
			 ## =========================== ## 
		    statniche <- c(statniche, sim=sim)
		    gc();   gc()
		    
		    etmp <- proc.time()[3]
		    cat(sprintf('Day %s done: elapsed time %0.2f sec',i,  etmp - stmp))
	}
#}

## ==================================================== ##
fl <- list.files(sprintf('%s', dir$rdsdir), full.names=T)

df <- (lapply(fl[grep('_MTS-', fl)], readRDS))
## ==================================================== ##

tsbi <- ts
tsbi[tsbi>0] <- 1
jac <- as.matrix(vegdist(tsbi, method='jaccard', binary=TRUE))
bray <- as.matrix(vegdist(ts, method='bray'))

abrupt <- function(x){
	abruptness <- c()
	for(i in 2:ncol(x)){
	    abruptness <- c(abruptness, x[i-1, i])
	}
	abruptness <- c(abruptness, NA)
	return(abruptness)
}

abrupts <- cbind(jaccard=abrupt(jac), 'Community change'=abrupt(bray))
lap_score <- data.frame( time=c(1,10,20,24,seq(30,110,10)), Similarity=statniche, abrupts, check.names = FALSE)

#nichelap <- do.call(rbind, df)
#colnames(nichelap) <- c('Similarity', 'NODF', 'H2', 'NODFc')
#lap_score <- data.frame( time=c(1,10,20,24,seq(30,110,10)), nichelap, abrupts, check.names = FALSE)

# nicheOverlap <- c()#lapply(df, data.frame, time=c(1,10,20,24,seq(30,110,10)),richness=rowSums(ts>0), 'Abruptness'=abruptness)
# for(i in 1:10){
# 	x= df[[i]]; names(x)[4]=c('nodfc')
# 	nicheOverlap <- rbind(nicheOverlap, data.frame(th=i/10, time=c(1,10,20,24,seq(30,110,10)),richness=rowSums(ts>0), abrupts,x))
# }

lap_score_Scale <- lap_score
lap_score_Scale[,-1] <- apply(lap_score_Scale[,-1], 2, scales::rescale, to=c(0,1))
lf <- gather(lap_score_Scale, key, value, -c(1,3))

g1 <- ggplot(lf, aes(x=time, y=value))+
        geom_line(aes(color=key, group=key))+
        scale_color_manual(values=c(NODFc='slateblue3', Similarity='firebrick3',
                                    'Community change'='black'))+
        labs(y='Rescaled scores', x='Time', color='Score')+
        theme_text(bsize=6, family = 'Arial')+
        theme(legend.key = element_rect(fill=NA))

g2 <- ggplot(lap_score, aes(y=`Community change`, x= Similarity, color=time))+
        geom_point()+
        labs(y='Community change', x='Jaccard similarity')+
      scale_color_gradientn(colors=brewer.pal(11, 'Spectral'))+
      theme_text(family = 'Arial')+
      theme(legend.key = element_rect(fill=NA),
            strip.placement = 'outside')

ggsave(plot=g1, sprintf('%s/score_ts.pdf', dir$figdir),device = cairo_pdf, w=9, h=4, unit='cm')
ggsave(plot=g2, sprintf('%s/corerelation_with_abruptness.pdf', dir$figdir),device = cairo_pdf, w=14, h=7, unit='cm')


# fitting <- function(x){
#     
#     fml <- formula(sprintf("`Community change`~%s", x))
#     
#     glmmres <- c()
#         
#         sub <- lap_score
#         model <- lm(fml, data=sub)	
#         coef <- summary(model)$coefficients
#         tmp <- data.frame(key =i, pval=summary(model)$coefficients[2,4], 
#         				  tval=sprintf('T value=%.2f', summary(model)$coefficients[2,3]),  
#         				  a=coef[1,1], b=coef[2,1], 
#                           fml=sprintf('y = %.2fx + %.2f',coef[2,1], coef[1,1]),
#                           x=max(sub[,x], na.rm=TRUE),  y=max(sub$Abruptness, na.rm=TRUE)+0.05)
#         glmmres <- rbind(glmmres, tmp)
#     
#     glmmres $pval_lab <- sprintf('<i>P</i> = %.2fx10<sup>-2</sup>', glmmres $pval*100)
# 	glmmres $pval_lab[glmmres $pval<0.001] <- '<i>P</i> < 1x10<sup>-3</sup>'
#     glmmres[glmmres $pval>0.001, c('a', 'b')] <- NA
#     
#     return(glmmres)
# }

# nicheOverlap <- lapply(df, data.frame, time=c(1,10,20,24,seq(30,110,10)),richness=rowSums(ts>0), abrupts)
# jac_nodfc <- sapply(nicheOverlap, function(x){ cor(x[-nrow(x),2], x[-nrow(x),'jaccard'], method='spearman') })
# jac_nodf <- sapply(nicheOverlap, function(x){ cor(x[-nrow(x),1], x[-nrow(x),'jaccard'], method='spearman') })
# bray_nodfc <- sapply(nicheOverlap, function(x){ cor(x[-nrow(x),2], x[-nrow(x),'bray'], method='spearman') })
# bray_nodf <- sapply(nicheOverlap, function(x){ cor(x[-nrow(x),1], x[-nrow(x),'bray'], method='spearman') })
# 
# pdf(sprintf('%s/pairs.pdf', dir$figdir))
# pairs(cbind(th=seq(0.1, 1, 0.1), jac_nodf, jac_nodfc, bray_nodfc, bray_nodf))
# dev.off()
# 
# pdf(sprintf('%s/nodfc_bray.pdf', dir$figdir))
# par(mfrow=c(3,4))
# for(i in 1:10){
# 	x= nicheOverlap[[i]]
# 	plot(x[-nrow(x),2], x[-nrow(x),'bray'],
# 	     main=sprintf('th%s',i/10), xlab='nodfC', ylab='abruptness')
# }
# dev.off()
# 
# pdf(sprintf('%s/nodfc_jac.pdf', dir$figdir))
# par(mfrow=c(3,4))
# for(i in 1:10){
#     x= nicheOverlap[[i]]
#     plot(x[-nrow(x),2], x[-nrow(x),'jaccard'],
#          main=sprintf('th%s',i/10), xlab='nodfC', ylab='abruptness')
# }
# dev.off()
# 
# lmres <- c()
# for(i in c('Similarity','NODFc')){
# 	lmres <- rbind(lmres, fitting(x=i))
# }
# 
# 
# lf <- gather(as.data.frame(lap_score), key, Scores, c(2,5))
# g <- ggplot(lf)+
# 	 geom_point(size=1.2, aes(x= Scores, y= `Community change`, color=time))+
# 	 geom_abline(data= lmres, aes(intercept=a, slope=b), color='royalblue', linetype=2)+
#      geom_text(data= lmres,  aes(x=x, y=0,label= pval), size=2, hjust=1, fill = NA, label.color = NA)+
# 	 facet_wrap(~key, scales='free',ncol=2,dir='v',
# 	 		strip.position='bottom')+
#      scale_color_gradientn(colors=brewer.pal(11,'Spectral'))+
# 	 theme_bw(base_size=8)+
# 	 theme(text=element_text(family='Arial'),
# 	 		strip.background=element_blank(),
# 	 		strip.placement='outside')+
# 	 labs(x='', y='Observed community change ')		
# 
# svglite::svglite(filename = sprintf('%s/niche_overlap_indices.svgz', dir$figdir), width=8/2.5, height=7/2.5)
# plot(g); dev.off()
# ggsave(plot=g, filename=sprintf('%s/niche_overlap_indices.pdf', dir$figdir),device = cairo_pdf, height=8, width=14, unit='cm')

## ==================================================== ##
genemat=bi
freq <- apply(genemat, 2, function(x)sum(x>0))/nrow(genemat)
time <- rownames(ts)

agg <- c()
for(i in unique(time)){
    
    present <- colnames(ts)[which(ts[i,]>0.001)]
    agg <- rbind(agg ,colSums(genemat[present, ]))
    
}

agg[agg>0] <- 1
pca <- cmdscale(vegdist(agg, 'jaccard'))

pca2 <- data.frame(time=time, lap_score, pca[,1:2])
pca2$time <- as.numeric(gsub('MTS-', '', pca2$time))
ggpcasim <- ggplot(pca2)+
         geom_path(aes(x=X1, y= X2), show.legend=FALSE,size=0.4,
                  arrow = arrow(type='closed',length = unit(0.015, "npc")), linetype=3)+
         geom_point(aes(x=X1, y= X2, fill= Similarity), shape=21, size=3)+         
         geom_text_repel(data=pca2, aes(x= X1, y= X2, label= sprintf('Day %s', time)), size=2)+		 
         viridis::scale_fill_viridis(guide = guide_colourbar(direction = "horizontal"))+
         scale_color_gradientn(colors=c( rev(brewer.pal(9,'Greys'))[1:2] ) )+
    theme_text()+
         theme(text=element_text(family='Arial'),
               legend.text=element_text(size=6),
               legend.title=element_text(size=6, vjust=0.8),
               legend.key.width=unit(0.4,"cm"),
               legend.key.height=unit(0.2,"cm"),
               legend.position=c(1, 0),
               legend.justification = c(1,2),
               plot.margin= unit(c(1, 1, 9, 1), "lines"))+
         labs(x='Axis 1',y='Axis 2',fill="Niche overlap score (Similarity)")+
    coord_fixed( max(pca2$X1)/max(pca2$X2))



ggsave(plot=ggpcasim, filename=sprintf('%s/com_assembly_similarity.pdf', dir$figdir),device = cairo_pdf, height=12, width=10, unit='cm')
