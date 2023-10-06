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
load.lib( c('ggplot2', 'RColorBrewer', 'tidyr', 'vegan', "ggrepel", 'extrafont','bipartite','maxnodf', 'ggtext', 'FactoClass', 'scatterplot3d', 'gridExtra'))

# -- make directory to save results
dir <- make.dir('MAG_metagenome/Niche_overlap')

# -- Load data table
sml <- readRDS('Table/03_04_sample_info_2.rds')
smltarget <- sml[sml$treat1=='Water/medium A' & 
				 sml$replicate.id==5 &
				 sml$time%in%paste('Day_', formatC(c( 1, 10, 20, 24, seq(30, 110, 10)), width=3, flag='0'), sep=''),]

gene.count.mat <- readRDS('Table/04_03_gene_count.matrix_represent.rds')
mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds') 
mag[mag ==''] <- 'Unidentified'

repmag <- mag[unique(mag$represent),]

ts <- readRDS('Table/04_MAG_timeSeries.rds')
ts <- ts/rowSums(ts)
cols <- readRDS('Table/04_color_palette.rds')

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
    ggsave(filename=sprintf('%s/Niche_spapce_PCoA%s_rv.pdf', dir$figdir, i),
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

pdf(sprintf('%s/Niche_spapce_day_rv.pdf', dir$figdir), width=18, height=11)
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
		
	simlist <- statniche <- c();
	for(i in as.character(unique(mag$sample))){ #i=as.character(unique(mag $sample))[5]
		    stmp <- proc.time()[3]
		    ## --------------------------------- ##
		    present <- colnames(ts)[which(ts[i,]>0.001)]
			  sub=bi[present,]
			  sub= sub[,colSums(sub)>0]
			
		    ## =========================== ## 
  			simlist[[i]] <- gene_dist <- as.matrix(vegdist(sub, method='jaccard', binary=TRUE))
  			diag(gene_dist) <- NA

  			mngd <- 1-mean(apply(gene_dist, 1, min, na.rm=TRUE))
  			mpd <- mean(1-gene_dist, na.rm=TRUE)
        
  			mngd.w <- 1-weighted.mean(apply(gene_dist, 1, min, na.rm=TRUE),
  			                          ts[i, colnames(gene_dist)])
  			
  			#nodf <- nestednodf(sub)$statistic[3]
  			#h2 <- H2fun(sub, H2_integer = FALSE)[1]
  			#nodfc <- NODFc(sub, quality=2)			
			  ## =========================== ## 
		    statniche <- rbind(statniche, data.frame(MPD=mpd, MNTD=mngd, `Weighted MNTD`=mngd.w))
		    gc();   gc()
		    
		    etmp <- proc.time()[3]
		    cat(sprintf('Day %s done: elapsed time %0.2f sec\n',i,  etmp - stmp))
	}
#}

## ==================================================== ##

richness <- apply(ts, 1, function(x){ sum(x>0) })
shannon <- diversity(ts, index="shannon")

#red.niche=lm(statniche ~ shannon)$residual

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

abrupts <- cbind('Amplicon based abruptness'= smltarget$abruptness_rel_tw5_tp1, 
                 'Shotgun based abruptness'= abrupt(bray))
lap_score <- data.frame( time=c(1,10,20,24,seq(30,110,10)), 
                         statniche,
                         Richness=richness,
                         Shannon=shannon,
                         abrupts, check.names = FALSE)

## =========================================================== ##
ts_plt <- function(data, exp, abrupt, title){
    
  
    cols <- c("firebrick3", "slateblue", "skyblue", "tan", "darkorange")[1:length(exp)]
    names(cols) <- colnames(data)[exp]
    
    colnames(data)[abrupt] <- "Community\nchange"
    data_Scale <- data
    data_Scale[,-1] <- apply(data_Scale[,-1], 2, scales::rescale, to=c(0,1))
    
    lf <- gather(data_Scale, key, value, c(exp, abrupt))
  
    g <- ggplot(lf, aes(x=time, y=value))+
        geom_line(aes(color=key, group=key, linetype=key))+
        scale_linetype_manual(values=c(1:(length(exp)+1)))+
        scale_color_manual(values=c(cols,
                                    'Community\nchange'='black'),
                           guide=guide_legend(override.aes = 
                                                list(linetype=1:(length(exp)+1),
                                                     linewidth=rep(0.4, length(exp)+1))))+
        labs(y='Rescaled scores', x='Time', color='Score', title=title)+
        #theme_text(bsize=6, family = 'Arial',unitsize = 0.8)+
        theme(legend.key = element_rect(fill=NA))+
        guides(linetype="none")
    return(g)
}    

g1_asv <- ts_plt(lap_score, 2:6, 7, "ASV based")
g1_mag <- ts_plt(lap_score, 2:6, 8, "MAG based")

ggsave(plot=cowplot::plot_grid(g1_asv, g1_mag, ncol=1), sprintf('%s/score_ts_rv.pdf', dir$figdir), w=14, h=12, unit='cm')

#ggsave(plot=cowplot::plot_grid(g1_asv, g1_mag, ncol=1), sprintf('%s/score_ts_rv.pdf', dir$figdir),
       device = cairo_pdf, w=12, h=8, unit='cm')

lf <- gather(lap_score, key, value, c(2:6))
lf2 <- subset(lf, lf$key=='MPD' | lf$key=='MNTD'| lf$key=="Weighted.MNTD")

scat_plt <- function(data, abrupt, title){
  
  colnames(data)[abrupt] <- "Community change"
  g <- ggplot(data, aes(y=`Community change`, x=value, color=time))+
      geom_smooth(method="lm")+
      geom_point()+
      labs(y='Community change', x='', title=title)+
      scale_color_gradientn(colors=brewer.pal(11, 'Spectral'))+
      #theme_text(family = 'Arial')+
      theme(legend.key = element_rect(fill=NA),
            strip.placement = 'outside')+
      facet_wrap(~key, strip.position = "bottom",
                 scales="free")
  return(g)
}

g2_asv <- scat_plt(data=lf2, 2, "ASV based")
g2_mag <- scat_plt(data=lf2, 3, "MAG based")

ggsave(plot=cowplot::plot_grid(g2_asv, g2_mag, ncol=1), 
        sprintf('%s/corerelation_with_abruptness_rv.pdf', dir$figdir), w=19, h=14, unit='cm')

#ggsave(plot=cowplot::plot_grid(g2_asv, g2_mag, ncol=1), 
        sprintf('%s/corerelation_with_abruptness_rv.pdf', dir$figdir),
       device = cairo_pdf, w=14, h=14, unit='cm')

## ====================================== ##
# Shannon vs. niche overlap scores

g1 <- ggplot(lap_score, aes(x=Shannon, y= MPD, color=time))+
        geom_point()+
      scale_color_gradientn(colors=brewer.pal(11, 'Spectral'))+
      geom_smooth(method="lm")+
      theme(legend.key = element_rect(fill=NA),
            strip.placement = 'outside')
            
g2 <- ggplot(lap_score, aes(x=Shannon, y= MNTD, color=time))+
        geom_point()+
      scale_color_gradientn(colors=brewer.pal(11, 'Spectral'))+
      geom_smooth(method="lm")+
      theme(legend.key = element_rect(fill=NA),
            strip.placement = 'outside')

g3 <- ggplot(lap_score, aes(x=Shannon, y=Weighted.MNTD, color=time))+
        geom_point()+
      scale_color_gradientn(colors=brewer.pal(11, 'Spectral'))+
      geom_smooth(method="lm")+
      theme(legend.key = element_rect(fill=NA),
            strip.placement = 'outside')
                   
ggsave(plot=grid.arrange(g1, g2, g3, nrow=1), sprintf('%s/Shannon_vs_nicheoverlap_rv.pdf', dir$figdir), w=22, h=5, unit='cm')

sink(sprintf('%s/Shannon_vs_nicheoverlap_rv.txt', dir$tabledir))
summary(lm(MPD ~ Shannon, data= lap_score))
summary(lm(MNTD ~ Shannon, data= lap_score))
summary(lm(Weighted.MNTD ~ Shannon, data= lap_score))
sink()

## ====================================== ##
lap_score_Scale <- lap_score
lap_score_Scale[,-c(1,7,8)] <- apply(lap_score[,-c(1,7,8)], 2, scale)

cor(lap_score[,-c(1)], use="pairwise.complete.obs")

lm_list <- function(data, res, exp){
  
  fmllist <- c()
  for(i in exp){
    
    fml <- formula( sprintf("`%s` ~ `%s`", colnames(data)[res],
                                       colnames(data)[i]))
    model <- lm(fml, data)
    fmllist[[colnames(data)[i]]] <- model
  }
   
  return(fmllist)
  
}

asv_lmres <- lm_list(lap_score_Scale, 7, 2:6)
mag_lmres <- lm_list(lap_score_Scale, 8, 2:6)

summary_list <- lapply(mag_lmres, summary)
summary_mat <- do.call(rbind, sapply(summary_list, "[", 4))
summary_mat <- summary_mat[grep("Intercept", rownames(summary_mat), invert = TRUE), ]

summary_mat_with_AIC <- cbind(summary_mat, AIC=sapply(mag_lmres, AIC))

write.csv(summary_mat_with_AIC, sprintf("%s/regression_result.csv", dir$tabledir))


###
summary_out <- function(data){
 
  summary_list <- lapply(data, summary)
  summary_mat <- do.call(rbind, sapply(summary_list, "[", 4))
  summary_mat <- summary_mat[grep("Intercept", rownames(summary_mat), invert = TRUE), ]
 
  summary_mat_with_AIC <- cbind(summary_mat, AIC=sapply(mag_lmres, AIC))
 
}
write.csv( summary_out (mag_lmres) , sprintf("%s/regression_result_mag.csv", dir$tabledir))
write.csv( summary_out (asv_lmres) , sprintf("%s/regression_result_asv.csv", dir$tabledir))

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

## ==================================================== ##
# MPD

ggpcasim <- ggplot(pca2)+
         geom_path(aes(x=X1, y= X2), show.legend=FALSE,size=0.4,
                  arrow = arrow(type='closed',length = unit(0.015, "npc")), linetype=3)+
         geom_point(aes(x=X1, y= X2, fill= MPD), shape=21, size=3)+         
         geom_text_repel(data=pca2, aes(x= X1, y= X2, label= sprintf('Day %s', time)), size=2)+		 
         viridis::scale_fill_viridis(guide = guide_colourbar(direction = "horizontal"))+
         scale_color_gradientn(colors=c( rev(brewer.pal(9,'Greys'))[1:2] ) )+
    #theme_text()+
         theme(text=element_text(family='Arial'),
               legend.text=element_text(size=6),
               legend.title=element_text(size=6, vjust=0.8),
               legend.key.width=unit(0.4,"cm"),
               legend.key.height=unit(0.2,"cm"),
               legend.position=c(1, 0),
               legend.justification = c(1,2),
               plot.margin= unit(c(1, 1, 9, 1), "lines"))+
         labs(x='Axis 1',y='Axis 2',fill="Niche overlap")+
    coord_fixed( max(pca2$X1)/max(pca2$X2))

ggsave(plot=ggpcasim, filename=sprintf('%s/PCoA_MPD_rv.pdf', dir$figdir), height=5, width=5)
## ==================================================== ##


## ==================================================== ##
# MNTD


ggpcasim <- ggplot(pca2)+
         geom_path(aes(x=X1, y= X2), show.legend=FALSE,size=0.4,
                  arrow = arrow(type='closed',length = unit(0.015, "npc")), linetype=3)+
         geom_point(aes(x=X1, y= X2, fill= MNTD), shape=21, size=3)+         
         geom_text_repel(data=pca2, aes(x= X1, y= X2, label= sprintf('Day %s', time)), size=2)+		 
         viridis::scale_fill_viridis(guide = guide_colourbar(direction = "horizontal"))+
         scale_color_gradientn(colors=c( rev(brewer.pal(9,'Greys'))[1:2] ) )+
    #theme_text()+
         theme(text=element_text(family='Arial'),
               legend.text=element_text(size=6),
               legend.title=element_text(size=6, vjust=0.8),
               legend.key.width=unit(0.4,"cm"),
               legend.key.height=unit(0.2,"cm"),
               legend.position=c(1, 0),
               legend.justification = c(1,2),
               plot.margin= unit(c(1, 1, 9, 1), "lines"))+
         labs(x='Axis 1',y='Axis 2',fill="Niche overlap")+
    coord_fixed( max(pca2$X1)/max(pca2$X2))

ggsave(plot=ggpcasim, filename=sprintf('%s/PCoA_MNTD_rv.pdf', dir$figdir), height=5, width=5)
## ==================================================== ##


## ==================================================== ##
# Weighted MNTD

ggpcasim <- ggplot(pca2)+
         geom_path(aes(x=X1, y= X2), show.legend=FALSE,size=0.4,
                  arrow = arrow(type='closed',length = unit(0.015, "npc")), linetype=3)+
         geom_point(aes(x=X1, y= X2, fill= Weighted.MNTD), shape=21, size=3)+         
         geom_text_repel(data=pca2, aes(x= X1, y= X2, label= sprintf('Day %s', time)), size=2)+		 
         viridis::scale_fill_viridis(guide = guide_colourbar(direction = "horizontal"))+
         scale_color_gradientn(colors=c( rev(brewer.pal(9,'Greys'))[1:2] ) )+
    #theme_text()+
         theme(text=element_text(family='Arial'),
               legend.text=element_text(size=6),
               legend.title=element_text(size=6, vjust=0.8),
               legend.key.width=unit(0.4,"cm"),
               legend.key.height=unit(0.2,"cm"),
               legend.position=c(1, 0),
               legend.justification = c(1,2),
               plot.margin= unit(c(1, 1, 9, 1), "lines"))+
         labs(x='Axis 1',y='Axis 2',fill="Niche overlap")+
    coord_fixed( max(pca2$X1)/max(pca2$X2))

ggsave(plot=ggpcasim, filename=sprintf('%s/PCoA_Weighted.MNTD_rv.pdf', dir$figdir), height=5, width=5)
## ==================================================== ##

write.table(lap_score, file="lap_score.txt", quote=F, sep='\t')
