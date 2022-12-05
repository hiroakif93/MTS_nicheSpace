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
source('functions/TFO.R')
library(AnalysisHelper)
load.lib( c('ggplot2', 'RColorBrewer', 'ggnetwork',  'igraph', 'cowplot','graphlayouts', 
            'extrafont', 'vegan', 'tidyr', 'plot3D'))

# -- make directory to save results
dir <- make.dir('MAG_metagenome/SMETANA_network')

# -- Load data table
crosfeed <- read.table('Table/ShotgunMetagenome/CarveMe2Smetana/smetana/detailed.tsv', header=T)
crosfeed <- apply(crosfeed, 2, function(x){ gsub('_rename', '', x) } )

mag <- readRDS('Table/ShotgunMetagenome/MAG_infomation.rds') 
mag[mag ==''] <- 'Unidentified'
repmag <- mag[unique(mag$represent),]

ts <- readRDS('Table/04_MAG_timeSeries.rds')
ts <- ts/rowSums(ts)

cols <- readRDS('Table/04_color_palette.rds')
coord <- readRDS('Table/04_06_coord.rds')

############################################################################

stat <- central <- c()

for(i in sort(as.character( unique(crosfeed[,1]) )) ){ #i=as.character( unique(crosfeed[,1]) )[2]
    
    crossub <- crosfeed[ crosfeed[,1]==i, ]
    
    el <- as.matrix(unique(crossub[ ,c(4,3, 6)])) 
    el <- gsub('_prokka', '', el) 
    
    graph <- graph_from_edgelist(el[,1:2])
    E(graph)$weight <- as.numeric(el[,3])
    topology= c(Treeness=T_FUNC(graph), Orderability=O_FUNC(graph), Feedfowardness=F_FUNC(graph))
    
    stat <- rbind(stat, c(time=as.numeric(gsub('MTS-','',i)),
                          richness=length(V(graph)), 
                          topology,
                          mean.SCS= mean( as.numeric(crossub[,6])) ) )
    day <- as.numeric( strsplit(i, '-')[[1]][2])
	
    adj <- as.matrix(as_adj(graph, attr='weight'))
    central <- rbind(central, cbind(time=paste('Day', day),
    				data.frame(bin=V(graph)$name,
    					  deg.in=degree(graph, mode='in', normalized =TRUE), deg.out=degree(graph, mode='out', normalized =TRUE), 
    					  deg.in.weight= colSums(adj), deg.out.weight= rowSums(adj), 
    					  deg.in.by.out=degree(graph, mode='in')/(degree(graph, mode='out')+0.000001),
      					  close.in=closeness(graph, mode="in", normalized =TRUE), close.out=closeness(graph, mode="out", normalized =TRUE),
      					  eigen=eigen_centrality(graph, directed=T, weights=NA)$vector,
      					  betweennes=betweenness(graph, directed=TRUE, normalized =TRUE),
      					  influence=log(influence(as.matrix(as_adj(graph)), mode='positive')[,1]),
      					  pagerank=log(page.rank(graph, weights=NA)$vector ),
      					  influence.weght=(influence(as.matrix(as_adj(graph, attr='weight')), mode='positive')[,1]),
      					  pagerank.weght=(page.rank(graph)$vector )) ))
    
    write.table(as.matrix(as_adj(graph, attr='weight')), sprintf('%s/adj_%s.txt', dir$tabledir, day), 
    sep='\t', quote=FALSE)

}

shape=c(21:25, 4)
names(shape)=sort(unique(mag$Class))

cent2 <- cbind(mag[as.character(central$bin),],central)
cent2$time <- factor(cent2$time, levels=unique(cent2$time))

centsub <- (cent2[,c(2,5:20, 33,34)])
write.table(centsub, 'Table/centrality.txt', row.names=FALSE, sep='\t', quote=FALSE)

centsub$influence <- round(centsub$influence, digits=10)
write.table(centsub, 'Table/centrality_digit10.txt', row.names=FALSE, sep='\t', quote=FALSE)

g <- ggplot(cent2)+
    geom_point(aes(x= (ifluence.weght), y= (pagerank.weght), fill=Genus, shape=Class), size=4)+
    facet_wrap(~time, scales='free', nrow=3)+
    theme_text(bsize=6)+
    scale_fill_manual(values= cols[as.character(cent2 $Genus)] )+
    theme(text=element_text(family='Arial'),
          strip.background = element_blank())+  
    scale_shape_manual(values= shape)+
    labs(x='Influence (logarithm)', y='Pagerank (logarithm)')+
    facet_wrap(~time, scales='free')
        
ggsave(plot=g, w=18,h=7,
       filename = sprintf('%s/infulence_pagerank.tiff', dir$figdir))
       
g <- ggplot(cent2)+
    geom_point(aes(x= ifluence.weght, y= pagerank.weght, fill=Genus, shape=Class), size=4)+
    theme_bw(base_size=14)+
    scale_fill_manual(values= cols[as.character(cent2 $Genus)] )+
    theme(text=element_text(family='Arial'))+  
    scale_shape_manual(values= shape)+
    labs(x='Influence (logarithm)', y='Pagerank (logarithm)')
        
ggsave(plot=g, w=7,h=4,
       filename = sprintf('%s/infulence_pagerank2.tiff', dir$figdir))

g <-ggplot(as.data.frame(stat))+
	geom_line(aes(x=time, y=Orderability), size=1, color=brewer.pal(11,'RdBu')[11])+
    geom_point(aes(x=time, y=Orderability), size=3, color=brewer.pal(11,'RdBu')[11])+
    geom_line(aes(x=time, y=Treeness), size=1, color=brewer.pal(11,'RdBu')[2])+
    geom_point(aes(x=time, y=Treeness), size=3, color=brewer.pal(11,'RdBu')[2])+
    geom_line(aes(x=time, y=Feedfowardness), size=1, color=brewer.pal(11,'RdBu')[4])+
    geom_point(aes(x=time, y=Feedfowardness), size=3, color=brewer.pal(11,'RdBu')[4])+
    theme_minimal(base_size=18)+
    theme(text=element_text(family='Arial'),
          #axis.title=element_text(size=15),
          legend.position = 'bottom',
          panel.background = element_rect(fill=NA, size=0.8))+
    labs(x='Day', y='Orderability', color='')

scatter3D(x=stat[,'Treeness'], z=stat[,'Orderability'], y=stat[,'Feedfowardness'],
		  colvar=stat[,'time'], theta=30, phi=0, type='b', pch=19)
    
ggsave(plot=g, w=6, h=4,
       filename=sprintf('%s/topology.tiff', dir$figdir),)

############################################################################

glist <- c()
for(i in as.character( unique(crosfeed[,1]) ) ){ #i=as.character( unique(crosfeed[,1]) )[1]
    
    crossub <- crosfeed[ crosfeed[,1]==i, ]
    
    el <- as.matrix(unique(crossub[ ,c(4,3, 6)])) 
    
    graph <- graph_from_edgelist(el[,1:2])
    E(graph)$weight <- as.numeric(el[,'scs'])

    ## ================================== ##
    adj <- as.matrix(as_adj(graph, attr='weight'))
    
    #gnet <- ggnetwork(graph, layout= layout )
    
    gnet <- ggnetwork(graph, arrow.gap=0.05, layout= coord[V(graph)$name,1:2])
    
    ggdf <- cbind(gnet, repmag[as.character(gnet$name), ])
    day<- as.numeric(strsplit(as.character(i), '-')[[1]][2])
	
    logweight <- log(as.numeric(el[,3])) 
    d <- ggplot(ggdf, aes(x=x, y=y))+
        geom_edges(aes(xend=xend, yend=yend, color=weight),alpha=0.6,size=0.8,
        		   arrow = arrow(length = unit(8, "pt"), type = "closed"), show.legend= FALSE)+
        geom_nodes(aes(fill=Genus, shape=Class), size=5, show.legend = FALSE)+
        scale_color_gradientn(colors=rev(brewer.pal(11,'RdBu')[-c(5:11)]), na.value=NA, 
        					   limits=c(0,1))+
        scale_fill_manual(values= cols[as.character(ggdf$Genus)] )+
        scale_size_area(max_size = 6)+
        theme_void(base_size=21)+  
        theme(text=element_text(family='Arial'))+  
        xlim(c(-0.1, 1.1))+ylim(c(-0.1,1))+
        labs(subtitle=paste('Day', day))+
        scale_shape_manual(values= shape)
    
    ## ================================== ##
    
    glist[[i]] <- d
}

ggsave(filename=sprintf('%s/metabolite_network_layout_pca.tiff', dir$figdir),
       plot=plot_grid(plotlist=glist, nrow=2), h=8, w=20)

############################################################################
