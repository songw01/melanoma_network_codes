#### check siRNA signature again
rm(list = ls())

root.dir <- "C:/Users/songw01/Documents/Melanoma_NatComm_Revision/revision_submission_Oct2020";setwd(root.dir)
library(RRHO)

###############################################
source("codes/sources/limma_DEG_functions.R")
source("codes/sources/network_validation_functions_v2.1.R")

##### 
data.file <- "processed_data/SKmel147_RNAseq/processed_expression.txt"
annot.cols = 1:2

# network validation parameters
net.files <- "processed_data/MEGENA_Results/pSKCM/MEGENA_Network.txt"
net.id <- c("primary")
names(net.files) <- net.id

# DEG result file
deg.file = "siRNA_Validation/siRNA_signatures.txt"

# DEG results parameters
pval.cutoff = 0.05;
fc.cutoff = 1.2

out.dir <- "siRNA_Validation";dir.create(out.dir)

##################################################
#### call expressions
data.df <- read.delim(file = data.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE,
                      nrows= 5)
data.mat <- as.matrix(data.df[,-annot.cols]);rownames(data.mat) <- data.df[[1]]

# create annotation table for samples
am <- do.call('rbind',strsplit(colnames(data.mat),"_"));colnames(am) <- c("submitter","siRNA.target","replicate.id")
adf <- data.frame(sid = colnames(data.mat),as.data.frame(am,stringsAsFactors = FALSE),stringsAsFactors = FALSE);
adf$siRNA.target <- gsub("^si","",adf$siRNA.target)


#### call networks
# get networks
networks <- lapply(net.files,function(x) graph.data.frame(read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE),directed = FALSE))
networks <- lapply(networks,function(x) {out <- x;V(out)$name <- gsub("\\|(.*)","",V(x)$name);return(out)})

#### per network, get shortest path distance from target node for a measure to test
require(reshape2)
spdf <- data.frame()
for (i in 1:length(networks))
{
  dg <- networks[[i]];
  E(dg)$weight = 1-E(networks[[i]])$weight
  dmat <- distances(graph = dg,v = setdiff(unique(adf$siRNA.target),"NTC"),to = V(dg),algorithm = "dijkstra")
  dmatdf <- melt(dmat)
  colnames(dmatdf) <- c("target.gene","downstream.gene","distance")
  dmatdf$network <- rep(names(networks)[i],nrow(dmatdf))
  spdf <- rbind.data.frame(spdf,dmatdf)
  rm(dg,dmat,dmatdf)
}

#### check target gene expressions &run differential expression analysis
complete.deg <- read.delim(file = deg.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)

######################### RUN RRHO per target per network
siRNA = setdiff(unique(adf$siRNA.target),"NTC")
res <- list()
for (ng in 1:length(networks))
{
  for (ns in 1:length(siRNA))
  {
    gl.net <- subset(spdf,target.gene == siRNA[ns] & network == names(networks)[ng])[,2:3]
    gl.net <- gl.net[!duplicated(gl.net[[1]]),]
    
    gl.deg <- subset(complete.deg,target.gene == siRNA[ns])
    colnames(gl.deg)[1] <- "gene"
    
	gl.deg$RRHO.statistic = -log10(gl.deg$P.Value) * sign(gl.deg$logFC)
    gl.deg <- gl.deg[,c("gene","RRHO.statistic")]
    
    cmbdf <- merge(x = gl.net,y = gl.deg,by.x = "downstream.gene",by.y = "gene");
    cmbdf <- cmbdf[!is.na(cmbdf[[2]]) & !is.na(cmbdf[[3]]),]
    
    # run RRHO
    RRHO.res <-  RRHO(cmbdf[,c(1,2)], cmbdf[,c(1,3)], 
                      BY=TRUE, alternative='enrichment')
    ival <- RRHO.res$stepsize * (1:nrow(RRHO.res$hypermat) )
    jval <- RRHO.res$stepsize * (1:ncol(RRHO.res$hypermat) )
    ival[length(ival)] <- min(c(ival[length(ival)],nrow(cmbdf)))
    jval[length(jval)] <- min(c(ival[length(jval)],nrow(cmbdf)))
    ivec <- signif(sort(cmbdf[[2]],decreasing = TRUE)[ival],3)
    jvec <- signif(sort(cmbdf[[3]],decreasing = TRUE)[jval],3)
    
    require(reshape)
    require(ggplot2)
    hyper.df <- melt(RRHO.res$hypermat)
    hyper.df[[1]] <- as.character(ivec[hyper.df[[1]]])
    hyper.df[[2]] <- as.character(jvec[hyper.df[[2]]])
    pobj <- ggplot(data = hyper.df,aes(x = X1,y = X2,fill = value)) + geom_tile() + 
      scale_x_discrete(limits = as.character(ivec),breaks = as.character(ivec)[seq(0,length(ivec),10)]) + 
      scale_y_discrete(limits = as.character(jvec),breaks = as.character(jvec)[seq(0,length(jvec),10)]) + 
      labs(x = "shortest distance",y = "-log10(p-value) * sign(logFC)") + 
      scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = median(RRHO.res$hypermat)) + 
      theme_bw() + guides(fill = guide_colourbar(title = "-log10(FET P)",
                                                 label.theme = element_text(angle = 45,vjust = 1,hjust = 1))) + 
      theme(legend.position = "bottom",legend.direction = "horizontal",
            axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
    # calculate p-value
    pval.testing <- pvalRRHO(RRHO.res, 50)
    
    ### make ggplot heatmap for hypermat object
    require(reshape)
    
    out <- list(input.data = cmbdf,RRHO.res = RRHO.res,hypermat.heatmap = pobj,
                pval.res = pval.testing,target.gene = siRNA[ns],tested.network = names(networks)[ng])
    res <- c(res,list(out))
  }
}

saveRDS(res,file = paste(out.dir,"/RRHO_results.RDS",sep = ""))

### get some plotting done
get_ggheat <- function(RRHO.res,cmbdf)
{
  require(reshape)
  require(ggplot2)
  ival <- RRHO.res$stepsize * (1:nrow(RRHO.res$hypermat) )
  jval <- RRHO.res$stepsize * (1:ncol(RRHO.res$hypermat) )
  ival[length(ival)] <- min(c(ival[length(ival)],nrow(cmbdf)))
  jval[length(jval)] <- min(c(ival[length(jval)],nrow(cmbdf)))
  ivec <- signif(sort(cmbdf[[2]],decreasing = TRUE)[ival],3)
  jvec <- signif(sort(cmbdf[[3]],decreasing = TRUE)[jval],3)
  
  #rownames(RRHO.res$hyper.mat) <- 1:nrow(RRHO.res$hypermat)
  #colnames(RRHO.res$hyper.mat) <- 1:ncol(RRHO.res$hypermat)
  hyper.df <- melt(RRHO.res$hypermat)
  hyper.df[[1]] <- as.character(ivec[hyper.df[[1]]])
  hyper.df[[2]] <- as.character(jvec[hyper.df[[2]]])
  
  # conversion to log10
  hyper.df$value <- hyper.df$value/log(10)
  hyper.df$value[hyper.df$value > 300] <- 300
  pobj <- ggplot(data = hyper.df,aes(x = X1,y = X2,fill = value)) + geom_tile() + 
    scale_x_discrete(limits = as.character(ivec),breaks = as.character(ivec)[seq(0,length(ivec),10)]) + 
    scale_y_discrete(limits = as.character(jvec),breaks = as.character(jvec)[seq(0,length(jvec),10)]) + 
    labs(x = "shortest distance",y = "-log10(p-value) * sign(logFC)") + 
    scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = median(hyper.df$value)) + 
    theme_bw() + guides(fill = guide_colourbar(title = "-log10(FET P)",
                                               label.theme = element_text(angle = 45,vjust = 1,hjust = 1))) + 
    theme(legend.position = "bottom",legend.direction = "horizontal",
          axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
  
  return(pobj)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    # Get the i,j matrix positions of the regions that contain this subplot
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    
    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                    layout.pos.col = matchidx$col))
    
  }
}

print_several <- function(plot.list,per.pg = 6,cols = 2)
{
  # print out 6 plots per page
  pgn <- ceiling(length(plot.list)/per.pg)
  pg = 0
  layout <- matrix(seq(1, cols * ceiling(per.pg/cols)),
                   ncol = cols, nrow = ceiling(per.pg/cols))
  for (p in 1:pgn)
  {
    pg.f = pg + per.pg
    pg.f <- min(c(pg.f,length(plot.list)))
    pl <- plot.list[(pg+1):pg.f]
    pg = pg.f
    #pl <- lapply(km,function(x) x$kmplot)
    multiplot(plotlist=pl, cols=NULL, layout=layout)
  }
}


plst <- vector("list",length(res))
for (i in 1:length(res))
{
  pobj <- get_ggheat(RRHO.res = res[[i]]$RRHO.res,cmbdf = res[[i]]$input.data)
  pobj <- pobj + labs(title = paste(res[[i]]$tested.network,res[[i]]$target.gene,sep = ":"))
  plst[[i]] <- pobj
  rm(pobj)
}

pdf(file = paste(out.dir,"/RRHO_heatmap.primary.pdf",sep = ""),width = 10,height = 4)
print_several(plot.list = plst,per.pg = 3,cols = 3)
dev.off()

quit(save = "no",status = 0)
