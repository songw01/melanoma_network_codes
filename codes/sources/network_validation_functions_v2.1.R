anyNA <- function(x) any(is.na(x))
factorwise.pvalue.correct <- function(output.Df,fact.cols,adj.method = "bonferroni")
{
	##### inputs
	# output.Df = FET output (data.frame object) containing at least "FET_pvalue" and "corrected.FET.pvalue" along with fact.cols
	# fact.cols = vector of character strings specifying column names of respective factor variables
	# adj.method = parameter is passed to "p.adjust" function to specify "method". 
	
	##### correct p-values
	# split row indices according to the factors given in fact.cols
	fact.i <- match(fact.cols,colnames(output.Df))

	# check if columns are factors:
	is.fact <- sapply(output.Df[,fact.i],is.factor)
	if (!all(is.fact)) stop("some column(s) is (are) not factor variables.")

	index.lst <- list(1:nrow(output.Df))
	for (ff in fact.i)
	{
	 new.Lst <- list()
	 for (jj in 1:length(index.lst))
	 {
	  new.Lst <- c(new.Lst,split(index.lst[[jj]],factor(output.Df[[ff]][index.lst[[jj]]])))
	 }
	 index.lst <- new.Lst;rm(new.Lst)
	}

	# correct p-values within each set of indices of index.lst[[i]]
	index.lst <- index.lst[which(sapply(index.lst,length) > 0)]
	cor.pval <- lapply(index.lst,function(x,y,mth) p.adjust(y[x],method = mth),y = output.Df$FET_pvalue,mth = adj.method)
	ip <- cbind(do.call('c',index.lst),do.call('c',cor.pval))
	output.Df$corrected.FET.pvalue[ip[,1]] <- ip[,2]

	new.Df <- output.Df
	return(new.Df)
}

######
genewise.fc.foldSummary <- function(FET.tbl)
{
 plot.obj <- ggplot(data = FET.tbl,aes(x = neighbor.layer,y = enrichment.foldchange,fill = DEG.ID)) + geom_bar(position = "dodge",stat = "identity") + facet_grid(fold.change ~ exp.cond) +
 scale_fill_manual(values = c("DN" = "blue","UP" = "red","ALL" = "green")) + theme_bw() + labs(x = "neighborhood layer",y = "enrichment fold change",title = as.character(FET.tbl$source.gene[1])) +
 geom_hline(aes(yintercept = 2,colour = "red",linetype = "dashed")) + 
 theme(axis.text = element_text(size = 15),plot.title = element_text(size = 20),axis.title = element_text(size = 17),strip.text = element_text(size = 15),
 legend.text = element_text(size = 15),legend.title = element_text(size = 18))
 return(plot.obj)
}

genewise.FETpvalue.foldSummary <- function(FET.tbl)
{
 plot.obj <- ggplot(data = FET.tbl,aes(x = neighbor.layer,y = -log10(FET_pvalue),fill = DEG.ID)) + geom_bar(position = "dodge",stat = "identity") + 
 scale_fill_manual(values = c("DN" = "blue","UP" = "red","ALL" = "green")) + theme_bw() + labs(x = "neighborhood layer",y = "-log10(FET p-value)",title = as.character(FET.tbl$source.gene[1])) +
 #geom_hline(aes(yintercept = -log10(cut.off),colour = "red")) + 
 theme(axis.text = element_text(size = 15),plot.title = element_text(size = 20),axis.title = element_text(size = 17),strip.text = element_text(size = 15),
 legend.text = element_text(size = 15),legend.title = element_text(size = 18))
 return(plot.obj)
}

genewise.overlap.foldSummary <- function(FET.tbl)
{
 plot.obj <- ggplot(data = FET.tbl,aes(x = neighbor.layer,y = actual.overlap,fill = DEG.ID)) + geom_bar(position = "dodge",stat = "identity") + 
 scale_fill_manual(values = c("down" = "blue","up" = "red","all" = "green")) + theme_bw() + labs(x = "neighborhood layer",y = "overlap",title = as.character(FET.tbl$source.gene[1])) +
 #geom_hline(aes(yintercept = -log10(cut.off),colour = "red")) + 
 theme(axis.text = element_text(size = 15),plot.title = element_text(size = 20),axis.title = element_text(size = 17),strip.text = element_text(size = 15),
 legend.text = element_text(size = 15),legend.title = element_text(size = 18))
 return(plot.obj)
}


genewise.FETpvalue.ggplot <- function(FET.tbl)
{
 
 plot.obj <- ggplot(data = FET.tbl,aes(x = neighbor.layer,y = -log10(FET_pvalue),fill = DEG.ID)) + geom_bar(position = "dodge",stat = "identity") + 
 scale_fill_manual(values = c("DN" = "blue","UP" = "red","ALL" = "green")) + theme_bw() + labs(x = "neighborhood layer",y = "-log10(FET p-value)",title = as.character(FET.tbl$source.gene[1])) +
 #geom_hline(aes(yintercept = -log10(cut.off),colour = "red")) + 
 theme(axis.text = element_text(size = 15),plot.title = element_text(size = 20),axis.title = element_text(size = 17),strip.text = element_text(size = 15),
 legend.text = element_text(size = 15),legend.title = element_text(size = 18))
}

genewise.fc.ggplot <- function(FET.tbl)
{
 plot.obj <- ggplot(data = FET.tbl,aes(x = neighbor.layer,y = enrichment.foldchange,fill = DEG.ID)) + geom_bar(position = "dodge",stat = "identity") +  scale_fill_manual(values = c("DN" = "blue","UP" = "red","ALL" = "green")) + theme_bw() + labs(x = "neighborhood layer",y = "enrichment fold change",title = as.character(FET.tbl$source.gene[1])) +
 geom_hline(aes(yintercept = 2,colour = "red",linetype = "dashed")) + 
 theme(axis.text = element_text(size = 15),plot.title = element_text(size = 20),axis.title = element_text(size = 17),strip.text = element_text(size = 15),
 legend.text = element_text(size = 15),legend.title = element_text(size = 18))
}

genewise.overlap.ggplot <- function(FET.tbl)
{
 plot.obj <- ggplot(data = FET.tbl,aes(x = neighbor.layer,y = actual.overlap,fill = DEG.ID)) + geom_bar(position = "dodge",stat = "identity") + 
 scale_fill_manual(values = c("DN" = "blue","UP" = "red","ALL" = "green")) + theme_bw() + labs(x = "neighborhood layer",y = "overlap",title = as.character(FET.tbl$source.gene[1])) +
 #geom_hline(aes(yintercept = -log10(cut.off),colour = "red")) + 
 theme(axis.text = element_text(size = 15),plot.title = element_text(size = 20),axis.title = element_text(size = 17),strip.text = element_text(size = 15),
 legend.text = element_text(size = 15),legend.title = element_text(size = 18))
}

##########


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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

 if (numPlots==1) {
    print(plots[[1]])

  } else {
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
}


###########################################gene
require(reshape)
require(igraph)
require(ggplot2)

run.neighborhood.KDA <- function(networks,signatures,annot.df,out.dir = NULL,mode = "out",max.layer = 5,max.prop = 0.5,min.size = 10,background = NULL,output.progress = FALSE)
{
	# networks = list of igraph objects. Each element of list is properly named. 
	# signatures = list of gene signatures. Each element of list is properly named and corresponds to first column of annot.df.
	# annot.df = signature annotation table. data.frame object containing "id" (names of signatures) and "target.gene" (target gene of neighborhood expansion). 
	# out.dir = a subdirectory to output 
	# mode = direction of neighborhood expansion.
	
	if (is.null(background)) background <- vector("list",length(networks))

	if (is.null(out.dir)) out.dir <- "KDA_output"
	dir.create(out.dir)
	
	# make sure signatures and annot.df are aligned
	common.id <- intersect(names(signatures),as.character(annot.df$id))
	signatures <- signatures[common.id]
	annot.df <- annot.df[match(common.id,as.character(annot.df$id)),]
	annot.df$target.gene <- factor(annot.df$target.gene)
	
	# segregate gene sets per target gene by index
	sig.geneLst <- split(1:nrow(annot.df),annot.df$target.gene)
	
	global.table <- data.frame()
	for (n in 1:length(networks))
	{
	 
	 # prepare inputs 
	 cat(paste("# Processing:",names(networks)[n],"\n",sep = ""))
	 output.header <- paste(names(networks)[n],"_KDA",sep = "")
	 network <- networks[[n]]
	 int.genes <- intersect(levels(factor(annot.df$target.gene)),V(network)$name)
	 
	 # obtain neighborhood
	 cat("- neighborhood calculation...\n")
	 neigh.out <- get.neighborhood(network,int.genes,mode = mode,max.layer = max.layer,max.prop = max.prop)
	 neigh.out <- do.call(c,neigh.out);
	 names(neigh.out) <- gsub(".n.layer_","_l",names(neigh.out))
	 
	 # kick-out small neighborhoods < min.size
	 neigh.out <- neigh.out[which(sapply(neigh.out,length) >= min.size)]
	 
	 # match pairs 
	 neigh.Lst <- split(1:length(neigh.out),factor(gsub("_l(.*)$","",names(neigh.out))))
	 common.gene <- intersect(names(neigh.Lst),names(sig.geneLst))
	 
	 ij.pair <- matrix(0,nrow = 0,ncol = 2)
	 for (g in common.gene)
	 {
	  xy <- do.call('rbind',lapply(sig.geneLst[g][[1]],function(x,y) cbind(rep(x,length(y)),y),y = neigh.Lst[g][[1]]))
	  ij.pair <- rbind(ij.pair,xy);rm(xy)
	 }
	 
	 # output summary
	 
	 output.str <- paste("network:",names(networks)[n],"\n|V| = ",vcount(networks[[n]]),",|E| = ",ecount(networks[[n]]),"\n",
	 "genes with DEG signature:",paste(levels(annot.df$target.gene),collapse = ","),"\n",
	 "genes present in network:",paste(int.genes,collapse = ","),"\n",sep = "")
	 neigh.str <- paste(paste(names(neigh.out),sapply(neigh.out,length),sep = ":"),collapse = "\n")
	 signature.str <- paste(paste(names(signatures),sapply(signatures,length),sep = ":"),collapse = "\n")
	 
	 if (output.progress)
	 {
		 sink(paste("./",out.dir,"/",output.header,".inputSummary.txt",sep = ""))
		 cat(output.str)
		 cat("\n")
		 cat("Tested neighborhoods:\n")
		 cat(neigh.str)
		 cat("\n")
		 cat("Tested signatures:\n")
		 cat(signature.str)
		 sink()
	 }
	 
	 if (is.null(background[[n]])) 
	 {
	  #bg <- Reduce("union",lapply(networks,function(g) V(g)$name))
	  bg <- V(network)$name
	 }else{
	  bg <- background[[n]]
	 }
	 
	 
	 FET.output <- perform.ijPairs.FET(geneSets1 = signatures,geneSets2 = neigh.out,ij = ij.pair,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = 4) 
	 target.gene <- as.character(annot.df$target.gene[match(as.character(FET.output[[1]]),annot.df$id)])
	 source.gene <- gsub("_(.*)$","",as.character(FET.output[[2]]))
	 FET.output <- FET.output[which(target.gene == source.gene),];
	 FET.output$corrected.FET.pvalue <- p.adjust(FET.output$FET_pvalue,"BH");
	 colnames(FET.output)[1:2] <- c("signature.id","neighborhood.id")
	 FET.output <- data.frame(network = rep(names(networks)[n],nrow(FET.output)),FET.output)
	 
	 global.table <- rbind.data.frame(global.table,FET.output)
	}
	
	split.neigh.id <- as.data.frame(do.call(rbind,strsplit(as.character(global.table$neighborhood.id),"_")));colnames(split.neigh.id) <- c("source.gene","neighbor.layer");
	global.table <- cbind.data.frame(global.table[,1:3],split.neigh.id,global.table[,4:ncol(global.table)])
	global.table$neighbor.layer <- as.integer(gsub("^l","",as.character(global.table$neighbor.layer)));
	# output results
	#write.table(global.table,file = paste(".",out.dir,"KDA_Enrichment.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)

	return(global.table)
}


get.neighborhood <- function(g,nodes,mode = "all",max.layer = 10,max.prop = 0.5)
{
 ######## inputs
 # g = igraph object,
 # nodes = root nodes to expand neighborhoods from,
 # mode = choose from c("all","in","out"): specify the traversing direction
 # max.layer = the maximum number of layers to traverse. 
 neigh.layers <- vector("list",length(nodes));names(neigh.layers) <- nodes;
 for (i in 1:length(nodes))
 {
  neigh.out <- lapply(1:max.layer,function(n,g,v,mm) {out <- neighborhood(g,order = n,nodes = v,mode = mm)[[1]];out <- V(g)$name[out];return(out)},g = g,v = nodes[i],mm = mode)
  
  # identify maximum layer that traversed less than 90% of all nodes. 
  mx.layer <- max(c(1,max(which(sapply(neigh.out,length) <= (max.prop*vcount(g)))))) 
  neigh.out <- neigh.out[1:mx.layer]
  # identify maximum layer that reaches the maximum number of neighborhood
  neigh.size <- sapply(neigh.out,length)
  max.num <- min(which(neigh.size == max(neigh.size)));
  neigh.out <- neigh.out[1:max.num]
  
  neigh.layers[[i]] <- neigh.out;
  names(neigh.layers[[i]]) <- paste("n.layer",1:length(neigh.out),sep = "_")
 }
 return(neigh.layers)
} 

read.geneSet <- function(geneSet.file)
{
 gene.list <- readLines(geneSet.file)
 gene.list <- strsplit(gene.list,"\t")
 
 names(gene.list) <- sapply(gene.list,function(x) x[1])
 gene.list <- lapply(gene.list,function(x) x[2:length(x)])
 return(gene.list)
}

output.geneSet.file <- function(geneSet,outputfname)
{
 if (!is.list(geneSet)) stop("geneSet is not a list.")
 if (is.null(names(geneSet))) stop("names(geneSet) is not defined properly.")
  
 sink(outputfname)
 cat(paste(paste(names(geneSet),"\t",sapply(geneSet,function(x) paste(x,collapse = "\t")),sep = ""),collapse = "\n"))
 sink()

 return(0)
}

############# get tables
# all pairs
make.Pairwise.Tables <- function(geneSets1,geneSets2,background)
{
 mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 d <- t(mem1) %*% mem2;
 b <- abs(t(mem1) %*% (mem2-1))
 c <- abs(t(mem1-1) %*% (mem2))
 a <- t(mem1-1) %*% (mem2-1);

 ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

 pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}

# selected pairs by ij
make.paired.Tables <- function(geneSets1,geneSets2,ij,background)
{
 
 #mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 #mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 pairwise.tables <- mapply(FUN = function(s1,s2,z) {
                                         v1 <- rep(0,length(z));v1[which(z %in% s1)] <- 1; 
										 v2 <- rep(0,length(z));v2[which(z %in% s2)] <- 1;
										 # n11,n12,n21,n22
										 as.table(matrix(c(sum(abs(v1 - 1) * abs(v2 - 1)),sum(abs(v2 - 1) * v1),sum(abs(v1 - 1) * v2),sum(v1 * v2)),nrow = 2))
                                },s1 = geneSets1[ij[,1]],s2 = geneSets2[ij[,2]],MoreArgs = list(z = background),SIMPLIFY = FALSE)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}

do.FisherExactTest <- function(table.count,N_bg = NULL)
{
 if (is.null(N_bg)) N_bg = sum(rowSums(table.count))
 
 out <- fisher.test(x = table.count,or = 1,alternative = "greater")
 odds.ratio <- out$estimate
 p.value <- out$p.value;
 geneSet1.count <- rowSums(table.count)[2]
 geneSet2.count <- colSums(table.count)[2]
 expected.count <- geneSet1.count/N_bg * geneSet2.count
 overlap.count <- table.count[2,2];
 fold.change <- overlap.count/expected.count
 
 out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
 names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
 return(out)
}

require(parallel)
require(foreach)
require(iterators)
perform.AllPairs.FET <- function(geneSets1,geneSets2,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)
 
 if (do.multicore)
 {
  #registerDoMC(n.cores)
  set.parallel.backend(n.cores)
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  output <- foreach (tbl = split.tables,.combine = 'c') %dopar% {
            out <- lapply(tbl,do.FisherExactTest) 
			return(out)
  }
 
 }else{
  output <- lapply(pairwise.tables,do.FisherExactTest)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 return(output)
}

perform.ijPairs.FET <- function(geneSets1,geneSets2,ij,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 pairwise.tables <- make.paired.Tables(geneSets1,geneSets2,ij,background)
 
 if (do.multicore)
 {
  #registerDoMC(n.cores)
  set.parallel.backend(n.cores)
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  output <- foreach (tbl = split.tables,.combine = 'c') %dopar% {
            out <- lapply(tbl,do.FisherExactTest) 
			return(out)
  }
 
 }else{
  output <- lapply(pairwise.tables,do.FisherExactTest)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 return(output)
}

###### functions for parallelization
split.indices <- function(nr,dn)
{
	 split.id <- do.call(c,lapply(1:ceiling(nr/dn),function(n,dn) rep(n,dn),dn = dn))[1:nr]
	 split.id <- split(1:nr,factor(split.id))
	 
	 if (length(split.id[[length(split.id)]]) <= 1) 
	 {
	  split.id[[length(split.id)-1]] <- c(split.id[[length(split.id)-1]],split.id[[length(split.id)]])
	  split.id <- split.id[1:(length(split.id)-1)]
	 }
	 
	 return(split.id)
}

parallel.neighbor.KDA <- function(networks,geneSet,annot.df,out.dir,
mode = "all",max.layer = 10,max.prop = 0.3,min.size = 10,background = NULL,
doPar = FALSE,n.core = 4,max.per.job = 400)
{
	if (doPar)
	{
	 require(doParallel)
	 
	 # register cores
	 cl <- makeCluster(n.core)
	 registerDoParallel(cl)
	 clusterEvalQ(cl , c(library(igraph),library(foreach)))
	 
	 # first phase of split: split jobs by max.per.job argument
	 njob.split <- split.indices(nr = nrow(annot.df),dn = max.per.job)
	 
	 cat(paste("- Total ",length(njob.split)," batches of ",max.per.job," gene sets for total ",nrow(annot.df)," gene sets\n",sep = ""))
	 # perform parallel computaion per njob.split[[n]]
	 
	 output.df <- data.frame()
	 for (n in 1:length(njob.split))
	 {
	  cat(paste("Processing: #",n," batch of ",length(njob.split[[n]])," gene sets\n",sep = ""))
	  split.id <- lapply(split.indices(nr = length(njob.split[[n]]),dn = ceiling(length(njob.split[[n]])/n.core)),function(ii,nn) nn[ii],nn = njob.split[[n]])
	  cat(paste("Breakdown of gene sets in each core:",paste(sapply(split.id,length),collapse = "/"),"\n",sep = ""))
	  
	  export.func <- c("run.neighborhood.KDA","get.neighborhood","perform.AllPairs.FET","do.FisherExactTest","make.Pairwise.Tables","perform.ijPairs.FET","make.paired.Tables")
	  tr.arg <- system.time(
	  out.df <- foreach(ii = split.id,.combine = 'rbind.data.frame',.export = export.func) %dopar% {
				   out <- run.neighborhood.KDA(networks = networks,signatures = geneSet[ii],annot.df = annot.df[ii,],out.dir = out.dir,mode = mode,max.layer = max.layer,max.prop = max.prop,min.size = min.size,background = background)
				   return(out)
	  }
	  )
	  
	  output.df <- rbind.data.frame(output.df,out.df)
	  print(tr.arg)
	  rm(split.id,out.df)
	 }
	  
	 
	 stopCluster(cl) 
	 
	}else{
	 ###### 3.2 perform analysis
	 output.df <- run.neighborhood.KDA(networks = networks,signatures = geneSet,annot.df = annot.df,out.dir = out.dir,mode = mode,max.layer = max.layer,max.prop = max.prop,min.size = min.size,background = background)
	}
	return(output.df)
}