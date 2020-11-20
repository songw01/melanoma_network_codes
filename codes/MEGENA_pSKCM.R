rm(list = ls())

library(MEGENA)

root.dir <- "/sc/orga/projects/zhangb03a/songw01/Melanoma";setwd(root.dir)
data.file <- "./processed_data/SKCM_primary__gene_expression.txt"
wkdir <- paste("MEGENA",gsub("\\.txt","",gsub("^(.*)/","",data.file)),sep = "_")

# input parameters
n.cores <- 12; ## number of cores
doPar <- TRUE; ## parallelize? 
method = "pearson" ## correlation method
FDR.cutoff = 0.05 ## correlation significance cutoff
module.pval = 0.05 ## module p-value
hub.pval = 0.05 ## hub p-value

annot.col <- 1 ## how many columns in the data are annotationa for genes? 
###########
# create output folder
dir.create(wkdir)

# load data
data.Df <- read.delim(file = data.file,sep = "\t",header = T)
datExpr <- as.matrix(data.Df[,-annot.col]);
rownames(datExpr) <- as.character(data.Df[[1]]);rm(data.Df)

# set working folder
setwd(wkdir)

#####
# calculate correlation
ijw <- calculate.correlation(datExpr = datExpr,doPerm = 10,method = method,FDR.cutoff = FDR.cutoff)

if (getDoParWorkers() == 1 & n.cores > 1 & doPar)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}

# calculate PFN
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores)
write.table(el,file = "MEGENA_Network.txt",sep = "\t",row.names = F,col.names = T,quote = F)
rm(ijw)

# do clustering
g <- graph.data.frame(el,directed = F)

MEGENA.output <- do.MEGENA(g,
 mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
 min.size = 10,max.size = vcount(g)/2,
 doPar = TRUE,num.cores = n.cores,n.perm = 100,
 save.output = TRUE)

save(MEGENA.output,file = "MEGENA_output.RData")

quit(save = "no",status = 0)