rm(list = ls())

library(MEGENA)

###########
root.dir <- "C:/Users/songw01/Documents/Melanoma_NatComm_Revision/revision_submission_Oct2020";setwd(root.dir)
source("codes/sources/enrichment_functions.R")

#### parameters
gmt.folder <- "processed_data/MSigDB_v5";
wkdir <- "processed_data/MEGENA_Results/pSKCM"
MEGENA.file <- "MEGENA_output.RData"
net.file <- "MEGENA_Network.txt"

# annotation attributes
annot.file <- NULL
symbol.col <- 2
id.col <- 1

# module/hub significance
mod.pvalue = 0.05
hub.pvalue = 0.05
min.size = 50
max.size = 5000

## indicators to perform the following down-stream analyses
# MsigDB enrichment analyses
do.MSigDB <- TRUE

file.pttrn <- ".v5.0.symbols_FET-Table.txt"
fpval.cutoff <- 0.05 # significance threshold
fpval.col <- "FET_pvalue" # column name assigning significance p-value for each test.

######################
# load gmt 
gmt.files <- list.files(path = gmt.folder,pattern = "\\.gmt$",full.names = T)

# work on wkdir
cat(paste("Processing:",wkdir,"\n",sep = ""))

# load MEGENA results
load(paste0(wkdir,"/",MEGENA.file))
PFN <- graph.data.frame(read.delim(file = paste0(wkdir,"/",net.file),sep = "\t",header = T),directed = FALSE)

# summarize module level results and get significant output
cat("Summarize modules...\n")
output <- MEGENA.ModuleSummary(MEGENA.output,
                               mod.pvalue = mod.pvalue,hub.pvalue = hub.pvalue,
                               min.size = min.size,max.size = max.size,
                               annot.table = NULL,id.col = 1,symbol.col = 2,
                               output.sig = TRUE)

summary.table <- output$module.table
names(output$modules) =gsub("comp1_","M",names(output$modules))

###### output modules
output.geneSet.file(output$modules,paste0(wkdir,"/multiscale_significant.modules.txt"))

###### perform MSigDB analysis
if (do.MSigDB & !is.null(gmt.folder))
{
  cat("MSigDB signature enrichments...\n")
  saveto <- paste0(wkdir,"/MSigDB");dir.create(saveto)
  
  gmt.files <- list.files(path = gmt.folder,pattern = "\\.gmt$",full.names = T)
  
  bg <- V(PFN)$name
  bg <- unique(gsub("\\|(.*)","",bg))
  
  gs.modules <- lapply(output$modules,function(x) unique(gsub("\\|(.*)","",x)))
  
  output.status <- run.MSigDB(module.geneSets = gs.modules,GeneSymbol = bg,MSigDB.files = gmt.files,saveto = saveto,
                              annot.table = NULL,id.col = id.col,symbol.col = symbol.col,do.multicore = F,n.cores = NULL)
  
  MSigDB.table <- format.FET.output(FET.folder = saveto,pvalue.cutoff = 0.05,fold.cutoff = 2)
  
  # make finalized table with MSigDB summary
  summary.table <- combine.table(abl = summary.table,bbl = MSigDB.table)
  write.table(summary.table,file = paste0(wkdir,"/multiscale_module_summary.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
}else{
  summary.table = read.delim(file = paste0(wkdir,"/multiscale_module_summary.txt"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  # apply simpler module names
  summary.table$id = gsub("^comp1_","M",summary.table$id)
  summary.table$module.parent = gsub("^comp1_","M",summary.table$module.parent)
  
  # subset by sizes
  summary.table = subset(summary.table,module.size >= min.size & module.size <= max.size)
}


###### combine summary table with Module Differential Connectivity Results
# load MDC results (see codes in folder MDC for running this analysis)
MDC.res = read.delim(file = "processed_data/MDC/primary.MDC_module_VS_GTEX_MDC_wFDR.xls",sep= "\t",header = TRUE,stringsAsFactors = FALSE)
MDC.res[[1]] = gsub("^comp1_","M",MDC.res[[1]])

summary.table$GTEx.MDC = MDC.res$MDC[match(summary.table[[1]],MDC.res[[1]])]
summary.table$GTEx.MDC.FDR = MDC.res$FDR[match(summary.table[[1]],MDC.res[[1]])]

###### Run survival signature enrichments
if (TRUE)
{
  surv.files <- c("processed_data/Survival/survival_stat.pSKCM.txt",
                  "processed_data/Survival/survival_stat.mSKCM.txt")
  
  # get common survival signature from both of metastatic and primary tumors
  surv.tbl <- lapply(surv.files,function(x) read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE))
  sigset <- list(common_poor_both = Reduce("intersect",lapply(surv.tbl,function(x) x[[1]][which(x$KM.p.value < 0.05)])),
                 primary_poor_up = subset(surv.tbl[[1]],coef > 0 & KM.p.value < 0.05)[[1]],
                 primary_poor_down = subset(surv.tbl[[1]],coef < 0 & KM.p.value < 0.05)[[1]],
                 metastatic_poor_up = subset(surv.tbl[[2]],coef > 0 & KM.p.value < 0.05)[[1]],
                 metastatic_poor_down = subset(surv.tbl[[2]],coef < 0 & KM.p.value < 0.05)[[1]])
  
  sigset$common_poor_up <- intersect(sigset$primary_poor_up,sigset$metastatic_poor_up)
  sigset$common_poor_down <- intersect(sigset$primary_poor_down,sigset$metastatic_poor_down)
  
  # get background
  bg <- V(PFN)$name
  #bg <- unique(gsub("\\|(.*)","",bg))
  
  # run enrichments
  res.df <- perform.AllPairs.FET(geneSets1 = output$modules,geneSets2 = sigset,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
  res.df <- res.df[order(res.df$FET_pvalue),]
  
  # write results to file
  write.table(res.df,file = paste0(wkdir,"/survival_signature.FET.txt"),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
}else{
  res.df = read.delim(file = paste0(wkdir,"/survival_signature.FET.txt"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
}

# update enrichment results into the module table
library(reshape2)
pmat = acast(data= res.df,formula= set1_Name ~ set2_Name,value.var = "FET_pvalue",fun.aggregate = function(x) min(x,na.rm = T))
pmat[is.infinite(pmat)] = 1
colnames(pmat) = paste0("FET.P_",colnames(pmat))
summary.table = cbind.data.frame(summary.table,as.data.frame(pmat[match(summary.table[[1]],rownames(pmat)),]))

##### Now, rank the table
order.table = summary.table[order(summary.table$FET.P_common_poor_both),]
order.table$module.rank <- 1:nrow(order.table)

##### Refine 
screen.skin <- function(module.table,cutoff = 0.01)
{
  ii <- which(module.table$GTEx.MDC.FDR < cutoff)
  new.table <- module.table[ii,]
  return(new.table)
}
refine.rank.v2 <- function(module.table)
{
  require(igraph)
  module.table$id <- gsub("^comp[0-9]_","M",module.table$id)
  module.table$module.parent <- gsub("^comp[0-9]_a","M",module.table$module.parent)
  mtbl.o = module.table
  # extract module hierarchy
  module.hierarchy <- graph.data.frame(module.table[,c("module.parent","id")],directed = TRUE)
  V(module.hierarchy)$label <- V(module.hierarchy)$name
  
  run.next = TRUE
  r = 1
  while (run.next)
  {
    cat("# number of modules:");cat(nrow(module.table));cat("\n")
    # extract children modules and ranking
    row.rank <- module.table$module.rank[r]
    child.module <- V(module.hierarchy)$name[neighborhood(module.hierarchy,order = 5,nodes = which(V(module.hierarchy)$name == module.table[[1]][r]),mode = "out")[[1]]];
    child.module <- setdiff(child.module,module.table[[1]][r])
    child.rank <- module.table$module.rank[match(child.module,module.table[[1]])]
    
    parent.module <- V(module.hierarchy)$name[neighborhood(module.hierarchy,order = 5,nodes = which(V(module.hierarchy)$name == module.table[[1]][r]),mode = "in")[[1]]];
    parent.module <- setdiff(parent.module,module.table[[1]][r])
    parent.rank <- module.table$module.rank[match(parent.module,module.table[[1]])]
    
    eliminate.module <- c(child.module[which(child.rank > row.rank)],parent.module[which(parent.rank > row.rank)])
    cat(paste("-- eliminating:",paste(eliminate.module,collapse = ","),"\n",sep = ""))
    module.table <- rbind.data.frame(module.table[1:r,],subset(module.table[(r+1):nrow(module.table),],!(id %in% eliminate.module)))
    
    r = r + 1
    
    if (r >= (nrow(module.table)-1)) run.next = FALSE
  }
  
  list(ranked = mtbl.o,refined = module.table)
  #return(module.table)
}


order.table <- screen.skin(order.table,cutoff = 0.05)
rank.res <- refine.rank.v2(order.table)

# output the results to a file
write.table(rank.res$refined,file = paste0(wkdir,"/module_summary_table.refine_ranked.txt"),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

