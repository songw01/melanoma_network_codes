rm(list = ls())

###################a
root.dir <- "C:/Users/songw01/Documents/Melanoma_NatComm_Revision/revision_submission_Oct2020";setwd(root.dir)
meth.file <- "CIT/ExprMethyl.matched.Primary.RData"
sig.file <- "CIT/TCGA.SKCM.Primary.Methyl_Gene.spearman.Bonferroni_0.05.tsv"
out.dir <- "CIT"

pval.cutoff = 0.05

#### make sure to install RCit codes first
install.packages("codes/sources/RCit_0.0.1.tar.gz",type = "source",repos = NULL)

#### extract methylation signatures: by correlations
sig.df <- read.delim(file = sig.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)

colnames(sig.df)[1:2] <- c("Variant","Gene")
sig.df$COR.ID[sig.df$r > 0] <- "POS"
sig.df$COR.ID[sig.df$r < 0] <- "NEG"
sig.df$set.id <- paste(sig.df$Variant,sig.df$Type,sig.df$COR.ID,sep = "_")
sig.df$set.id.concensus <- paste(sig.df$Variant,sig.df$COR.ID,sep = "_")
sig.df$set.id.type <- paste(sig.df$Variant,sig.df$Type,sep = "_")

####
load(meth.file)

##################### CIT portion
##### create trio based on the following chain types: segment -> cis -> trans, segment -> trans -> trans

all.segs <- unique(sig.df$Variant)
trio <- matrix(NA,nrow = 0,ncol = 3);
duo <- matrix(NA,nrow = 0,ncol = 2);
colnames(trio) <- c("L","G","T")
reg.type = c()

#### do CIT, anchored on methylations
for (si in 1:length(all.segs))
{
  cis.sdf <- subset(sig.df,Variant == all.segs[si] & Type == "cis")
  trans.sdf <- subset(sig.df,Variant == all.segs[si] & Type == "trans")
  trans.sdf <- subset(trans.sdf,!(Gene %in% cis.sdf$Gene))
  
  # segregate significant couplets
  if (nrow(cis.sdf) == 1)
  {
    vec <- c(all.segs[si],cis.sdf$Gene)
    duo <- rbind(duo,vec)
    rm(vec)
  }
  
  # account cis -> cis chains
  if (TRUE)
  {
    if (nrow(cis.sdf) > 1)
    {
      for (ci in 1:(nrow(cis.sdf)-1))
      {
        for (ti in (ci+1):nrow(cis.sdf))
        {
          vec <- c(all.segs[si],cis.sdf$Gene[ci],cis.sdf$Gene[ti])
          trio <- rbind(trio,vec)
          trio <- rbind(trio,vec[c(1,3,2)])
          reg.type <- c(reg.type,"CC","CC")
          rm(vec)
        }
      }
    }
  }
  # account cis -> trans chains
  if (TRUE)
  {
    if (nrow(cis.sdf) > 0 & nrow(trans.sdf) > 0)
    {
      for (ci in 1:nrow(cis.sdf))
      {
        for (ti in 1:nrow(trans.sdf))
        {
          vec <- c(all.segs[si],cis.sdf$Gene[ci],trans.sdf$Gene[ti])
          trio <- rbind(trio,vec)
          reg.type <- c(reg.type,"CT")
          
          rm(vec)
        }
      }
    }
  }
  
  # account trans -> trans chains
  if (FALSE)
  {
    if (nrow(trans.sdf) > 1)
    {
      for (ci in 1:(nrow(trans.sdf)-1))
      {
        for (ti in (ci+1):nrow(trans.sdf))
        {
          vec <- c(all.segs[si],trans.sdf$Gene[ci],trans.sdf$Gene[ti])
          trio <- rbind(trio,vec)
          trio <- rbind(trio,vec[c(1,3,2)])
          reg.type <- c(reg.type,"TT","TT")
          rm(vec)
        }
      }
    }
  }
  
  rm(cis.sdf,trans.sdf)
}

# get indexed trio
trio.index <- cbind(match(trio[,1],colnames(L)),match(trio[,2],colnames(G)),match(trio[,3],colnames(G)))
trio.m <- trio;
trio.m.index <- trio.index;
reg.type.m <- paste("M",reg.type,sep="")
rm(trio.index,trio)

##### get trans exp -> methylation -> cis
trio <- matrix(NA,nrow = 0,ncol = 3);
colnames(trio) <- c("L","G","T")
reg.type = c()
for (si in 1:length(all.segs))
{
  cis.sdf <- subset(sig.df,Variant == all.segs[si] & Type == "cis")
  trans.sdf <- subset(sig.df,Variant == all.segs[si] & Type == "trans")
  trans.sdf <- subset(trans.sdf,!(Gene %in% cis.sdf$Gene))
  
  # account cis -> trans chains
  if (TRUE)
  {
    if (nrow(cis.sdf) > 0 & nrow(trans.sdf) > 0)
    {
      for (ci in 1:nrow(cis.sdf))
      {
        for (ti in 1:nrow(trans.sdf))
        {
          vec <- c(trans.sdf$Gene[ti],all.segs[si],cis.sdf$Gene[ci])
          trio <- rbind(trio,vec)
          reg.type <- c(reg.type,"TMC")
          
          rm(vec)
        }
      }
    }
  }
}

# get indexed trio
trio.g.index <- cbind(match(trio[,1],colnames(G)),match(trio[,2],colnames(L)),match(trio[,3],colnames(G)))
trio.g <- trio
reg.type.g <- reg.type

##### run CIT on trio
library(RCit)

# run analysis with methylation anchored data
citres <- citfun(L = L,G = G,T = G,trios=trio.m.index,maxit=10000)
citres <- data.frame(as.data.frame(trio.m),as.data.frame(citres[,-c(1:3)]),reg.type = reg.type.m)
citres.m <- citres[order(citres$p_cit),]

# run analysis with gene anchored data
citres <- citfun(L = G,G = L,T = G,trios=trio.g.index,maxit=10000)
citres <- data.frame(as.data.frame(trio.g),as.data.frame(citres[,-c(1:3)]),reg.type = reg.type.g)
citres.g <- citres[order(citres$p_cit),]

# merge results
citres <- rbind.data.frame(citres.m,citres.g)

# output CIT results
write.table(citres,file = paste(out.dir,"/CIT_results.methylation.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#write.table(subset(citres,p_cit < 0.2),file = paste(out.dir,"/CIT_results.methylation.signif_0.2.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#write.table(subset(citres,p_cit < 0.1),file = paste(out.dir,"/CIT_results.methylation.signif_0.1.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(subset(citres,p_cit < 0.05),file = paste(out.dir,"/CIT_results.methylation.signif_0.05.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

quit(save = "no",status = 0)
