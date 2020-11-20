
library(limma)
read.data <- function(x,annot.cols = 1:2) 
{
 data.Df <- read.delim(file = x,sep = "\t",header = T)
 data.mat <- as.matrix(data.Df[,-annot.cols])
 rownames(data.mat) <- paste(as.character(data.Df[[2]]),as.character(data.Df[[1]]),sep = "|")
 return(data.mat)
}

run.limma <- function(m,degm)
{
	fit <- lmFit(m,degm)
	fit <- eBayes(fit)
	DEG.table <- topTable(fit,coef=2,number = nrow(m))
	return(DEG.table)
}

count.FD <- function(rho,thresh.vec)
{
	 # combine threshold and correlation values for efficiency
	 n.rho <- length(rho)
	 label.vec <- c(rep(1,length(rho)),rep(0,length(thresh.vec)));# label correlation value and threshold value
	 rho <- c(rho,thresh.vec)# combine correlation and threhold values for efficiency.
	 # sort correlation values
	 i <- order(rho)
	 rho <- rho[i]
	 label.vec <- label.vec[i]
	 
	 # count permuted correlation values above thresholds
	 j <- which(label.vec == 0)
	 output <- cbind(rho[j],1-(n.rho - cumsum(label.vec)[j])/n.rho)
	 colnames(output) <- c("rho.cutoff","FDR")
	 return(output)
} 

limma.DEG <- function(data.mat1,data.mat2,comp.id,
do.permute = TRUE,n.perm = 100)
{
	#data.mat1 <- read.data(data.file1)
	#data.mat2 <- read.data(data.file2)
	comb.mat <- cbind(data.mat1,data.mat2)
	design <- cbind(rep(1,ncol(comb.mat)),c(rep(1,ncol(data.mat1)),rep(0,ncol(data.mat2))));
	colnames(design) <- c("const",comp.id)
	
	# run one-off limma
	cat("-- commence Limma\n")
	DEG.table <- run.limma(m = comb.mat,degm = design)
	DEG.table <- data.frame(processed.probe = rownames(DEG.table),probe = gsub("^(.*)\\|","",rownames(DEG.table)),gene.symbol = gsub("\\|(.*)","",rownames(DEG.table)),DEG.table)
	
	if (do.permute)
	{
	 #min.pwr <- floor(log10(min(DEG.table$P.Value,na.rm = T)))-2
	 #thresh.vec <- c(10^seq(min.pwr,-1,0.005),seq(0.1,1,0.005))
	 
	 thresh.vec <- DEG.table$P.Value
	 cat("-- commence permuting:");cat(n.perm);cat(" times\n");
	 perm.ind <- lapply(1:n.perm,function(i,n) sample(1:n,n),n = ncol(comb.mat))
	 perm.DEG <- lapply(perm.ind,function(i,mat,degg) run.limma(m = mat[,i],degm = degg),mat = comb.mat,degg = design)
	 perm.DEG <- lapply(perm.DEG,function(x) x$P.Value)
	 
	 cat("summarizing stats...\n")
	 z <- rank(thresh.vec);

	 PR <- z/length(thresh.vec);
	 #PR = count.FD(DEG.table$P.Value,thresh.vec);PR = PR[,2]
	 
	 count.out <- lapply(perm.DEG,function(x,vec) count.FD(rho = x,thresh.vec = vec),vec = thresh.vec)
	 FPR = Reduce("+",lapply(count.out,function(x) x[,2]))/n.perm; 
	 FDR = FPR/PR;FDR[which(FPR == 0)] <- 0;FDR[which(FDR > 1)] <- 1;
	 FDR.table <- data.frame(threshold = thresh.vec,FPR = FPR,PR = PR,FDR = FDR)
	 DEG.table <- data.frame(DEG.table,FDR.perm = FDR)
	 rm(perm.DEG,PR,count.out,FPR,FDR)
	}
	
	return(DEG.table)
}
