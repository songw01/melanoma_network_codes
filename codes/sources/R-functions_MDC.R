###############################################################################################
#
#                    Modular Differential Connectivity (MDC)
#
#    Function: Compute ratio of the mean intra-module connectivity in one network to that of 
#                the same set of genes in another network
#
#    Author: Bin Zhang
#            binzhang.ucla@gmail.com or bin.zhang@mssm.edu
#
#    Copyright: @MNML at Moutn Sinai
#
#    Reference: Zhang B et al. Cell 153(3):707-720. PMID: 23622250 
#
#    Please keep the code to your group as the code will be published along with a method paper. 
#
#    A weighted interaction network analysis (WINA) and a R package for computing the modular 
#     differential connectivity (MDC) will be on the website:
#       http://research.mssm.edu/multiscalenetwork/Resources.html
#
#
#    Input:
#         1) inputfnameB       -- mRNA data (tab-delimited text file) for one disease state B
#         2) inputfnameBModule -- module assignment file as output from WINA or WGCNA
#         3) inputfname        -- mRNA data (tab-delimited text file) for another disease state A
#         4) shortnames        -- shortnames for B and A, respectively. Order is important
#         5) headCol           -- # of gene information columns in inputfname
#         6) headCol2          -- # of gene information columns in inputfnameB
#         7) corrlpower        -- the exponent of the power function, determined by WINA or WGCNA
#
#   Output:
#         1) "*_wFDR.xls"      -- MDC with false discovery rates
#         2) "*.png"           -- MDC plot
#
#
#   Last Modified: Bin Zhang, 09/14/2013 
#
################################################################################################

# to split "abc|123", use sep="\\|", "abc.123" use "\\."
splitString =function(mystring, separator="; "){
  splitted = NULL
  for (each in mystring){
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a)
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
}


#get the filename without extension
#
getFileExtension=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    
    if( length(splitted) >1){
      return (splitted[length(splitted)])
    } else{
      return ("")
    }
}

#get the filename without extension
getFileName=function(fullfname){
    ext=getFileExtension(fullfname)
    if(ext ==""){
       return (fullfname)
    }
    extd = paste(".", ext, sep="")
    splitted=splitString(fullfname, extd)

    splitted[1]
}

#get the filename without extension
getFileNames=function(fullfnames){

  final = NULL
  for(ef in fullfnames) {
     fn = getFileName(ef)
     final = c(final, fn)
  }
  return (final)
}

#get the filename without extension and path information
getFileNameNopath=function(fullfname){
   res = NULL
   for(each in fullfname) {
    myfilename = getFileName(each)
    splitted=unlist( strsplit(myfilename, "/") )
     res= c(res, splitted[length(splitted) ])
   }
   return (res)
}

# assume that find indices of the common components in the two sets
#
getMatchedIndex2way=function(cvector, dvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(dvect))
  subIdex   = cbind(dvect, index2)

  merged    = merge(orgIdx, subIdex, by.x=1, by.y=1, all=F)
  merged    = as.matrix(merged)
  
  outIndex  = cbind(as.integer(merged[,2]), as.integer(merged[,3]) )
  mo = order(outIndex[,1])

  outIndex = outIndex[mo,]

  return (outIndex)
}

permuteVect = function(myvect){
  pidx = sample(c(1:length(myvect)), length(myvect), replace=F)
  return(myvect[pidx])
}

# a matrix with three cols: MDC random_MDC_mean random_MDC_sd
#
MDC_FDR = function(mdcvect){
  if(mdcvect[1] < 1) {
     fdr = sum(mdcvect[1]>mdcvect[-1])/length(mdcvect[-1])
  } else {
     fdr = sum(mdcvect[1]<mdcvect[-1])/length(mdcvect[-1])
  }
  return (fdr)
}

# a matrix with three cols: MDC random_MDC_mean random_MDC_sd
#
MDC_FDR_by_normal_distr = function(mdcvect){
  if(mdcvect[1] < 1) {
     fdr = pnorm(mdcvect[1], mean=mdcvect[2], sd=mdcvect[3], lower.tail=TRUE)
  } else {
     fdr = pnorm(mdcvect[1], mean=mdcvect[2], sd=mdcvect[3], lower.tail=FALSE)
  }
  return (fdr)
}

# xaxisLabel_rotate = 45 # slant x-labels
#
barplot_by_oneGroup = function(datMatrix, maskMatrix=NULL, 
                  rowcol=c("red", "green", "blue", "grey"), vlines=NULL,
                  title="", ylab="", imgfile=NULL, imgwid=800, imghei=400, 
                  mcex=1.6, xlab_cex=1, xaxisLabel_rotate=NULL, usr_delta=2, 
                  legendcols=1, legendxDelta=0, legendCex=0.8, 
                  ipointsize=12, iunits="px", ires=72, icompression="lzw",
                  show_value=T, showLegend=T, show_xaxisLabel=T, show_mask=T,
                  xmargin=4, ymargin=4, zmargin=1, xlabel_orientation=0, ValueXdelta=0, 
                  ValueYdelta=NULL, valuecex=0.8, yMAX=NULL, yMIN=0) {

   xmatrixMj   = datMatrix
   xmaskMatrix = maskMatrix
 
   if(!is.null(maskMatrix)){
      bar_y = paste(datMatrix[,1], maskMatrix[,1] , sep="")
   }

   tab_names  = rownames(xmatrixMj)
   no.groupsX = dim(xmatrixMj)[1]; no.tabs=dim(xmatrixMj)[2]

      interv = 1; space4grp = interv+1
      xcoord =  0.5+ seq(1, space4grp*no.groupsX, space4grp)
      
   gene_max = max(xmatrixMj)
   which_min= as.integer(which.min(gene_max))

   if( is.null(yMAX) ) {
      ymax = (as.integer(max(xmatrixMj)/200)+1)*200; #max(ymeans);
      ymax = (as.integer(max(xmatrixMj)/20)+1)*20; #max(ymeans);
   } else{
      ymax = yMAX
   }

   if(is.null(yMIN) ) {ymin=min(xmatrixMj)
   } else {ymin=yMIN}

   mycolor = rowcol[c(1:no.groupsX)]

   if(!is.null(imgfile)) {
     openImgDev(imgfile, iwidth =imgwid, iheight = imghei, ipointsize=ipointsize, iunits=iunits, ires=ires, icompression=icompression)
     par(mfrow=c(1,1), mar=c(xmargin, ymargin, zmargin, 0) + 0.1, cex=mcex,las=1)#, srt=90)
   }

   barplot(datMatrix[,1], beside = TRUE, space=interv,
        col = mycolor, axisnames=F, ylab=ylab, ylim = c(ymin, ymax),
        legend =F)#shortnames) #legend =F

   if(title!="") {
       title(main = title, font.main = 1)
   }

   #err.bp(t(ymeans), t(ysds), two.side=F) 
   if(show_xaxisLabel) {
       if(is.null(xaxisLabel_rotate) ) {
          axis(1, at =xcoord, labels =rownames(xmatrixMj), las=xlabel_orientation, 
               col = "black", col.axis="black", lwd = 0, tick = T, line =1)
       } else{
          axis(1, at =xcoord, labels =rep("", nrow(datMatrix)), las=xlabel_orientation, 
               col = "black", col.axis="black", lwd = 0, tick = T, line =1)
          text(xcoord, par("usr")[3] - usr_delta, srt = xaxisLabel_rotate, adj = 1, labels = rownames(xmatrixMj), 
          xpd = TRUE, cex=xlab_cex)
          par(srt=0) 
       }
   }

   ValueYdelta2 = ValueYdelta
   if(is.null(ValueYdelta)){
      ValueYdelta2 = ymax/20
   }

   if(show_value) {
      text(x=xcoord, y=datMatrix[,1]+ValueYdelta2, srt=0, labels=bar_y, cex=valuecex, col="black")
   }

   if(show_mask) {
      text(x=xcoord, y=datMatrix[,1]+ValueYdelta2, srt=0, labels=maskMatrix[,1], cex=valuecex, col="black")
   }


   if(!is.null(vlines)){
        for(xx in vlines){abline(h=xx,col="gray", lty=3)}
   }

   if(!is.null(imgfile)) {
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, srt=0)
     dev.off()
   }

}

closest_integer_bound = function(val){
  logv = log10(val)
  ii   = as.integer(logv)
  idec = val/(10^ii) # 2.50
  iremaind = idec - as.integer(idec)

  ix = ifelse(iremaind<0.5, 0.5, 1)
  return ( (as.integer(idec)+ix)*(10^ii) )
}

#iunits: "px", "cm", "mm"
#
# Letter Size 8.5 x 11 inches
# images widths: 8.3cm ~ 3.26in, 12.35cm ~ 4.86in, 17.35cm ~ 6.83in
#        heights: max 23.35cm ~ 9.19in
#
openImgDev=function(imgname, iwidth = 1024, iheight=1024, ipointsize=12, iunits="px", ires=72, icompression="lzw")
{
  imgtype = getFileExtension(imgname)
  
  if (imgtype=="ps"){
     postscript(imgname,width=iwidth, height=iheight, pointsize=ipointsize)
  }else if (imgtype=="png"){
     png(imgname, width=iwidth, height=iheight, pointsize=ipointsize, units=iunits, res=ires)
  }else if ( (imgtype=="jpg") | (imgtype=="jpeg") ){
     jpeg(imgname, width=iwidth, height=iheight, pointsize=ipointsize,units=iunits, res=ires,quality =100)
  }else if ( (imgtype=="tif")|(imgtype=="tiff") ){
     tiff(imgname, width=iwidth, height=iheight, pointsize=ipointsize,units=iunits, res=ires,compression=icompression)
  }else if ( (imgtype=="pdf") | (imgtype=="PDF") ){
     pdf(imgname)
     return
  }else{
     png(imgname, width=iwidth, height=iheight, pointsize=ipointsize, units=iunits, res=ires)
  }
  trellis.device(new = FALSE, col = TRUE) 
}

collect_garbage=function(){
    while (gc()[2,4] != gc()[2,4]){}
}

SetTextContrastColor <- function(color)
{
  ifelse( mean(col2rgb(color)) > 127, "black", "white")
}

repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}
