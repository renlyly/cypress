
plotFDC<-function(simulation_results,sample_size=10){
  ss<-unique(simulation_results@FDC_bio_smry$ss)
  lfc_mean<-unique(simulation_results@FDC_bio_smry$lfc)

  if(!(sample_size %in% ss)) stop("Sample size should be one of your design sample size set.")
  if(length(lfc_mean)<2) stop("The length of lfc_set input should be greater than 2")

  ###FDC result for ct_smry
  ct_n<-ncol(simulation_results@ct_FDC_bio_smry)-3
  ct_smry_tmp<-simulation_results@ct_FDC_bio_smry
  ct_smry<-ct_smry_tmp[which(ct_smry_tmp$ss==sample_size),seq_len(ct_n)]

  ###FDC result for ss_smry
  ss_FDC<-simulation_results@FDC_bio_smry$FDC
  ss_smry<-matrix(ss_FDC,nrow =length(lfc_mean),ncol=length(ss),byrow = TRUE)
  rownames(ss_smry)<-lfc_mean
  colnames(ss_smry)<-ss


  col<-RColorBrewer::brewer.pal(6,"Set1")
  pch<-c(16,17,18,2,8,15)
  ylab <-c('FDC')
  par(mfrow=c(1,2))

  ##1
  matplot(x=lfc_mean,ct_smry, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0, max(ct_smry,na.rm = TRUE)),xlab='Effect size',
          ylab = ylab,xaxt="n")
  axis(1, at=lfc_mean,labels=lfc_mean, las=0)
  legend("topright",c(paste0('Celltype', seq_len(ct_n))),
         col=col, pch=pch,cex=0.7,ncol = 2,
         y.intersp=0.6,bty = 'n',x.intersp = 0.5,text.width = 0.5,horiz = FALSE)

  ##2
  matplot(x=lfc_mean,ss_smry, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,max(ss_smry,na.rm = TRUE)),xlab='Effect size',
          ylab = ylab,xaxt="n")
  axis(1, at=lfc_mean,labels=lfc_mean, las=0)
  legend("bottomleft",c(paste0('N=',ss)),col=col, pch=pch,pt.cex=1,cex=0.7,
         y.intersp=0.8,bty = 'n',x.intersp = 0.6,text.width = 0.3,horiz = FALSE)
}


