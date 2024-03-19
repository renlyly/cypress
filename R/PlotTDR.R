plotTDR<-function(simulation_results,effect.size=1,sample_size=10){
  if (!is(simulation_results, "cypress_out")) {
    stop("simulation_results must be an S4 object of class 'cypress_out'.")
  }
  if (!is.numeric(sample_size) || sample_size <= 5) {
    stop("sample_size must be a numeric value greater than 5.")
  }

  ss <- unique(getcypress(simulation_results, "FDC_bio_smry")$ss)
  lfc_mean <- unique(getcypress(simulation_results, "FDC_bio_smry")$lfc)

  if(!(sample_size %in% ss)) stop("Sample size should be one of your design sample size set.")
  if(!(effect.size %in% lfc_mean)) stop("Effect size should be one of your design effect size set.")
  if(length(lfc_mean)<2) stop("The length of lfc_set input should be greater than 2")
  top_rank<-seq(50,600,100)

  ct_n <- length(unique(getcypress(simulation_results, "ct_TDR_bio_smry")$ct))
  ct_smry_tmp <- getcypress(simulation_results, "ct_TDR_bio_smry")

  ct_smry<-ct_smry_tmp[which(ct_smry_tmp$ss==sample_size&ct_smry_tmp$lfc==effect.size ),seq_along(top_rank)]
  ct_smry<-t(ct_smry)

  TDR_smry <- getcypress(simulation_results, "TDR_bio_smry")

  ef_smry_tmp<-TDR_smry[which(TDR_smry$ss==sample_size),seq_along(top_rank)]
  ef_smry<-t(ef_smry_tmp)
  ss_smry_tmp<-TDR_smry[which(TDR_smry$lfc==effect.size),seq_along(top_rank)]
  ss_smry<-t(ss_smry_tmp)
  Top_TDR<-TDR_smry[,"350"]
  ef_ss_smry<-matrix(Top_TDR,nrow =length(lfc_mean),ncol=length(ss),byrow = TRUE)
  par(mfrow=c(2,2))
  ylab <- c('TDR')
  col<-RColorBrewer::brewer.pal(6,"Set1")
  pch <- c(16,17,18,2,8,15)
  matplot(x=top_rank,ct_smry, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Top rank genes', ylab = ylab,xaxt="n")
  axis(1, at=top_rank,labels=top_rank, las=0)
  legend("bottomleft",c(paste0('Celltype',seq_len(ct_n))),col=col,
         pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.6,bty = 'n',x.intersp = 0.5,horiz = FALSE,ncol=2)
  matplot(x=top_rank,ef_smry, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Top rank genes', ylab = ylab,xaxt="n")
  axis(1, at=top_rank,labels=top_rank, las=0)
  legend("topright",c(paste0('Effect size=',lfc_mean)),
         col=col, pch=pch,pt.cex=1,cex=0.8,y.intersp=0.5,bty = 'n',
         x.intersp = 0.5,horiz = FALSE)
  matplot(x=top_rank,ss_smry, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Top rank genes', ylab = ylab,xaxt="n")
  axis(1, at=top_rank,labels=top_rank, las=0)
  legend("bottomleft",c(paste0('N=',ss)),
         col=col, pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.6,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)
  matplot(x=lfc_mean,ef_ss_smry, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Effect size', ylab = ylab,xaxt="n")
  axis(1, at=lfc_mean,labels=lfc_mean, las=0)
  legend("topleft",c(paste0('N=',ss)),
         col=col, pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.6,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)
}



