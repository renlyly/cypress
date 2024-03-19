
plotPower<-function(simulation_results,effect.size=1,sample_size=10){
  if (!is(simulation_results, "cypress_out")) {
    stop("simulation_results must be an S4 object of class 'cypress_out'.")
  }
  if (!is.numeric(sample_size) || sample_size <= 5) {
    stop("sample_size must be a numeric value greater than 5.")
  }
  if (!is.numeric(effect.size) || effect.size <= 0) {
    stop("sample_size must be a numeric value greater than 0.")
  }

  ss <- unique(getcypress(simulation_results, "FDC_bio_smry")$ss)
  lfc_mean <- unique(getcypress(simulation_results, "FDC_bio_smry")$lfc)


  if(!(sample_size %in% ss)) stop("Sample size should be one of your design sample size set.")
  if(!(effect.size %in% lfc_mean)) stop("Effect size should be one of your design effect size set.")
  if(length(lfc_mean)<2) stop("The length of lfc_set input should be greater than 2")

  ## power result for ss_smry (effect_size x sample)
  ss_pwr <- getcypress(simulation_results, "PWR_bio_smry")$PWR

  ss_smry<-matrix(ss_pwr,nrow =length(lfc_mean),ncol=length(ss),byrow = TRUE)
  rownames(ss_smry)<-lfc_mean
  colnames(ss_smry)<-ss

  ## power result for ct_lfc_smry (effect_size x celltype)

  ct_n <- ncol(getcypress(simulation_results, "ct_PWR_bio_smry")) - 3
  ct_smry_tmp <- getcypress(simulation_results, "ct_PWR_bio_smry")

  ct_smry<-ct_smry_tmp[which(ct_smry_tmp$ss==sample_size),seq_len(ct_n)]

  ## power result for ct_ss_smry (sample_size x celltype)
  ct_ss_smry<-ct_smry_tmp[which(ct_smry_tmp$lfc==effect.size),seq_len(ct_n)]

  ## power result for strata_smry (strata x sample_size | strata x effect_size)
  strata<- seq_len(9)

  strata_smry_tmp <- getcypress(simulation_results, "PWR_strata_bio_smry")

strata_smry_ss<-t(strata_smry_tmp[which(strata_smry_tmp$lfc==effect.size),strata])
strata_smry_eff<-t(strata_smry_tmp[which(strata_smry_tmp$ss==sample_size),strata])

  ## power result for strata_ct_smry
strata_smry_ct_tmp <- getcypress(simulation_results, "PWR_strata_ct_bio_smry")

strata_smry_ct<-t(strata_smry_ct_tmp[which(strata_smry_ct_tmp$lfc==effect.size &
                                      strata_smry_ct_tmp$ss==sample_size),strata])


  ##. draw power figures
  par(mfrow=c(2,3))
  ylab <- c('Power')
  col<-RColorBrewer::brewer.pal(9,"Set1")
  pch<-c(16,17,18,2,8,15)
  ## 1
  matplot(lfc_mean,ss_smry,type='b',pch=pch,col=col,lwd=1.5,ylim=c(0,1),
          xlab='Effect size',ylab=ylab,xaxt="n")
  axis(1, at=lfc_mean,labels=lfc_mean, las=0)
  legend("topleft",c(paste0('N=',ss)),col=col, pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.5,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)
  ## 2
  matplot(lfc_mean,ct_smry,type='b',pch=pch,col=col,lwd=1.5,ylim=c(0,1),
          xlab='Effect size',ylab=ylab,xaxt="n")
  axis(1, at=lfc_mean,labels=lfc_mean, las=0)
  legend("topleft",c(paste0('Celltype',seq_len(ct_n))),col=col,
         pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.5,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)

  ## 3
  matplot(ss,ct_ss_smry,type="b",cex=1,pch=pch,col=col,lwd=1.5,ylim=c(0,1),
          xlab='Sample size',ylab=ylab,xaxt="n")
  axis(1, at=ss,labels=ss, las=0)
  legend("topleft",c(paste0('Celltype',seq_len(ct_n))),
         col=col, pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.5,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)

  ## 4
  matplot(x=strata,y=strata_smry_ct, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Strata', ylab = ylab,xaxt="n")
  axis(1, at=seq_len(9),labels=seq_len(9), las=0)
  legend("topleft",c(paste0('Celltype',seq_len(ct_n))),col=col, pch=pch,cex=0.8,
         y.intersp=0.5,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)

  ## 5
  matplot(x=strata,y=strata_smry_ss, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Strata', ylab = ylab,xaxt="n")
  axis(1, at=seq_len(9),labels=seq_len(9), las=0)
  legend("topleft",c(paste0('N=',ss)),col=col, pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.5,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)

  ## 6
  matplot(x=strata,strata_smry_eff, type="b",pch=pch,cex=1,lwd=1.5,col=col,
          ylim=c(0,1),xlab='Strata', ylab = ylab,xaxt="n")
  axis(1, at=seq_len(9),labels=seq_len(9), las=0)
  legend("topleft",c(paste0('Effect size=',lfc_mean)),
         col=col, pch=pch,pt.cex=1,cex=0.8,
         y.intersp=0.5,bty = 'n',x.intersp = 0.5,text.width = 0.3,horiz = FALSE)

}

