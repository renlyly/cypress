#########################################################
### cypress_wrapper returns a S4 object with 7 elements.
###
# 'exp_design': a S4 object returned by 'exp_design' function.
# 'lfc_target': target effect size, should be greater than or equal to 0.
# The absolute LFC lower than this values will be treated as None-DEGs.
# Default to 0.5
# 'fdr_thred': adjusted p value threshold, should not below 0 or exceed 1.
#
# Returns:
# 'ct_TDR_bio_smry': Cell type specific target TDR for all scenarios.
# Columns of 'ct', 'ss', 'lfc', and 'lfc_sd' respectively indicates cell type
# index,
# sample size, effect size, and effect size standard deviation.
#
# 'TDR_bio_smry': Target TDR for all scenarios averaged across cell types.
#
# 'ct_PWR_bio_smry': Cell type specific target power for all scenarios. Columns
# of 'ss', 'lfc', and  'lfc_sd' respectively indicates sample size, effect size,
# and effect size standard deviation.
#
# 'PWR_bio_smry': Target power for all scenarios averaged across cell types.
#
# 'PWR_strata_bio_smry': Target power for each strata under simulation scenarios
# averaged across cell types.
#
# 'ct_FDC_bio_smry': Cell type specific target FDC for all scenarios.
#
# 'FDC_bio_smry': Target FDC for all scenarios averaged across cell types.



cypress_wrapper <- function(exp_design,
                            lfc_target, fdr_thred,BPPARAM=bpparam()){

  # simulation setting
  if(!(lfc_target >= 0)) stop("Target LFC should be non-negative")
  if(fdr_thred < 0 || fdr_thred > 1) stop("FDR threshold should be between 0 and 1")

  scenarios <- getcypress(exp_design, "scenarios")
  n_scenarios <- getcypress(exp_design, "n_scenarios")
  ncell_type <- getcypress(exp_design, "ncell_type")
  DE_pct <- getcypress(exp_design, "DE_pct")
  n_sim <- getcypress(exp_design, "n_sim")
  n_gene <- getcypress(exp_design, "n_gene")
  n_strata <- getcypress(exp_design, "n_strata")
  # distribution parameter
  health_lmean_m <- getcypress(exp_design, "health_lmean_m")
  health_lmean_d <- getcypress(exp_design, "health_lmean_d")
  lod_m <- getcypress(exp_design, "lod_m")
  lod_d <- getcypress(exp_design, "lod_d")
  # proportion parameter
  health_alpha <- getcypress(exp_design, "health_alpha")
  case_alpha <- getcypress(exp_design, "case_alpha")

  ct_TDR_bio_smry <- NULL
  TDR_bio_smry <- NULL
  ct_PWR_bio_smry <- NULL
  PWR_bio_smry <- NULL
  PWR_strata_bio_smry <- NULL
  PWR_strata_ct_bio_smry <- NULL
  ct_FDC_bio_smry <- NULL
  FDC_bio_smry <- NULL
  lfc_sd <- 0.5
  top_g <- seq(50, 600, by = 100)
  for(i in seq_len(n_scenarios)){
    lfc_mean <- scenarios[i,'lfc']
    nsample_each_group <- scenarios[i,'ss_group']
    message("LFC: ",lfc_mean," Sample Size: ",nsample_each_group)
    add_columns <- data.frame(
      ct = seq_len(ncell_type),
      ss = nsample_each_group,
      lfc = lfc_mean,
      lfc_sd = lfc_sd
    )
    # Implement simulation
   csRNA_results <- BiocParallel::bplapply(seq_len(n_sim),csRNA_seq_sim,
                              n_gene,DE_pct,ncell_type,lfc_mean,lfc_sd,
                              health_lmean_m,health_lmean_d,lod_m,lod_d,
                              nsample_each_group,lfc_target,fdr_thred,top_g,
                              health_alpha,case_alpha,n_strata,BPPARAM=BPPARAM)
    # cell type specific target TDR
    ct_TDR_bio <- lapply(csRNA_results, function(x) x$ct_TDR_bio)
    ct_TDR_bio <- abind(ct_TDR_bio, along = 3)
    ct_TDR_bio <- apply(ct_TDR_bio,seq_len(2),mean,na.rm=TRUE)
    ct_TDR_bio <- cbind(ct_TDR_bio,add_columns)
    ct_TDR_bio_smry <- rbind(ct_TDR_bio_smry,ct_TDR_bio)
    # target TDR
    TDR_bio <- lapply(csRNA_results, function(x) x$TDR_bio)
    TDR_bio <- abind(TDR_bio, along = 0)
    TDR_bio <- apply(TDR_bio,2,mean,na.rm=TRUE)
    TDR_bio <- cbind(t(TDR_bio),
                data.frame(ss=nsample_each_group,lfc=lfc_mean,lfc_sd=lfc_sd))
    TDR_bio_smry <- rbind(TDR_bio_smry,TDR_bio)
    # cell type specific target power
    ct_PWR_bio <- lapply(csRNA_results, function(x) x$ct_PWR_bio)
    ct_PWR_bio <- abind(ct_PWR_bio, along = 0)
    ct_PWR_bio <- apply(ct_PWR_bio,2,mean,na.rm=TRUE)
    ct_PWR_bio <- cbind(t(ct_PWR_bio),
                   data.frame(ss = nsample_each_group,lfc=lfc_mean,lfc_sd=lfc_sd))
    ct_PWR_bio_smry <- rbind(ct_PWR_bio_smry,ct_PWR_bio)
    # target power
    PWR_bio <- lapply(csRNA_results, function(x) x$PWR_bio)
    PWR_bio <- abind(PWR_bio, along = 0)
    PWR_bio <- mean(PWR_bio,na.rm = TRUE)
    PWR_bio <- data.frame(t(c( PWR_bio,nsample_each_group,lfc_mean,lfc_sd )))
    colnames(PWR_bio) <- c("PWR","ss","lfc","lfc_sd")
    PWR_bio_smry <- rbind(PWR_bio_smry,PWR_bio)
    # target power by gene expression stratification
    PWR_strata_bio <- lapply(csRNA_results, function(x) x$PWR_strata_bio)
    PWR_strata_bio <- abind(PWR_strata_bio, along = 1)
    PWR_strata_bio <- apply(PWR_strata_bio,2,mean,na.rm=TRUE)
    PWR_strata_bio <- cbind(t(PWR_strata_bio),
                   data.frame(ss=nsample_each_group,lfc=lfc_mean,lfc_sd=lfc_sd))
    PWR_strata_bio_smry <- rbind(PWR_strata_bio_smry,PWR_strata_bio)
    # cell type specific target power by gene expression stratification
    PWR_strata_ct_bio <- lapply(csRNA_results, function(x) x$PWR_strata_ct_bio)
    PWR_strata_ct_bio <- abind(PWR_strata_ct_bio, along = 3)
    PWR_strata_ct_bio <- apply(PWR_strata_ct_bio,seq_len(2),mean,na.rm=TRUE)
    PWR_strata_ct_bio <- cbind(PWR_strata_ct_bio,
                               add_columns)
    PWR_strata_ct_bio_smry <- rbind(PWR_strata_ct_bio_smry, PWR_strata_ct_bio)
    # cell type specific target FDC
    ct_FDC_bio <- lapply(csRNA_results, function(x) x$ct_FDC_bio)
    ct_FDC_bio <- abind(ct_FDC_bio, along = 0)
    ct_FDC_bio <- apply(ct_FDC_bio,2,mean,na.rm=TRUE)
    ct_FDC_bio <- cbind(t(ct_FDC_bio),
                   data.frame(ss=nsample_each_group,lfc=lfc_mean,lfc_sd=lfc_sd))
    ct_FDC_bio_smry <- rbind(ct_FDC_bio_smry,ct_FDC_bio)
    # target FDC
    FDC_bio <- lapply(csRNA_results, function(x) x$FDC_bio)
    FDC_bio <- abind(FDC_bio, along = 1)
    FDC_bio <- mean(FDC_bio,na.rm = TRUE)
    FDC_bio <- data.frame(t(c(FDC_bio,nsample_each_group,lfc_mean,lfc_sd)))
    colnames(FDC_bio) <- c("FDC","ss","lfc","lfc_sd")
    FDC_bio_smry <- rbind(FDC_bio_smry,FDC_bio)
  }

  cypress_smry <- cypress_out(ct_TDR_bio_smry = ct_TDR_bio_smry,
                              TDR_bio_smry = TDR_bio_smry,
                              ct_PWR_bio_smry = ct_PWR_bio_smry,
                              PWR_bio_smry = PWR_bio_smry,
                              PWR_strata_bio_smry = PWR_strata_bio_smry,
                              PWR_strata_ct_bio_smry = PWR_strata_ct_bio_smry,
                              ct_FDC_bio_smry = ct_FDC_bio_smry,
                              FDC_bio_smry = FDC_bio_smry)
  message("Complete simulation for CYPRESS.")
  return(cypress_smry)
}
