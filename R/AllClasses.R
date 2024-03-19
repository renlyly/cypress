design_in <- setClass("design_in", slots=c(n_sim="numeric",
                                           n_gene="numeric",
                                           DE_pct="numeric",
                                           ncell_type="numeric",
                                           scenarios="data.frame",
                                           sample_CT_prop = "matrix",
                                           n_scenarios="numeric",
                                           n_strata="numeric",
                                           health_lmean_m="numeric",
                                           health_lmean_d="matrix",
                                           lod_m="numeric",
                                           lod_d="matrix",
                                           health_alpha = "numeric",
                                           case_alpha = "numeric"
)
)


cypress_out <- setClass("cypress_out", slots=c(ct_TDR_bio_smry="data.frame",
                                               TDR_bio_smry="data.frame",
                                               ct_PWR_bio_smry="data.frame",
                                               PWR_bio_smry="data.frame",
                                               PWR_strata_bio_smry="data.frame",
                                          PWR_strata_ct_bio_smry = "data.frame",
                                               ct_FDC_bio_smry="data.frame",
                                               FDC_bio_smry="data.frame"
)
)

est_out<-setClass("est_out",slots = c(health_alpha = "numeric",
                                      case_alpha = "numeric",
                                      health_lmean_m = "numeric",
                                      health_lmean_d = "matrix",
                                      lod_m = "numeric",
                                      lod_d = "matrix",
                                      sample_CT_prop = "matrix",
                                      genename = "character",
                                      samplename = "character",
                                      CTname = "character",
                                      dimensions_Z_hat_ary = "integer"
)
)


