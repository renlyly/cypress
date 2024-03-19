design_in <- setClass("design_in", slots=c(n_sim="numeric",
                                           n_gene="numeric",
                                           DE_pct="numeric",
                                           ncell_type="numeric",
                                           scenarios="data.frame",
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




setValidity("cypress_out", function(object) {
  msg <- NULL
  slotNames <- c("ct_TDR_bio_smry", "TDR_bio_smry", "ct_PWR_bio_smry", "PWR_bio_smry",
                 "PWR_strata_bio_smry", "PWR_strata_ct_bio_smry", "ct_FDC_bio_smry", "FDC_bio_smry")

  for (slotName in slotNames) {
    if (!is.data.frame(slot(object, slotName))) {
      msg <- c(msg, paste("'", slotName, "' must be a data.frame.", sep = ""))
    }
  }

  if (length(msg) == 0) TRUE else msg
})



setValidity("est_out", function(object) {
  msg <- NULL
  if (!is.numeric(object@health_alpha)) msg <- c(msg, "'health_alpha' must be numeric.")
  if (!is.numeric(object@case_alpha)) msg <- c(msg, "'case_alpha' must be numeric.")
  if (!is.numeric(object@health_lmean_m)) msg <- c(msg, "'health_lmean_m' must be numeric.")
  if (!is.matrix(object@health_lmean_d)) msg <- c(msg, "'health_lmean_d' must be a matrix.")
  if (!is.numeric(object@lod_m)) msg <- c(msg, "'lod_m' must be numeric.")
  if (!is.matrix(object@lod_d)) msg <- c(msg, "'lod_d' must be a matrix.")
  if (!is.matrix(object@sample_CT_prop)) msg <- c(msg, "'sample_CT_prop' must be a matrix.")
  if (!is.character(object@genename)) msg <- c(msg, "'genename' must be character.")
  if (!is.character(object@samplename)) msg <- c(msg, "'samplename' must be character.")
  if (!is.character(object@CTname)) msg <- c(msg, "'CTname' must be character.")
  if (!is.integer(object@dimensions_Z_hat_ary)) msg <- c(msg, "'dimensions_Z_hat_ary' must be integer.")
})



