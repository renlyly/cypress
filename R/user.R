###### Below are function for user used


simFromData<- function(INPUTdata = NULL, CT_index = NULL,CT_unk=FALSE,
                        n_sim = 3, n_gene = 30000, DE_pct = 0.05,
                        ss_group_set = c(10, 20, 50, 100),
                        lfc_set = c(0, 0.5, 1, 1.5, 2),
                        lfc_target = 0.5, fdr_thred = 0.1,BPPARAM=bpparam()){

  if (!is(INPUTdata, "SummarizedExperiment")) {
    stop("INPUTdata must be a SummarizedExperiment object.")
  }

  if (!is.null(CT_index) && !is.numeric(CT_index)) {
    stop("CT_index must be numeric.")
  }

  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 2) {
    stop("n_sim must be a single numeric value and at least 2.")
  }

  if (!is.numeric(n_gene) || length(n_gene) != 1 || n_gene < 500) {
    stop("n_gene must be a single numeric value and at least 500.")
  }

  if (!is.numeric(fdr_thred) || length(fdr_thred) != 1 || fdr_thred <= 0) {
    stop("fdr_thred must be a single numeric value greater than 0.")
  }


  if (!is.numeric(lfc_target) || lfc_target < 0) {
    stop("lfc_target must be a positive numeric value.")
  }


  estimate_all <- cypress_est(INPUTdata,CT_index,CT_unk)

  K <- ncol(getcypress(estimate_all, "sample_CT_prop"))
  power_long <-cypressEmbedded(n_sim = n_sim, n_gene = n_gene, DE_pct = DE_pct,
                             ncell_type = K, ss_group_set = ss_group_set,
                             lfc_set = lfc_set,
                             sim_param = estimate_all,
                             lfc_target = lfc_target, fdr_thred = fdr_thred,
                            BPPARAM=BPPARAM)

  return(power_long)
}



simFromParam<- function(n_sim = 3, n_gene = 30000, DE_pct = 0.05,
                         ss_group_set = c(10, 20, 50, 100),
                         lfc_set = c(0, 0.5, 1, 1.5, 2),
                         sim_param="IAD",lfc_target = 0.5, fdr_thred = 0.1,
                        BPPARAM=bpparam()){


  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 2) {
    stop("n_sim must be a single numeric value and at least 2.")
  }
  if (!is.numeric(n_gene) || length(n_gene) != 1 || n_gene < 500) {
    stop("n_gene must be a single numeric value and at least 500.")
  }
  if (!is.numeric(fdr_thred) || length(fdr_thred) != 1 || fdr_thred <= 0) {
    stop("fdr_thred must be a single numeric value and greater than 0.")
  }



  if (sim_param == "IAD") {

    GSE60424_param <- NULL
    data(list = 'quickParaGSE60424', envir = environment())

    K <- ncol(getcypress(GSE60424_param, "health_lmean_d"))
    power_short<- cypressEmbedded(n_sim = n_sim, n_gene = n_gene, DE_pct = DE_pct,
                                ncell_type = K, ss_group_set = ss_group_set,
                                lfc_set = lfc_set,
                                sim_param = GSE60424_param,
                                lfc_target = lfc_target, fdr_thred = fdr_thred,
                                BPPARAM=BPPARAM)
  } else if (sim_param == "IBD") {
    ibd_prop_param <- NULL
    data(list = 'quickParaIBD', envir = environment())
    K <- ncol(getcypress(ibd_prop_param, "health_lmean_d"))
    power_short<- cypressEmbedded(n_sim = n_sim, n_gene = n_gene, DE_pct = DE_pct,
                                ncell_type = K, ss_group_set = ss_group_set,
                                lfc_set = lfc_set,
                                sim_param = ibd_prop_param,
                                lfc_target = lfc_target, fdr_thred = fdr_thred,
                                BPPARAM=BPPARAM)
  }
   else if (sim_param == "ASD") {
      asd_noprop_param <- NULL
      data(list = 'quickParaASD', envir = environment())
      K <- ncol(getcypress(asd_noprop_param, "health_lmean_d"))
      power_short <-cypressEmbedded(n_sim = n_sim, n_gene = n_gene, DE_pct = DE_pct,
                                  ncell_type = K, ss_group_set = ss_group_set,
                                  lfc_set = lfc_set,
                                  sim_param = asd_noprop_param,
                                  lfc_target = lfc_target, fdr_thred = fdr_thred,
                                  BPPARAM=BPPARAM)

  }
  else {
    stop("Invalid parameter value. Please choose 'GSE', 'IBD', or 'ASD'.")
  }
  return(power_short)
}




quickPower<- function(data="IAD"){

  if (data == "IAD") {
    GSE60424Power <- NULL
    data(list = 'quickPowerGSE60424', envir = environment())
    return(GSE60424Power)
  } else if (data == "IBD") {
    ibd_propPower <- NULL
    data(list = 'quickPowerIBD', envir = environment())
    return(ibd_propPower)
  } else if (data == "ASD") {
    asd_propPower <- NULL
    data(list = 'quickPowerASD', envir = environment())
    return(asd_propPower)
  } else {
  stop("Invalid parameter value. Please choose 'IAD', 'IBD', or 'ASD'.")
  }

}
