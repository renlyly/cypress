## Functions below are used to implement TOAST.
## The first function is to run TOAST model.
## The following two functions are to organize TOAST outputs, which are lists.
## The last function was used as lapply input for the second function, which
## serves as summarization purposes.
run_TOAST <- function(nsample_each_group,est_CT_prop,RNAseq_final_count) {
  ncell_type <- ncol(est_CT_prop)
  sim_design <- as.data.frame(factor(rep(c(1,2),nsample_each_group)))
  colnames(sim_design) <- "disease"
  sim_Design_out <- makeDesign(sim_design,est_CT_prop)
  sim_fitted_model_strata <- fitModel(sim_Design_out,
                                      as.matrix(RNAseq_final_count))
  sim_res_TOAST_strata <- csTest(sim_fitted_model_strata, coef = "disease",
                                 cell_type = NULL, verbose = FALSE, sort = TRUE)
  return(sim_res_TOAST_strata)
}
summary_TOAST <- function(sim_res_TOAST_strata,ncell_type){
  TOAST_out_name <- names(sim_res_TOAST_strata[seq_len(ncell_type)])
  TOAST_out <- lapply(seq_len(ncell_type), pull_TOAST,sim_res_TOAST_strata)
  TOAST_out <- do.call(cbind, TOAST_out)
  TOAST_out <- TOAST_out[c("gene_id",TOAST_out_name)]
  return(TOAST_out)
}
pull_TOAST <- function(i,sim_res_TOAST_strata){
  TOAST_each_cell <- sim_res_TOAST_strata[[i]]
  TOAST_each_cell$gene_id <- as.numeric(rownames(TOAST_each_cell))
  TOAST_each_cell <- TOAST_each_cell[c("gene_id",'fdr')]
  TOAST_each_cell <- TOAST_each_cell[order(TOAST_each_cell$gene_id),]
  names(TOAST_each_cell)[2] <- names(sim_res_TOAST_strata)[i]
  return(TOAST_each_cell)
}
