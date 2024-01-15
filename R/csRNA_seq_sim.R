
csRNA_seq_sim <- function(i,n_gene,DE_pct,ncell_type,lfc_mean,lfc_sd,
                          health_lmean_m,health_lmean_d,lod_m,lod_d,
                          nsample_each_group,lfc_target,fdr_thred,top_g,
                          health_alpha,case_alpha,n_strata){
  lfc_ct_mtx <- lfc_simulator(n_gene, DE_pct, ncell_type,
                              lfc_mean, lfc_sd)
  cs_ref_panel_all <- cs_ref_panel_simulator(ncell_type, n_gene,
                                             health_lmean_m, health_lmean_d,
                                             lod_m, lod_d,lfc_ct_mtx)
  cellT_expr_healthy <- cs_ref_panel_all[[1]]
  colnames(cellT_expr_healthy) <- c(paste0("Celltype",seq_len(ncell_type)))
  cellT_expr_case <- cs_ref_panel_all[[2]]
  colnames(cellT_expr_case) <- c(paste0("Celltype",seq_len(ncell_type)))
  # Here simulate RNAseq count data and true proportion.
  sample_results <- lapply(seq_len(nsample_each_group) , simulate_sample,
                           cellT_expr_healthy,cellT_expr_case,
                           health_alpha,case_alpha,
                           ncell_type,n_gene)
  # Extract 'cell_prop' matrices and 'cellT_expr_count' matrices
  cell_prop_list <- lapply(sample_results, function(sample) sample$cell_prop)
  cellT_expr_count_list <- lapply(sample_results,
                                  function(sample) sample$cellT_expr_count)
  # Row-bind 'cell_prop' matrices and column-bind 'cellT_expr_count' matrices
  sample_CT_prop <- do.call(rbind, cell_prop_list)
  RNAseq_final_count <- do.call(cbind, cellT_expr_count_list)
  # Deconvolution
  refinx <- findRefinx(cellT_expr_healthy, nmarker = 1000, sortBy = "cv")
  ctrl_est_prop <- perform_deconvolution(cellT_expr_healthy,"control",
                                         RNAseq_final_count,refinx,ncell_type)
  case_est_prop <- perform_deconvolution(cellT_expr_case,"case",
                                         RNAseq_final_count,refinx,ncell_type)
  est_CT_prop <- rbind(ctrl_est_prop, case_est_prop)
  est_CT_prop <- est_CT_prop[rownames(sample_CT_prop),]
  gene_CT_DE_connect <- gt_organizer(RNAseq_final_count, lfc_ct_mtx, ncell_type)
  # ##TOAST Implementation.
  sim_res_TOAST_strata <- run_TOAST(nsample_each_group,
                                    est_CT_prop,RNAseq_final_count)
  TOAST_out <- summary_TOAST(sim_res_TOAST_strata,ncell_type)
  gene_CT_FDR <- merge(TOAST_out, gene_CT_DE_connect, by = "gene_id")
  # Metric Implementations
  gene_CT_FDR <- cbind(gene_CT_FDR,
                       DE_CT_bio = ifelse(abs(gene_CT_FDR$LFC) < lfc_target, 0,
                                          gene_CT_FDR$DE_CT))
  # Target TDR
  ct_tdr_bio <- TDR_pull(gene_CT_FDR,top_g,ncell_type)
  tdr_bio <- apply(ct_tdr_bio,2,mean,na.rm=TRUE)
  # Target Power
  ct_pwr_bio <- POWER_pull(gene_CT_FDR,fdr_thred,ncell_type)
  ct_pwr <- ct_pwr_bio[[1]];  ct_pwr_m <- ct_pwr_bio[[2]]
  # Target Power by strata.
  POWER_strata <- lapply(seq_len(n_strata),POWER_strata_pull,gene_CT_FDR,
                         fdr_thred,ncell_type)
  POWER_strata_bio <- lapply(POWER_strata, function(x) x[[2]])
  POWER_strata_bio <- do.call(cbind, POWER_strata_bio)
  colnames(POWER_strata_bio) <- seq_len(n_strata)
  POWER_strata_ct <- lapply(POWER_strata, function(x) x[[1]])
  POWER_strata_ct <- do.call(cbind, POWER_strata_ct)
  colnames(POWER_strata_ct) <- seq_len(n_strata)
  rownames(POWER_strata_ct) <- NULL
  # Target FDC
  FDC_bio <- FDC_pull(gene_CT_FDR,fdr_thred,ncell_type)
  FDC_bio_ct <- FDC_bio[[1]];  FDC_bio_m <- FDC_bio[[2]]
  simulation_results <- list(ct_TDR_bio = ct_tdr_bio, TDR_bio = tdr_bio,
                             ct_PWR_bio = ct_pwr, PWR_bio = ct_pwr_m,
                             PWR_strata_bio = POWER_strata_bio,
                             PWR_strata_ct_bio = POWER_strata_ct,
                             ct_FDC_bio = FDC_bio_ct, FDC_bio = FDC_bio_m)
  return(simulation_results)
}
