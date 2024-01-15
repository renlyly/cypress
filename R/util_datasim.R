remove_top_n_rows <- function(data, n){
  celltypes <- colnames(data)
  # Use lapply to iterate through the columns (cell types)
  top_indices_list <- lapply(celltypes, function(ct) {
    top_indices <- head(order(-data[, ct]), n)
    return(top_indices)
  })
  # Combine the top indices from all cell types
  removed_gene <- unlist(top_indices_list)
  return(removed_gene)
}

# This function is to simulate effect sizes (log-fold-changes) matrix.
# LFC is drawing from a normal distribution with specified mean and sd.
# The dimension is G by K. If a feature on a specific cell types is
# up or down regulated, then it has a simu
lfc_simulator <- function(n_gene, DE_pct, ncell_type,
                          lfc_mean, lfc_sd){
  DE_n <- n_gene*DE_pct
  DE_total <- DE_n*ncell_type
  CT_code <- rep(seq_len(ncell_type), each = DE_n)
  lfc <- rnorm(n = DE_total, mean = lfc_mean, sd = lfc_sd)
  lfc_celltype <- cbind(CT_code,lfc)
  # Create DE change amount matrix for the first 1-1500 genes (first half up, remain down)
  lfc_ct_mtx <- apply(lfc_celltype,1, function(x){
    each_DE <- c(rep(0,ncell_type))
    each_DE[x[1]] <- x[2]
    return(each_DE)})
  lfc_ct_mtx <- rbind(t(lfc_ct_mtx),
                      matrix(0,nrow = (n_gene-DE_total),
                             ncol = ncell_type))
  return(lfc_ct_mtx)
}

# This function is to simulate cell type specific reference panel. Users should
# provide mean vector, var-cov matrix, dispesion vector, var-cov matrix, and
# the matrix of effect size. This function returns reference panels respectively
# for the case and control groups.
cs_ref_panel_simulator <- function(ncell_type, n_gene,
                                   health_lmean_m, health_lmean_d,
                                   lod_m, lod_d,lfc_ct_mtx){
sim_health_lmean_trim <-MASS::mvrnorm(n_gene*1.5,health_lmean_m,health_lmean_d)
  sim_lod_trim <- -abs(MASS::mvrnorm(n_gene*1.5, lod_m, lod_d))
  # select the extreme parameters
  gene_out1 <- unique(remove_top_n_rows(sim_health_lmean_trim, 50))
  gene_out2 <- unique(remove_top_n_rows(-sim_lod_trim, 50))
  gene_out <- unique(c(gene_out1,gene_out2))
  # filter out extreme parameters
  sim_health_lmean_trim <- sim_health_lmean_trim[-gene_out,]
  sim_lod_trim <- sim_lod_trim[-gene_out,]
  # dimensionality check
  sim_health_lmean_trim <- sim_health_lmean_trim[seq_len(n_gene),]
  sim_lod_trim <- sim_lod_trim[seq_len(n_gene),]
  # creat case group
  sim_case_lmean_trim <- sim_health_lmean_trim + lfc_ct_mtx
  # exponentiate back to original scale
  sim_case_mean_trim <- exp(sim_case_lmean_trim)
  sim_health_mean_trim <- exp(sim_health_lmean_trim)
  sim_d_trim <- exp(sim_lod_trim)
  # gamma process: create arrays to expedite simulations
  sim_case_mean_trim_array <- array(c(sim_case_mean_trim, sim_d_trim),
                                    dim = c(n_gene, ncell_type, 2))
  sim_health_mean_trim_array <- array(c(sim_health_mean_trim, sim_d_trim),
                                      dim = c(n_gene, ncell_type, 2))

  cellT_expr_healthy  <- apply(sim_health_mean_trim_array, seq_len(2), function(x)
    rgamma(1, shape = 1/x[2],scale = x[1]*x[2]))
  rownames(cellT_expr_healthy) <- seq_len(n_gene)
  cellT_expr_case <- apply(sim_case_mean_trim_array, seq_len(2), function(x)
    rgamma(1, shape = 1/x[2],scale = x[1]*x[2]))
  rownames(cellT_expr_case) <- seq_len(n_gene)
  cs_ref_panel_all <- list(cellT_expr_healthy, cellT_expr_case)
  return(cs_ref_panel_all)
}

# This function is to simulate cell type proportion, sum up cell type specific
# underlying expression panel, and simulate RNAseq count data.
simulate_sample <- function(i,cellT_expr_healthy,cellT_expr_case,
                            health_alpha,case_alpha,
                            ncell_type,n_gene) {
  cell_prop_health <- sirt::dirichlet.simul(alpha = t(matrix(health_alpha)))
  cell_prop_case <- sirt::dirichlet.simul(alpha = t(matrix(case_alpha)))
  cellT_expr_healthy_mix <- cellT_expr_healthy %*% t(cell_prop_health)
  cellT_expr_case_mix <- cellT_expr_case %*% t(cell_prop_case)
  cellT_expr_healthy_count <- vapply(cellT_expr_healthy_mix,
                                     function(x) rpois(1, x), integer(1))
  cellT_expr_case_count <- vapply(cellT_expr_case_mix,
                                  function(x) rpois(1, x), integer(1))

  cell_prop <- rbind(cell_prop_health, cell_prop_case)
  rownames(cell_prop) <- c(paste0("control", i), paste0("case", i))
  colnames(cell_prop) <- paste0("Celltype",seq( seq_len(ncell_type)))

  cellT_expr_count <- cbind(cellT_expr_healthy_count, cellT_expr_case_count)
  colnames(cellT_expr_count) <- c(paste0("control", i), paste0("case", i))
  rownames(cellT_expr_count) <- seq(seq_len(n_gene))

  sample_i <- list(cell_prop = cell_prop, cellT_expr_count = cellT_expr_count)
  return(sample_i)
}


perform_deconvolution <- function(cellT_expr, group, RNAseq_final_count,
                                  refinx,ncell_type) {
  selected_columns <- grep(paste0("^",group), colnames(RNAseq_final_count),
                           value = TRUE)
  sig_matrix_raw <- as.data.frame(cellT_expr[refinx,])
mixture_file_raw <- as.data.frame(RNAseq_final_count[refinx, selected_columns])
  result2 <- CIBERSORT(sig_matrix = sig_matrix_raw,
                       mixture_file = mixture_file_raw,
                       perm = 0, QN = TRUE, absolute = FALSE,
                       abs_method = 'sig.score')
  est_CT_prop <- result2[, seq_len(ncell_type)]
  return(est_CT_prop)
}



# This function summarize the ground true information used for later metric
# calculation.
gt_organizer <- function(RNAseq_final_count, lfc_ct_mtx, ncell_type){
  gene_CT_DE_connect <- data.frame(
    gene_id = as.numeric(rownames(RNAseq_final_count)),
    gene_mean =  rowMeans(RNAseq_final_count),
    DE_CT = ifelse(lfc_ct_mtx==0,0,1) %*% seq_len(ncell_type),
    LFC = lfc_ct_mtx %*% rep(1,ncell_type),
    strata = cut(rowMeans(RNAseq_final_count),
                 c(-Inf, 10, 20, 40, 80, 160, 320, 640, 1280, Inf),
                 right = TRUE, labels = c(seq_len(9))))
  return(gene_CT_DE_connect)
}

# cut(rnorm(1000,0,100),)
