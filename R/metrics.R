### TDR CT specific
### The first function is to conduct the TDR calculation.
### The second function is to organizing the TDR calculation and pull out the
### FDR results.
TDR_CT <- function(i,top_gene,gene_CT_FDR){
  iter_fdr <- gene_CT_FDR %>% arrange(gene_CT_FDR[,1+i])
  ct_tdr <- apply(top_gene,1,function(x) sum(iter_fdr$DE_CT_bio[seq_len(x)] == i)/x)
  return(ct_tdr)
}
TDR_pull <- function(gene_CT_FDR,top_g,ncell_type){
  top_gene <- data.frame(top_g)
  ct_tdr_bio <- lapply(seq_len(ncell_type),TDR_CT,top_gene, gene_CT_FDR)
  ct_tdr_bio <- do.call(rbind,ct_tdr_bio)
  colnames(ct_tdr_bio) <- top_g
  return(ct_tdr_bio)
}

### POWER CT specific
POWER_CT <- function(i,gene_CT_FDR,true_DEG_per_ct,fdr_thred){
  ct_pwr_count <- sum((gene_CT_FDR[,i+1] < fdr_thred) & (gene_CT_FDR$DE_CT_bio == i))
  ct_pwr_count <- ifelse(length(ct_pwr_count)==0, NA, ct_pwr_count)
  ct_pwr_count_denominator <- true_DEG_per_ct[which(true_DEG_per_ct$Celltype == i),'count']
  ct_pwr_count_denominator <- ifelse(length(ct_pwr_count_denominator)==0, NA, ct_pwr_count_denominator)
  ct_pwr <- ct_pwr_count / ct_pwr_count_denominator
  ct_all <- c(ct_pwr_count,ct_pwr)
  return(ct_all)
}
POWER_pull <- function(gene_CT_FDR,fdr_thred,ncell_type){
  true_DEG_per_ct <- as.data.frame(table(gene_CT_FDR$DE_CT_bio))
  colnames(true_DEG_per_ct) <- c("Celltype","count")
  pwr_bio <- lapply(seq_len(ncell_type),POWER_CT,gene_CT_FDR,true_DEG_per_ct,fdr_thred)
  pwr_bio <- do.call(rbind,pwr_bio)
  rownames(pwr_bio) <- c(paste0("Celltype",seq_len(ncell_type)))
  ct_pwr_bio <- pwr_bio[,2]
  pwr_bio_m <- sum(pwr_bio[,1])/sum(gene_CT_FDR$DE_CT_bio != 0)
  power_out <- list(ct_pwr_bio,pwr_bio_m)
  return(power_out)
}

POWER_strata_pull <- function(i,gene_CT_FDR,fdr_thred,ncell_type){
  strata <- gene_CT_FDR$strata
  gene_CT_FDR_strata <- gene_CT_FDR %>% filter(strata == i)
  if (nrow(gene_CT_FDR_strata) == 0) return(list(rep(NA, ncell_type), NA))
  POWER_strata <- POWER_pull(gene_CT_FDR_strata,fdr_thred,ncell_type)
  return(POWER_strata)
}


### FDC CT specific
FDC_CT <- function(i,gene_CT_FDR,fdr_thred){
false_dis <- sum((gene_CT_FDR[,i+1] < fdr_thred) & (gene_CT_FDR$DE_CT_bio != i))
true_dis <- sum((gene_CT_FDR[,i+1] < fdr_thred) & (gene_CT_FDR$DE_CT_bio == i))
  fdc_ct <- c(false_dis,true_dis)
  return(fdc_ct)
}
FDC_pull <- function(gene_CT_FDR,fdr_thred, ncell_type){
  FDC_bio <- lapply(seq_len(ncell_type),FDC_CT,gene_CT_FDR,fdr_thred)
  FDC_bio <- do.call(rbind,FDC_bio)
  colnames(FDC_bio) <- c("FD","TD")
  rownames(FDC_bio) <- c(paste0("Celltype",seq_len(ncell_type)))
  ct_FDC_bio <- ifelse(FDC_bio[,2] == 0,NA,FDC_bio[,1]/FDC_bio[,2])
  FDC_bio_m <- sum(FDC_bio[,1])/sum(FDC_bio[,2])
  FDC_pull_out <- list(ct_FDC_bio,FDC_bio_m)
  return(FDC_pull_out)
}
