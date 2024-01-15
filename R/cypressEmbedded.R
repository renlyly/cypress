
# This function wraps up cypress_design and cypress_wrapper functions to
# facilitate the user-friendness of our power evaluation pipelines.

cypressEmbedded <- function(n_sim = 30, n_gene = 30000, DE_pct = 0.05,
                            ncell_type = 6, ss_group_set = c(10,20,50),
                            lfc_set = c(0, 0.5, 1, 1.5),
                            sim_param ,
                            lfc_target = 0.5, fdr_thred = 0.1){

  # design experiment
  exp_design <- cypress_design(n_sim, n_gene, DE_pct,
                               ncell_type, ss_group_set,
                               lfc_set, sim_param)
  # conduct simulation and metric evaluation.
  simulation_results <- cypress_wrapper(exp_design, lfc_target, fdr_thred)

  return(simulation_results)

}
