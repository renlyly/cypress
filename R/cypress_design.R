#########################################################
### exp_design: this function is to design the experiment.
# 'n_sim': total number of iterations users wish to conduct. Default to 30.
#
# 'n_gene': total number of genetic features users with to conduct.
#  Default to 30000.
#
# 'DE_pct': percentage of DEG on each cell type. Default to 0.05.
#
# 'ncell_type': total number of cell types. Should match with simulation
#  parameters provided in 'sim_param'.
#
# 'lfc_set': effect sizes users wish to simulate. The length should be less than
#  or equal to 5. Default to 0,0.5,1,1.5.
#
# 'ss_group_set': sample sizes per group users wish to simulate. The length
#  should be less than or equal to 5. Default to 10,20,50.
#
# ‘sim_param’: A S4 object includes all the necessary simulation parameters.
#  'health_lmean_m', 'health_lmean_d', 'lod_m', 'lod_d',
#   'health_alpha' and  'case_alpha' should all be specified.
# exp_design  = cypress_design(n_sim = 3, n_gene = 30000, DE_pct = 0.05,
#                              ncell_type = 6, ss_group_set = c(10,20,50),
#                              lfc_set = c(0,0.5, 1, 1.5),
#                              sim_param = GSE60424_param)
# exp_design
# 1 healthy 2 case
cypress_design <- function(n_sim, n_gene, DE_pct,
                           ncell_type, ss_group_set,
                           lfc_set,
                           sim_param){
  if(!(n_sim > 0)) stop("Total number of iterations should be positive.")
  if(n_sim > 100) stop("Total number of iterations needs be below 100 to save computation time.")
  if(n_gene < 1000) stop("Total number of simulated genes be 1000 or above.")
  if(length(ss_group_set)>5 || length(lfc_set)>5)
    stop("The length of both ss_group_set and lfc_set input should not be greater than 5")
  if(max(ss_group_set)>100) stop("The maximum sample size should not exceed 100.")
  if(is.null(sim_param)) stop("Users need to specify simulation parameters.")
  if(!all(lfc_set >= 0)) stop("LFC should be non-negative")

  # lfc_sd is fixed
  lfc_sd <- 0.5
  # strata is fixed
  exp_strata <- c(-Inf, 10, 20, 40, 80,
                 160, 320, 640, 1280, Inf)

  scenarios <- expand.grid(ss_group = ss_group_set,
                           lfc = lfc_set)
  n_scenarios <- dim(scenarios)[1]
  n_strata <- length(exp_strata)-1
  # mean and var/cov of distribution mean
  health_lmean_m <- sim_param@health_lmean_m
  health_lmean_d <- sim_param@health_lmean_d
  # mean and var/cov of distribution dispersion
  lod_m <- sim_param@lod_m
  lod_d <- sim_param@lod_d
  # proportion alpha parameter.
  health_alpha <- sim_param@health_alpha
  case_alpha <- sim_param@case_alpha
  if(!(ncell_type == length(health_alpha) &
     ncell_type == length(case_alpha) &
     ncell_type == length(health_lmean_m) &
     ncell_type == nrow(health_lmean_d) &
     ncell_type == length(lod_m) &
     ncell_type == nrow(lod_d))){
 stop("The number of cell types is not match with dimensions of simulation parameters")
  }
  design_use <- design_in(n_sim = n_sim, n_gene = n_gene, DE_pct = DE_pct,
                          ncell_type = ncell_type,
                          scenarios = scenarios, n_scenarios =  n_scenarios,
                          n_strata = n_strata, health_lmean_m = health_lmean_m,
                          health_lmean_d = health_lmean_d,
                          lod_m = lod_m, lod_d = lod_d,
                          health_alpha = health_alpha, case_alpha = case_alpha)
  message("Complete experiment design for CYPRESS.")
  return(design_use)
}
