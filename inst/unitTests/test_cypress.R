##Unit Test using RUnit

#test quickPower and plot functions
test_quickPower <- function() {
  set.seed(123)
  Quick_power <- quickPower(data = "IAD")

  checkEquals(ncol(getcypress(Quick_power, "ct_TDR_bio_smry")), 10)
  checkEquals(ncol(getcypress(Quick_power, "TDR_bio_smry")), 9)
  checkEquals(ncol(getcypress(Quick_power, "ct_PWR_bio_smry")), 9)
  checkEquals(ncol(getcypress(Quick_power, "PWR_bio_smry")), 4)
  checkEquals(ncol(getcypress(Quick_power, "PWR_strata_bio_smry")), 12)
  checkEquals(ncol(getcypress(Quick_power, "PWR_strata_ct_bio_smry")), 13)
  checkEquals(ncol(getcypress(Quick_power, "ct_FDC_bio_smry")), 9)
  checkEquals(ncol(getcypress(Quick_power, "FDC_bio_smry")), 4)



  p1 <- plotPower(Quick_power,sample_size=10)

  p2 <- plotTDR(Quick_power,sample_size=10)

  p3 <- plotFDC(Quick_power,sample_size=10)


  checkTrue(all(c(length(p1),
                  length(p2),
                  length(p3)) > 0))

}





#test simulation functions
test_sim <- function() {

  data(ASD_prop_se)
  length_ct <-(seq_len(6) + 2)
  nsim<-2
  ss_groupset<-c(8,10)

  test <- simFromData(INPUTdata = ASD_prop, CT_index = length_ct, CT_unk = FALSE,
                        n_sim = nsim,n_gene = 1000, DE_pct = 0.05,
                        ss_group_set = ss_groupset,
                        lfc_set = c(1, 1.5))

  test2 <- simFromParam(sim_param="IAD",n_sim = 2,DE_pct = 0.05,n_gene = 1000,
                         ss_group_set = c(8, 10),
                         lfc_set = c(1, 1.5),
                         lfc_target = 0.5, fdr_thred = 0.1)


  checkEquals(ncol(slot(test,"ct_FDC_bio_smry")),length(length_ct)+3)
  checkEquals(ncol(slot(test,"PWR_strata_ct_bio_smry")),13)
  checkEquals(nrow(slot(test2,"ct_TDR_bio_smry")),length(length_ct)*nsim*2)
  checkEquals(unique(slot(test2,"ct_PWR_bio_smry")$ss),ss_groupset )



}




