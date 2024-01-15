##Unit Test using RUnit



#test quickPower and plot functions
test_quickPower <- function() {
  set.seed(123)
  Quick_power <- quickPower(data = "IAD")

  checkEquals(ncol(Quick_power@ct_TDR_bio_smry),10)
  checkEquals(ncol(Quick_power@TDR_bio_smry), 9)
  checkEquals(ncol(Quick_power@ct_PWR_bio_smry), 9)
  checkEquals(ncol(Quick_power@PWR_bio_smry), 4)
  checkEquals(ncol(Quick_power@PWR_strata_bio_smry), 12)
  checkEquals(ncol(Quick_power@PWR_strata_ct_bio_smry), 13)
  checkEquals(ncol(Quick_power@ct_FDC_bio_smry), 9)
  checkEquals(ncol(Quick_power@FDC_bio_smry), 4)


  p1 <- plotPower(Quick_power,sample_size=10)

  p2 <- plotTDR(Quick_power,sample_size=10)

  p3 <- plotFDC(Quick_power,sample_size=10)


  checkTrue(all(c(length(p1),
                  length(p2),
                  length(p3)) > 0))

}
