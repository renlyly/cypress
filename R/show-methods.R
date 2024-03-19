
setMethod("show", "cypress_out",
          function(object){
            cat("Cell type specific target TDR:", "\n")
            print(head(slot(object, "ct_TDR_bio_smry") ,4) )
            cat("\n")
            cat("Cell type specific target power:", "\n")
            print(head(slot(object, "ct_PWR_bio_smry") ,4) )
            cat("\n")
            cat("Cell type specific target power by gene expression stratification:", "\n")
            print(head( slot(object, "PWR_strata_ct_bio_smry") ,4) )
            cat("\n")
            cat("Cell type specific target FDC:", "\n")
            print(head(slot(object, "ct_FDC_bio_smry") ,4) )
          }
)

setMethod(f = "show", signature = "est_out",
          definition = function(object) {
            cat("CYPRESS Estimated Object\n")
            cat("-------------------------------\n")
            cat("Health Alpha: ", "\n")
            print(slot(object, "health_alpha"))
            cat("\n")
            cat("Case Alpha: ", "\n")
            print(slot(object, "case_alpha"))
            cat("\n")
            cat("Health LMean M: ",  "\n")
            print(slot(object, "health_lmean_m"))
            cat("\n")
            cat("LOD M: ", "\n")
            print(slot(object, "lod_m"))
            cat("\n")
          }
)
