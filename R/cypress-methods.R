## show function
setGeneric("getcypress", function(object, name) {standardGeneric("getcypress")})
setGeneric("setcypress", function(object, name, value) {standardGeneric("setcypress")})



setMethod("show", "cypress_out",
          function(object){
            cat("Cell type specific target TDR:", "\n")
            print(head(getcypress(object, "ct_TDR_bio_smry") ,4) )
            cat("\n")
            cat("Cell type specific target power:", "\n")
            print(head(getcypress(object, "ct_PWR_bio_smry") ,4) )
            cat("\n")
            cat("Cell type specific target power by gene expression stratification:", "\n")
            print(head( getcypress(object, "PWR_strata_ct_bio_smry") ,4) )
            cat("\n")
            cat("Cell type specific target FDC:", "\n")
            print(head(getcypress(object, "ct_FDC_bio_smry") ,4) )
          }
)

setMethod(f = "show", signature = "est_out",
          definition = function(object) {
            cat("CYPRESS Estimated Object\n")
            cat("-------------------------------\n")
            cat("Health Alpha: ", "\n")
            print(getcypress(object, "health_alpha"))
            cat("\n")
            cat("Case Alpha: ", "\n")
            print(getcypress(object, "case_alpha"))
            cat("\n")
            cat("Health LMean M: ",  "\n")
            print(getcypress(object, "health_lmean_m"))
            cat("\n")
            cat("LOD M: ", "\n")
            print(getcypress(object, "lod_m"))
            cat("\n")
          }
)

## accessor


setMethod("getcypress", "cypress_out", function(object, name) {
  slot(object, name)
})

setMethod("getcypress", "est_out", function(object, name) {
  slot(object, name)
})

setMethod("getcypress", "design_in", function(object, name) {
  slot(object, name)
})

## replace methods

setMethod("setcypress", "cypress_out", function(object, name, value) {
  checkmate::assertString(name)
  slot(object, name) <- value
  validObject(object)
  return(object)
})

setMethod("setcypress", "est_out", function(object, name, value) {
  checkmate::assertString(name)
  slot(object, name) <- value
  validObject(object)
  return(object)
})

setMethod("setcypress", "design_in", function(object, name, value) {
  checkmate::assertString(name)
  slot(object, name) <- value
  validObject(object)
  return(object)
})

## Validity

setValidity("cypress_out", function(object) {
  msg <- NULL
  slotin <- c("ct_TDR_bio_smry", "TDR_bio_smry", "ct_PWR_bio_smry", "PWR_bio_smry",
                 "PWR_strata_bio_smry", "PWR_strata_ct_bio_smry", "ct_FDC_bio_smry", "FDC_bio_smry")


  actualSlotNames <- slotNames(object)
  if (!all(slotin %in% actualSlotNames)) {
    missingSlots <- slotin[!slotin %in% actualSlotNames]
    msg <- c(msg, paste("Missing slots:", paste(missingSlots, collapse=", ")))
  }

  # Check for dimensions
  for (slotName in slotin) {
    if (is.null(getcypress(object, slotName))) next  # Skip if slot is NULL
    if (nrow(getcypress(object, slotName)) < 4 || ncol(getcypress(object, slotName)) < 4) {
      msg <- c(msg, paste("'", slotName, "' must have more than 3 rows and 3 columns.", sep=""))
    }
  }

  # Check for required column names
  requiredColumns <- c("ss", "lfc", "lfc_sd")
  additionalColumns <- list(
    ct_TDR_bio_smry = "ct",
    PWR_strata_ct_bio_smry = "ct"
  )

  for (slotName in slotin) {
    currentCols <- colnames(getcypress(object, slotName))
    # Check for general required columns
    if (!all(requiredColumns %in% currentCols)) {
      missingCols <- requiredColumns[!requiredColumns %in% currentCols]
      msg <- c(msg, paste("'", slotName, "' is missing columns:", paste(missingCols, collapse=", "), sep=""))
    }
    # Check for additional required columns in specific slots
    if (slotName %in% names(additionalColumns)) {
      additionalCol <- additionalColumns[[slotName]]
      if (!additionalCol %in% currentCols) {
        msg <- c(msg, paste("'", slotName, "' is missing the '", additionalCol, "' column.", sep=""))
      }
    }
  }

  if (length(msg) == 0) TRUE else msg
})


setValidity("est_out", function(object) {
  msg <- NULL

  requiredSlots <- c("health_alpha", "case_alpha", "health_lmean_m", "health_lmean_d",
                     "lod_m", "lod_d")
  actualSlotNames <- slotNames(object)

  if (!all(requiredSlots %in% actualSlotNames)) {
    missingSlots <- requiredSlots[!requiredSlots %in% actualSlotNames]
    msg <- c(msg, paste("Missing required slots:", paste(missingSlots, collapse=", ")))
  }

  # Check lengths for numeric variables
  numericSlots <- c("health_alpha", "case_alpha", "health_lmean_m", "lod_m")
  for (slotName in numericSlots) {
    if (length(getcypress(object, slotName)) < 3) {
      msg <- c(msg, paste("'", slotName, "' length must be at least 3.", sep=""))
    }
  }

  # Check matrix dimensions for equality
  matrixSlotsEqualDims <- c("lod_d", "health_lmean_d")
  for (slotName in matrixSlotsEqualDims) {
    matrix <- getcypress(object, slotName)
    if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
      msg <- c(msg, paste("'", slotName, "' must be a square matrix.", sep=""))
    }
  }
  if (length(msg) == 0) TRUE else msg
})


setValidity("design_in", function(object) {
  msg <- NULL

  requiredSlots <- c("health_alpha", "case_alpha", "health_lmean_m", "health_lmean_d",
                     "lod_m", "lod_d")
  actualSlotNames <- slotNames(object)

  if (!all(requiredSlots %in% actualSlotNames)) {
    missingSlots <- requiredSlots[!requiredSlots %in% actualSlotNames]
    msg <- c(msg, paste("Missing required slots:", paste(missingSlots, collapse=", ")))
  }

  # Check lengths for numeric variables
  numericSlots <- c("health_alpha", "case_alpha", "health_lmean_m", "lod_m")
  for (slotName in numericSlots) {
    if (length(getcypress(object, slotName)) < 3) {
      msg <- c(msg, paste("'", slotName, "' length must be at least 3.", sep=""))
    }
  }

  # Check matrix dimensions for equality
  matrixSlotsEqualDims <- c("lod_d", "health_lmean_d")
  for (slotName in matrixSlotsEqualDims) {
    matrix <- getcypress(object, slotName)
    if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
      msg <- c(msg, paste("'", slotName, "' must be a square matrix.", sep=""))
    }
  }


  if (length(msg) == 0) TRUE else msg
})






