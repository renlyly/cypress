

cypress_prop_unk <- function(Y, cell_type = NULL) {
  if (is.null(Y) | is.null(cell_type)) {
    stop("Y or cell_type is not provided!")
  }

  outT <- TOAST::myRefFreeCellMix(log(Y+1)/max(log(Y+1)),
        mu0=myRefFreeCellMixInitialize(log(Y+1)/max(log(Y+1)), K = cell_type),
                           verbose = FALSE)# in toast RefFree
  est_CT_prop_T <- abs(round(outT$Omega, 4))
  est_CT_prop_T <- sweep(est_CT_prop_T,1,1/apply(est_CT_prop_T,1,sum),FUN = "*")
  return(est_CT_prop_T)
}


cypress_prop_decon <- function(final_count, design1, sample_CT_prop) {
  if (missing(final_count) | missing(design1)) {  # Check input
    stop("final_count or design1 is not provided!")
  }
  if (is.null(sample_CT_prop)) {# Check sample_CT_prop
    stop("sample_CT_prop must be provided and cannot be  NULL!")
  } else if (is.numeric(sample_CT_prop) & length(sample_CT_prop) == 1) {
 #If sample_CT_prop is numeric, call cypress_prop_unk function
 sample_CT_prop <- cypress_prop_unk(Y = final_count, cell_type = sample_CT_prop)
  } else if (!is.matrix(sample_CT_prop)) {
    stop("sample_CT_prop must be a matrix or a single number!")
  }

  if (is.numeric(design1)) {  # check design
    design1 <- data.frame(disease = design1)
    if (missing(sample_CT_prop) || is.null(sample_CT_prop)) {
      rownames(design1) <- c(paste("C", seq_len(nrow(design1))  ))
    } else {
      rownames(design1) <- rownames(sample_CT_prop)
    }
  }
  if (is.data.frame(design1)) {
    if (any(is.na(match(c("disease"), colnames(design1)))) ||
        any(grepl("^[0-9]", colnames(design1)))) {
      colnames(design1) <- "disease"
    }
    if (any(is.na(match(rownames(design1), seq_len(nrow(design1)) ))) ||
        any(grepl("^[0-9]", rownames(design1)))) {
      if (missing(sample_CT_prop) || is.null(sample_CT_prop) ) {
        rownames(design1) <- c(paste("C", seq_len(nrow(design1))  ))
      } else {
        rownames(design1) <- rownames(sample_CT_prop)
      }
    }
  }
  if (any(rowSums(sample_CT_prop) < 0.995 | rowSums(sample_CT_prop) > 1.005)) {
    stop("The row sum for Cell type proportion must be between 0.99 and 1, please check the matrix")
  }
 tca.seq <- TCA::tca(X = final_count, W = sample_CT_prop,
                 refit_W = FALSE,   C1 = design1,   verbose = FALSE)
  Z_hat <- TCA::tensor(as.matrix(final_count), tca.seq, verbose = FALSE)
  nsample_each_group <- ncol(final_count)
  ncell_type <- ncol(sample_CT_prop)
  Z_hat_ary <- array(data = do.call(c, Z_hat),
                     dim = c(nrow(final_count), nsample_each_group, ncell_type))
  return(Z_hat_ary)
}

cypress_prop_est <- function(Z_hat_ary, genename = NULL, celltypename = NULL) {
  G <- dim(Z_hat_ary)[1]; N <- dim(Z_hat_ary)[2]
  K <- dim(Z_hat_ary)[3]
  # check the dim for names
  if (!is.null(genename) && length(genename) != G) {
    stop("Length of genename must match the first dimension of Z_hat_ary!")
  }
  if (!is.null(celltypename) && length(celltypename) != K) {
    stop("Length of celltypename must match the third dimension of Z_hat_ary!")
  }
  # initial output
  health_lmean <- matrix(NA, nrow = G, ncol = K)
  health_lod <- matrix(NA, nrow = G, ncol = K)
  rownames(health_lmean) <- genename
  colnames(health_lmean) <- celltypename
  rownames(health_lod) <- genename
  colnames(health_lod) <- celltypename
  vec <- seq_len(N)
  even_indices <- seq(2, length(vec), by = 2)

  results <- vapply(seq_len(K), function(i) {
    Z_hat_ary_mod <- pmax(Z_hat_ary[, even_indices - 1, i], 0) # negative to 0
    round_Z_hat <- round(Z_hat_ary_mod)
    Proper_param <- estParam1(as.matrix(round_Z_hat), type = 1)
    list(lmeans = Proper_param$lmeans, lOD = Proper_param$lOD)
  }, FUN.VALUE = list(lmeans = numeric(N), lOD = numeric(N)))
  health_lmean <- do.call(cbind, results[1,])
  health_lod <- do.call(cbind, results[2,])
  health_lmean[is.na(health_lmean) | health_lmean == -Inf] <- NA # deal with NA
  col_lmeans <- colMeans(health_lmean, na.rm = TRUE)
  health_lmean <- ifelse(is.na(health_lmean),
                         col_lmeans[col(health_lmean)], health_lmean)
  col_lod <- colMeans(health_lod, na.rm = TRUE)
    return(list(health_lmean = health_lmean, health_lod = health_lod))
}


cypress_prop_alpha <-function(final_count, sample_CT_prop = NA, cell_type = NA){
  if ((missing(final_count) & missing(sample_CT_prop)) |
      (!missing(final_count) & !missing(sample_CT_prop))) {  # input check
    stop("Only one of final_count or sample_CT_prop should be provided!")
  }
  if (!missing(final_count)) {
    if (missing(cell_type)) {
      stop("cell_type must be provided if final_count is used!")
    }
    sample_CT_prop <- cypress_prop_unk(Y = final_count, cell_type = cell_type)
  }
  health_alpha <- sirt::dirichlet.mle(x = sample_CT_prop)
  return(list(alpha = health_alpha$alpha, alpha0 = health_alpha$alpha0))
}

estParam1 <- function (X, type = 1)
{
  if (! ( (c( "matrix") %in% class(X)) & (c( "matrix") %in% class(X)) )   )
    stop("input X should be either an ExpressionSet or a matrix")
  if (!is.matrix(X))
    X <- Biobase::exprs(X)
  if (!(type == 1 | type == 2))
    stop("type must be either 1 or 2")
  if (type == 1) {
    res <- getDisp1(X)
  }
  if (type == 2) {
    res <- getDisp2(X)
  }
  res
}

getDisp1<-function(X){
  seqDepth <-colSums(X)
  k<-seqDepth/median(seqDepth)
  X2<-sweep(X, 2, k, FUN="/")##pseudo Counts

  m<-rowMeans(X2)
  v<-MatrixGenerics::rowVars(X2)
  phi.g <- (v-m)/m^2
  phi.g0  <- (v-m)/m^2
  ## only keep those with good coverage
  Good<-m>30 & rowMeans(X>0)>0.8
  phi.g0<-phi.g0[Good]
  phi.g0<-replace(phi.g0,is.na(phi.g0),.001)
  phi.g0<-replace(phi.g0,phi.g0<.001,0.001)
  ii<-(phi.g<=0.001) | is.na(phi.g)
  #set.seed(seed)
  phi.g[ii]<-sample(phi.g0, sum(ii), replace=TRUE)
  list(seqDepth=seqDepth,lmeans=log(m),lOD=log(phi.g))
}

getDisp2<-function(X){
  X<-edgeR::DGEList(counts=X,lib.size=colSums(X))
  y<-edgeR::estimateCommonDisp(X)
  y<- (estimateTrendedDisp(y))
  y<-edgeR::estimateTagwiseDisp(y)
  lmean<-y$AveLogCPM*log(2)+log(y$pseudo.lib.size/1e6)
  phi.g<-y$tagwise.dispersion
  list(seqDepth=y$sample$lib.size, lmean=lmean, lOD=log(phi.g))

}

cypress_prop_trim<- function(final_count,lower_prop=0.05,
                             upper_prop=0.95,lower_d=NULL){

quantile1 <- quantile(as.matrix(final_count), probs = c(lower_prop, upper_prop))

  Rowfilter <- which(apply(final_count, 1,
                           function(x) all(x < min(quantile1[[2]], 100000 )  ) ) )
  Rowfilter1<-findRefinx(as.matrix(final_count),
                         nmarker =  round(0.5*nrow(final_count) ))
  return(final_count[intersect(Rowfilter,Rowfilter1),])

}


cypress_prop <- function(Count_matrix, design , sample_CT_prop = NULL,
                         ncell_type = NULL) {
  Count_matrix<-cypress_prop_trim(Count_matrix)
  genename <- rownames(Count_matrix);  samplename <- colnames(Count_matrix)
  G <- nrow(Count_matrix);  N <- ncol(Count_matrix)
  method <- ifelse(is.null(ncell_type) & is.null(sample_CT_prop), 1,
                   ifelse(is.null(sample_CT_prop) & !is.null(ncell_type), 2, 3))
  if (is.null(design)) {
    stop("'design' must not be NULL")
  }
  if (is.vector(design)) {
    if (length(design) != N || !all(design %in% c(1, 2))) {
      stop("'design' must be a numeric vector of length N containing only 1 and 2")
    }
  } else if (is.data.frame(design)) {
    if (nrow(design) != N || !all(unlist(design) %in% c(1, 2))) {
      stop("'design' must be a data frame with N rows containing only 1 and 2")
    }
    design <- design[, 1]
  } else {
    stop("'design' must be a numeric vector or a data frame")
  }
  if (method == 1) {
    stop("Cell type information is missing")
  }
  if (method == 2) {  # method 2 needs to calculate sample_CT_prop
    K <- ncell_type
    sample_CT_prop <- cypress_prop_unk(Y = Count_matrix, cell_type = K)
    colnames(sample_CT_prop)<-paste0("Celltype",colnames(sample_CT_prop))
  }
  CTname <- colnames(sample_CT_prop);  K <- ncol(sample_CT_prop)
  Z_hat_ary <- cypress_prop_decon(final_count = Count_matrix,
             design1 = design, sample_CT_prop = sample_CT_prop)
  est_mean_lod <- cypress_prop_est(Z_hat_ary = Z_hat_ary,
                 genename = genename, celltypename = CTname)
  health_lmean_m <- colMeans(est_mean_lod$health_lmean)
  health_lmean_d <- cov(est_mean_lod$health_lmean)
lod_m <- colMeans(est_mean_lod$health_lod);lod_d <- cov(est_mean_lod$health_lod)
names(health_lmean_m)<-CTname; names(lod_m)<-CTname
rownames(lod_d)<-CTname;rownames(health_lmean_d)<-CTname
colnames(lod_d)<-CTname;colnames(health_lmean_d)<-CTname
  health_sample_CT_prop <- sample_CT_prop[which(design == 1),, drop = FALSE]
  case_sample_CT_prop <- sample_CT_prop[which(design == 2),, drop = FALSE]
  health_alpha <- cypress_prop_alpha(sample_CT_prop = health_sample_CT_prop)
  case_alpha <- cypress_prop_alpha(sample_CT_prop = case_sample_CT_prop)
result<-est_out(health_alpha = health_alpha$alpha,case_alpha = case_alpha$alpha,
  health_lmean_m=health_lmean_m, health_lmean_d=health_lmean_d,
  lod_m=lod_m,lod_d=lod_d,sample_CT_prop = sample_CT_prop,
  genename = genename,samplename = samplename, CTname = CTname,
  dimensions_Z_hat_ary = dim(Z_hat_ary) )
  return(result)}


# read S4 to list
S4tolist<-function(INPUT = NULL, CT_index = NULL,CT_unk= FALSE ){
  if (is.null(INPUT) || !inherits(INPUT, "SummarizedExperiment")) {
    stop("'INPUT' must be a SummarizedExperiment object")
  }
  Count_matrix <- assay(INPUT)
  design <- colData(INPUT)$disease
    if (CT_unk) {
      sample_CT_prop <- NULL
    } else {
      if (is.null(CT_index)) {
        sample_CT_prop <- as.data.frame(SummarizedExperiment::colData(INPUT))[, seq_len(ncol(SummarizedExperiment::colData(INPUT)) )]
        sample_CT_prop <- as.matrix(sample_CT_prop)
        } else {
        if (length(CT_index) < 3) {
          stop("Error: length of 'CT index' must be >= 3")
        }
        sample_CT_prop <- as.data.frame(SummarizedExperiment::colData(INPUT))[, CT_index]
        sample_CT_prop <- as.matrix(sample_CT_prop)
      }
    }
  ncell_type <- NULL
    if (CT_unk) {
      if (is.numeric(CT_index) && length(CT_index) == 1) {
        ncell_type <- CT_index
      } else {
        ncell_type <- length(CT_index)
      }
    } else {
      ncell_type <- ncol(sample_CT_prop)
    }
  if (ncell_type < 3) {
    stop("Cell type number must >= 3")
  }
    return(list(Count = Count_matrix, design = design,
                sample_CT_prop = sample_CT_prop, ncell_type = ncell_type))

 }


cypress_est <- function(INPUT = NULL, CT_index = NULL,CT_unk= FALSE ) {
  List0<-S4tolist(INPUT, CT_index,CT_unk)
  design1 <- as.data.frame(List0$design)
  colnames(design1) <- "disease"
  rownames(design1) <- colnames(List0$Count)
output<-cypress_prop(Count_matrix = as.matrix(List0$Count),
                     design = List0$design,
               sample_CT_prop=List0$sample_CT_prop,ncell_type=List0$ncell_type)

  return(output)
}


