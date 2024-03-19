
## License: http://cibersort.stanford.edu/CIBERSORT_License.txt





CoreAlg <- function(X, y, absolute, abs_method){


  svn_itor <- 3
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=FALSE)
    model
  }
    out <- parallel::mclapply(seq_len(svn_itor), res)#, mc.cores=svn_itor
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  t <- 1
  while(t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- stats::cor(k, y)
    t <- t + 1
  }

  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)] <-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q))
  if(absolute && abs_method == 'no.sumto1') w <- q
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}


doPerm <- function(perm, X, Y, absolute, abs_method){

  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    ## random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    ## standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    ## run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)

    mix_r <- result$mix_r

    ## store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}



CIBERSORT <- function(sig_matrix, mixture_file, perm, QN = TRUE, absolute,
                      abs_method='sig.score'){

  if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 0){
    stop("None identical gene between eset and reference had been found.
      Check your eset using: intersect(rownames(eset), rownames(reference))")
  }

  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score')
    stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")

  X<- sig_matrix

  Y <- tibble::rownames_to_column(mixture_file,var = "symbol")
  ## to prevent crashing on duplicated gene symbols, add unique numbers to name
  dups <- dim(Y)[1] - length(unique(Y[,1]))
  if(dups > 0) {
    warning(dups," duplicated gene symbol(s) found in mixture file!",sep="")
    rownames(Y) <- make.names(Y[,1], unique=TRUE)
  }else {rownames(Y) <- Y[,1]}
  Y <- Y[,-1]
  X <- data.matrix(X)
  Y <- data.matrix(Y)

  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm

  ## anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  ## quantile normalization of mixture file
    if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  ## store original mixtures
  Yorig <- Y
  Ymedian <- max(median(Yorig),1)

  ## intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  ## standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  ## empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',
                                         sep=""))

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  ## iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    ## standardize mixture
    y <- (y - mean(y)) / sd(y)

    ## run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)

    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    if(absolute && abs_method == 'sig.score') {
      w <- w * median(Y[,itor]) / Ymedian
    }

    ## calculate p-value
    if(P > 0) {pval <- 1 -(which.min(abs(nulldist - mix_r)) / length(nulldist))}

    ## print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) out <- c(out, sum(w))
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }


  ## return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
  else{
    colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",
                          paste('Absolute score (',abs_method,')',sep=""))}
  obj
}
