# Author: Andrew Chen, andrewac@pennmedicine.upenn.edu
# Adapted from code by Jean-Philippe Fortin, fortin946@gmail.com, available at
# https://github.com/Jfortin1/ComBatHarmonization
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license.

# Correcting Covariance Batch Effects: CovBat
# Performs original ComBat to residualize and harmonize observations followed by
# PCA step to harmonize covariance across sites
covbat <- function(x, # input data
                   bat, # vector of batch numbers
                   mod = NULL, # model matrix
                   percent.var = 80, # PC number selected to explain this % of variation
                   n.pc = NULL, # if specified, use this number of PCs instead
                   mean.only = FALSE, # scale parameter in initial ComBat, works better with scaling
                   score.parametric = TRUE, # parametric ComBat for scores
                   score.eb = FALSE, # empirical Bayes for scores
                   eb=TRUE, verbose=TRUE, parametric=TRUE,
                   resid=FALSE, # leave residualized
                   train = NULL # labels for train set used to determine betas
                   ) {
  # apply ComBat to remove mean/scale, keep intercept/linear pred on the side
  dat <- as.matrix(x)
  
  .checkConstantRows <- function(dat){
    sds <- rowSds(dat)
    ns <- sum(sds==0)
    if (ns>0){
      message <- paste0(ns, " rows (features) were found to be constant across samples. Please remove these rows before running ComBat.")
      stop(message)
    }
  }
  .checkConstantRows(dat)
  if (eb){
    if (verbose) cat("[combat] Performing ComBat with empirical Bayes\n")
  } else {
    if (verbose) cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
  }
  # make batch a factor and make a set of indicators for batch
  batch <- as.factor(bat)
  batchmod <- model.matrix(~-1+batch)  
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')
  
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- lapply(levels(batch), function(x)which(batch==x))
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  # check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])
  
  # Number of covariates or covariate levels
  if (verbose) cat("[combat] Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  # Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n.batch+1)){
      stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
  }
  
  ## Standardize Data across features
  if (verbose) cat('[combat] Standardizing Data across features\n')
  
  # Estimate coefficients using training set if specified, otherwise use full data
  if (!is.null(train)) {
    design_tr <- design[train,]
    
    B.hat1 <- solve(crossprod(design_tr))
    B.hat1 <- tcrossprod(B.hat1, design_tr)
    B.hat <- tcrossprod(B.hat1, dat[,train])
  } else {
    B.hat1 <- solve(crossprod(design))
    B.hat1 <- tcrossprod(B.hat1, design)
    B.hat <- tcrossprod(B.hat1, dat)
  }
  
  # Standardization Model
  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
  stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  
  if(!is.null(design)){
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp%*%B.hat)
  }	
  s.data <- (dat-stand.mean)/(tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  
  ##Get regression batch effect parameters
  if (eb){
    if (verbose) cat("[combat] Fitting L/S model and finding priors\n")
  } else {
    if (verbose) cat("[combat] Fitting L/S model\n")
  }
  batch.design <- design[,1:n.batch]
  gamma.hat <- tcrossprod(solve(crossprod(batch.design, batch.design)), batch.design)
  gamma.hat <- tcrossprod(gamma.hat, s.data)
  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,rowVars(s.data[,i], na.rm=TRUE)) # fixed error
  }
  
  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb){
    ##Find Priors
    #gamma.bar <- apply(gamma.hat, 1, mean)
    #t2 <- apply(gamma.hat, 1, var)
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apriorMat(delta.hat)
    b.prior <- bpriorMat(delta.hat)
    
    ##Find EB batch adjustments
    if (parametric){
      if (verbose) cat("[combat] Finding parametric adjustments\n")
      for (i in 1:n.batch){
        temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    } else {
      if (verbose) cat("[combat] Finding non-parametric adjustments\n")
      for (i in 1:n.batch){
        temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),gamma.hat[i,], delta.hat[i,])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    }
    
  } 
  
  if (mean.only) {
    delta.star <- array(1, dim = dim(delta.star))
  }
  
  ### Normalize the Data ###
  if (verbose) cat("[combat] Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    if (eb){
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/tcrossprod(sqrt(delta.star[j,]), rep(1,n.batches[j]))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/tcrossprod(sqrt(delta.hat[j,]), rep(1,n.batches[j]))
    }
    j <- j+1
  }
  
  comdata <- bayesdata
  x_pc <- prcomp(t(bayesdata)) # PC on ComBat-adjusted data
  
  # Subset scores based on percent of variance explained
  npc <- which(cumsum(x_pc$sdev/sum(x_pc$sdev)) > percent.var/100)[1]
  if (!is.null(n.pc)) {npc <- n.pc}
  # print(npc)
  scores <- x_pc$x[,1:npc]
  
  # ComBat without covariates to remove site effect in score mean/variance
  scores_com <- combat(t(scores), bat, eb = score.eb, parametric = score.parametric)
  full_scores <- x_pc$x
  full_scores[,1:npc] <- t(scores_com$dat.combat)
  
  # Project scores back into observation space
  x.covbat <- t(full_scores %*% t(x_pc$rotation))
  
  if (resid == FALSE) {
    x.covbat <- x.covbat * (tcrossprod(sqrt(var.pooled), rep(1,n.array)))+stand.mean
  } else {
    x.covbat <- x.covbat * (tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  }
  
  return(list(dat.covbat=x.covbat, 
              s.data=s.data,
              com.data=comdata,
              gamma.hat=gamma.hat, delta.hat=delta.hat, 
              gamma.star=gamma.star, delta.star=delta.star, 
              gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, b.prior=b.prior, batch=batch, mod=mod, 
              stand.mean=stand.mean, stand.sd=sqrt(var.pooled)[,1],
              scores.gamma=scores_com$gamma.hat, scores.delta=scores_com$delta.hat,
              npc=npc, x.pc = x_pc, com.scores = scores_com))
}