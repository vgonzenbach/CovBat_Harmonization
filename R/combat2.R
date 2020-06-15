# Author: Jean-Philippe Fortin, fortin946@gmail.com
# This is a modification of the ComBat function code from the sva package that 
# can be found at https://bioconductor.org/packages/release/bioc/html/sva.html 
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license. 
# Code optimization improved by Richard Beare 

# Modified by Andrew Chen
# Added functionality to use only training data as input and to have 
# residualized observations as output

#' Combatting batch effects when combining batches of gene expression microarray 
#' data
#' 
#' \code{combat} is a modified version of the ComBat code written by 
#' Jean-Philippe Fortin available at 
#' \url{https://github.com/Jfortin1/ComBatHarmonization/}. The function
#' harmonizes the mean and variance of observations across sites under an
#' empirical Bayes framework. \code{combat} additionally includes options
#' to output residualized observations, estimate coefficients using a training
#' subset, and regress out unwanted confounders. \code{combat2} further includes
#' options to control for covariate effects in variance using 
#' \link[lmvar]{lmvar}
#'
#' @param dat 
#' @param batch 
#' @param mod An optional design matrix to preserve, usually output of 
#'    \link[stats]{model.matrix}.
#' @param var.mod An optional design matrix to preserve in the variance, usually
#'    output of \link[stats]{model.matrix}.
#' @param nuisance.mod An optional design matrix to regress out, usually output 
#'    of \link[stats]{model.matrix} without intercept.
#' @param train 
#' @param resid 
#' @param eb 
#' @param parametric 
#' @param mean.only 
#' @param verbose 
#'
#' @return
#' 
#' @import lmvar
#' @export
#'
#' @examples
#' 
#' @seealso Modification to ComBat to regress out unwanted confounders 
#' proposed by Wachinger et al. (2020), \url{https://arxiv.org/abs/2002.05049}
combat2 <- function(dat, batch, mod = NULL, var.mod = NULL, nuisance.mod = NULL, 
                   train = NULL, resid = FALSE, eb = TRUE, 
                   parametric = TRUE, mean.only = FALSE, verbose = FALSE)
{
  dat <- as.matrix(dat)
  
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
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')
  
  # add nuisance mod to mod, if specified
  if (!is.null(nuisance.mod)) {
    mod <- cbind(mod, nuisance.mod)
  }
  
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
  
  if (!is.null(var.mod)) {
    design.sigma <- cbind(batchmod,var.mod)
    # check for intercept in covariates, and drop if present
    check <- apply(design.sigma, 2, function(x) all(x == 1))
    design.sigma <- as.matrix(design.sigma[,!check])
  }
  
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
  
  varB.hat <- NULL
  
  # Standardization Model
  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
  sd.pooled <- (tcrossprod(sqrt(var.pooled), rep(1,n.array))) # rearrange as p x n
  stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  
  if(!is.null(design)){
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp%*%B.hat)
  }
  
  # Estimate variance effects if controlling for covariate effect in variance
  # then remove predicted values (including intercept) before harmonization
  p <- dim(dat)[1]
  if (!is.null(var.mod)) {
    # used to residualize out all effects except for site
    # tailored to fit output of lmvar
    design_bat <- design[,-n.batch] # drop last batch since is reference
    design_bat[,1:(n.batch-1)] <- 0 # remove batches
    design_bat <- cbind(`intercept` = rep(1, n.array), design_bat)
    
    design_sig_bat <- design.sigma[,-n.batch] # drop last batch since is reference
    design_sig_bat[,1:(n.batch-1)] <- 0 # remove batches
    design_sig_bat <- cbind(`intercept` = rep(1, n.array), design_sig_bat)
    
    var_fit_all <- list() # for debugging
    stand.mean <- matrix(0, p, n.array)
    sd.pooled <- matrix(0, p, n.array)
    
    # store mean/variance effects
    B.hat <- matrix(0, ncol(design_bat), p, 
                    dimnames = list(colnames(design_bat), rownames(dat)))
    varB.hat <- matrix(0, ncol(design_sig_bat), p, 
                       dimnames = list(colnames(design_sig_bat), rownames(dat)))
    
    for (i in 1:p) {
      var_fit <- lmvar(dat[i,], X_mu = design, X_sigma = design.sigma)
      B.hat[,i] <- var_fit$coefficients_mu
      varB.hat[,i] <- var_fit$coefficients_sigma

      stand.mean[i,] <- design_bat %*% var_fit$coefficients_mu
      sd.pooled[i,] <- exp(design_sig_bat %*% var_fit$coefficients_sigma)
      var_fit_all[[i]] <- var_fit # store for debugging
      
      if (!is.null(train)) {
        var_fit <- lmvar(dat[i, train], X_mu = design[train,],
                         X_sigma = design.sigma[train,])
        B.hat[,i] <- var_fit$coefficients_mu
        varB.hat[,i] <- var_fit$coefficients_sigma
        
        stand.mean[i,] <- design_bat %*% var_fit$coefficients_mu
        sd.pooled[i,] <- exp(design_sig_bat %*% var_fit$coefficients_sigma)
        
        var_fit_all[[i]] <- var_fit # store for debugging
      }
    }
  }
  
  s.data <- (dat-stand.mean)/sd.pooled
  
  ## Get regression batch effect parameters
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
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/
        tcrossprod(sqrt(delta.star[j,]), rep(1,n.batches[j]))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/
        tcrossprod(sqrt(delta.hat[j,]), rep(1,n.batches[j]))
    }
    j <- j+1
  }
  
  # If using nuisance.mod, reintroduce wanted covariates
  all.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  if (!is.null(nuisance.mod)) {
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    tmp[,(n.batch+dim(mod)[2]-dim(nuisance.mod)[2]):(dim(design)[2])] <- 0
    wanted.mean <- all.mean+t(tmp%*%B.hat)
  }
  
  # Reintroduce covariate effects and intercept
  if (resid) {
    bayesdata <- bayesdata*sd.pooled
  } else if (!is.null(nuisance.mod)) {
    bayesdata <- bayesdata*sd.pooled + wanted.mean
  } else {
    bayesdata <- bayesdata*sd.pooled + stand.mean
  }
  
  return(list(dat.combat=bayesdata,
              s.data=s.data,
              gamma.hat=gamma.hat, delta.hat=delta.hat,
              gamma.star=gamma.star, delta.star=delta.star,
              gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, b.prior=b.prior, 
              batch=batch, mod=mod,
              stand.mean=stand.mean, stand.sd=sd.pooled,
              B.hat = B.hat,
              varB.hat = varB.hat))
}
