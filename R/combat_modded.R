# Author: Jean-Philippe Fortin, fortin946@gmail.com
# This is a modification of the ComBat function code from the sva package that 
# can be found at https://bioconductor.org/packages/release/bioc/html/sva.html 
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license. 
# Code optimization improved by Richard Beare 

# Modified by Andrew Chen for covbat.R
# Added functionality to use only training data as input and to have 
# residualized observations as output

#' Combatting batch effects when combining batches of gene expression microarray 
#' data
#' 
#' \code{combat_modded} is a modified version of the ComBat code written by 
#' Jean-Philippe Fortin available at 
#' \url{https://github.com/Jfortin1/ComBatHarmonization/}. The function
#' harmonizes the mean and variance of observations across sites under an
#' empirical Bayes framework. \code{combat_modded} additionally includes options
#' to output residualized observations, estimate coefficients using a training
#' subset, and regress out unwanted confounders.
#'
#' @param dat 
#' @param batch 
#' @param mod An optional design matrix to preserve, usually output of 
#'    \link[stats]{model.matrix}.
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
#' @export
#'
#' @examples
#' 
#' @seealso Modification to ComBat to regress out unwanted confounders 
#' proposed by Wachinger et al. (2020), \url{https://arxiv.org/abs/2002.05049}
combat_modded <- function(dat, batch, mod=NULL, nuisance.mod = NULL, 
                          train = NULL, resid = FALSE, eb = TRUE, 
                          parametric = TRUE, mean.only = FALSE, verbose = TRUE)
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
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/tcrossprod(sqrt(delta.star[j,]), rep(1,n.batches[j]))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/tcrossprod(sqrt(delta.hat[j,]), rep(1,n.batches[j]))
    }
    j <- j+1
  }
  
  # Reintroduce wanted covariates
  all.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  if (!is.null(nuisance.mod)) {
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    tmp[,(n.batch+dim(mod)[2]-dim(nuisance.mod)[2]):(dim(design)[2])] <- 0
    wanted.mean <- all.mean+t(tmp%*%B.hat)
  }
  
  if (!is.null(nuisance.mod)) {
    bayesdata <- (bayesdata*(tcrossprod(sqrt(var.pooled), rep(1,n.array)))) +
      wanted.mean
  } else if (resid == FALSE) {
    bayesdata <- (bayesdata*(tcrossprod(sqrt(var.pooled), rep(1,n.array)))) +
      stand.mean
  } else {
    bayesdata <- bayesdata*(tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  }
  
  return(list(dat.combat=bayesdata,
              s.data=s.data,
              gamma.hat=gamma.hat, delta.hat=delta.hat,
              gamma.star=gamma.star, delta.star=delta.star,
              gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, b.prior=b.prior, 
              batch=batch, mod=mod,
              stand.mean=stand.mean, stand.sd=sqrt(var.pooled)[,1]),
  )
}