# Author: Andrew Chen, andrewac@pennmedicine.upenn.edu
# Adapted from code by Jean-Philippe Fortin, fortin946@gmail.com, available at
# https://github.com/Jfortin1/ComBatHarmonization
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license.

# Correcting Covariance Batch Effects: CovBat
# Performs original ComBat to residualize and harmonize observations followed by
# PCA step to harmonize covariance across sites

#' Correcting Covariance Batch Effects: CovBat
#' 
#' \code{covbat2} harmonizes covariance of observations across sites. It
#' first applies ComBat to harmonize mean and variance then adjusts variance of
#' PC scores to harmonize covariance. This version includes an experimental
#' option to control for covariates in the harmonization of PC scores. Also, 
#' this version of CovBat uses easier-to-follow code relying on 
#' \link[CovBat]{combat2}
#'
#' @param x A \emph{p x n} matrix (or object coercible by \link[base]{as.matrix}
#'   to a numeric matrix) of observations where \emph{p} is the number of
#'   features and \emph{n} is the number of subjects.
#' @param bat Factor (or object coercible by \link[base]{as.factor} to a 
#'    factor) designating batch IDs.
#' @param mod Optional design matrix of covariates to preserve, usually from 
#'    \link[stats]{model.matrix}.
#' @param score.mod Optional design matrix of covariates to preserve in the
#'    PC score variance, usually from \link[stats]{model.matrix}.
#' @param percent.var Numeric. The number of harmonized principal component 
#'    scores is selected to explain this proportion of the variance.
#' @param n.pc Optional numeric. If specified, this number of principal
#'    component scores is harmonized. Overrides \code{percent.var}.
#' @param train Optional logical vector specifying subset of observations used
#'    to estimate coefficients.
#' @param mean.only If \code{TRUE}, ComBat step does not harmonize variance.
#' @param std.var If \code{TRUE}, scales variances to be equal to 1 before PCA.
#' @param resid Whether to leave intercept and covariates regressed out.
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization.
#' @param parametric If \code{TRUE}, uses parametric adjustments in ComBat step
#' @param score.eb If \code{TRUE}, uses ComBat model with empirical Bayes for
#'   harmonization of scores.
#' @param score.parametric If \code{TRUE}, uses parametric ComBat is used
#'   for harmonization of scores.
#' @param verbose Whether additional details are printed to console.
#'
#' @return \code{covbat2} returns a list containing the following components:
#' \item{dat.covbat}{Harmonized data as a matrix with same dimensions as 
#'    \code{x}}
#' \item{combat.out}{List output of \link[CovBat]{combat} from the ComBat step}
#' \item{combat.scores}{List output of \link[CovBat]{combat} from the score
#'    harmonization step}
#' \item{n.pc}{Number of principal components harmonized}
#' \item{prcomp.out}{PCA results after ComBat, output of \link[stats]{prcomp}}
#' 
#' @examples
covbat2 <- function(x, bat, mod = NULL, score.mod = NULL, percent.var = 0.95, 
                    n.pc = NULL, train = NULL, mean.only = FALSE, 
                    std.var = TRUE, resid = FALSE, eb = TRUE, parametric = TRUE,
                    score.eb = FALSE, score.parametric = TRUE, verbose = FALSE)
{
  # Use ComBat to remove mean/variance effects
  com_out <- combat2(x, bat, mod = mod, mean.only = mean.only, eb = eb,
                      parametric = parametric, resid = TRUE, verbose = verbose)
  x_com <- com_out$dat.combat
  
  # PC on ComBat-adjusted data
  x_pc <- prcomp(t(x_com), center = TRUE, scale. = std.var)
  
  # Subset scores based on percent of variance explained
  npc <- which(cumsum(x_pc$sdev^2/sum(x_pc$sdev^2)) > percent.var)[1]
  if (!is.null(n.pc)) {npc <- n.pc}
  # print(npc)
  scores <- x_pc$x[,1:npc]
  
  # ComBat without covariates to remove site effect in score mean/variance
  scores_com <- combat2(t(scores), bat, eb = score.eb, var.mod = score.mod,
                        parametric = score.parametric, verbose = verbose)
  full_scores <- x_pc$x
  full_scores[,1:npc] <- t(scores_com$dat.combat)
  
  # Project scores back into observation space
  if (std.var) {
    x_covbat <- t(full_scores %*% t(x_pc$rotation)) * 
      matrix(x_pc$scale, dim(x_com)[1], dim(x_com)[2]) + 
      matrix(x_pc$center, dim(x_com)[1], dim(x_com)[2])
  } else {
    x_covbat <- t(full_scores %*% t(x_pc$rotation)) + 
      matrix(x_pc$center, dim(x_com)[1], dim(x_com)[2])
  }
  # Return to original centering
  if (resid == FALSE) {
    x_covbat <- x_covbat + com_out$stand.mean
  }
  
  return(list(dat.covbat = x_covbat, 
              combat.out = com_out,
              combat.scores = scores_com,
              npc=npc, x.pc = x_pc))
}
