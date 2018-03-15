# Purpose: Descriptions of example data
# Updated: 180314

#' Simulated Data with Bivariate Normal Outcomes
#' 
#' Example data for 1000 subjects. The outcomes are bivariate normal with unit
#' variances and correlation \eqn{\rho=0.5}. Each outcome depends on a design matrix containing three independent 
#' N(0,1) covariates. \code{Beta} contains regression coefficients for the target outcome, and \code{Alpha} contains
#' regression coefficients for the secondary ("surrogate") outcome.
#' 
#' @format A list of vectors and matrices. 
#' \describe{ 
#'    \item{y.t}{Numeric vector containing the target outcome.}
#'    \item{y.s}{Numeric vector containing the surrogate outcome.}
#'    \item{D.t}{Design matrix for the target outcome.}
#'    \item{D.s}{Design matrix for the surrogate outcome.}
#'    \item{Beta}{Regression coefficients used to generate the target outcome.}
#'    \item{Alpha}{Regression coefficients used to generate the surrogate outcome.}
#' }
"D.bvr"