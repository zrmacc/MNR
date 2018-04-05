# Purpose: Descriptions of example data
# Updated: 180315

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
#'    \item{Z.t}{Model matrix for the target outcome.}
#'    \item{Z.s}{Model matrix for the surrogate outcome.}
#'    \item{Beta}{Regression coefficients used to generate the target outcome.}
#'    \item{Alpha}{Regression coefficients used to generate the surrogate outcome.}
#' }
"D.bnr"

#' Simulated Data with Triivariate Normal Outcomes
#' 
#' Example data for 1000 subjects. The outcomes are trivariate normal with unit 
#' variances and an exchangeable correlation structure, \eqn{\rho=0.5}. The 
#' target and first surrogate outcomes each depend on three independent N(0,1)
#' covariates, while the second surrogate outcome depends on four such
#' covariates. \code{Beta} contains regression coefficients for the target
#' outcome, and \code{Alpha} contains regression coefficients for the secondary
#' ("surrogate") outcomes.
#' 
#' @format A list of vectors and matrices. 
#' \describe{ 
#'   \item{y.t}{Numeric vector containing the target outcome.} 
#'   \item{y.s}{Numeric matrix containing the surrogate outcomes.} 
#'   \item{Z.t}{Model matrix for the target outcome.} 
#'   \item{L.s}{List of model matrices for the surrogate outcomes.} 
#'   \item{Beta}{Regression coefficients used to generate the target outcome.} 
#'   \item{Alpha}{Regression coefficients used to generate the surrogate 
#'   outcome.} }
"D.mnr"