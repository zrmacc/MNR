# Purpose: Generate data for multivarate normal regression package
# Updated: 180910

#' Simulate Multivariate Normal Data
#' 
#' Function to simulate outcomes from a multivariate normal regression model.
#' 
#' @param X List of design matrices, one for each outcome.
#' @param b List of regression coefficients, one for each outcome. 
#' @param S Outcome covariance structure. 
#'   
#' @importFrom mvnfast rmvn
#' @export
#' @return Numeric \eqn{n\times k} matrix, where \eqn{n} is the number of rows in each
#' design matrix, and \eqn{k} is the number of rows in the covariance structure. 
#' @examples 
#' \dontrun{
#' set.seed(100);
#' # Observations
#' n = 1e3;
#' ## Design matrices
#' X1 = cbind(1,matrix(rnorm(2*n),nrow=n));
#' colnames(X1) = c("int",paste0("x0",seq(1:2)));
#' X2 = cbind(1,matrix(rnorm(3*n),nrow=n));
#' colnames(X2) = c("int",paste0("x1",seq(1:3)));
#' X3 = cbind(1,matrix(rnorm(4*n),nrow=n));
#' colnames(X3) = c("int",paste0("x2",seq(1:4)));
#' X = list(X1,X2,X3);
#' # Target Parameter
#' b1 = c(-1,0.1,-0.1);
#' b2 = c(1,-0.1,0.1,0);
#' b3 = c(0,0.1,-0.1,0.1,-0.1);
#' b = list(b1,b2,b3);
#' # Exchangeable covariance structure
#' S = array(0.5,dim=c(3,3)) + 0.5*diag(3);
#' # Generate data
#' Y = rMNR(X=X,b=b,S=S);
#' }

rMNR = function(X,b,S){
  ## Input check
  if(!is.list(X)||(!all(unlist(lapply(X,is.matrix))))){
    stop("A list of numeric model matrices is expected for X.")};
  if(!is.list(b)||(!all(unlist(lapply(b,is.vector))))){
    stop("A list of numeric vectors is expected for b.")};
  if(!is.matrix(S)||(!isSymmetric(S))){
    stop("A numeric covariance matrix is expected for S.")};
  # Observations
  n = nrow(X[[1]]);
  # Outcome dimension
  k = ncol(S);
  # Linear predictors
  H = array(0,dim=c(n,0));
  # Surrogate predictors
  for(j in 1:k){
    H = cbind(H,MMP(X[[j]],b[[j]]));
  }
  # Residuals 
  E = rmvn(n=n,mu=rep(0,k),sigma=S);
  # Outcomes
  Y = H+E;
  ## Output
  colnames(Y) = paste0("y",seq(1:k));
  rownames(Y) = seq(1:n);
  return(Y);
}
