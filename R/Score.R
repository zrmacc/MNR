# Purpose: Score test for multivariate outcome regression
# Updated: 180910

########################
# Multivariate case
########################

#' Multivariate Regression Score Test
#' 
#' For an outcome of interest among \code{Y}, tests the hypothesis that a subset
#' of the regression coefficients for that outcome are fixed at the reference 
#' value \code{b10}. In particular, suppose \eqn{\beta} denotes the regression 
#' coefficient for the target outcome. Partition 
#' \eqn{\beta=(\beta_{1},\beta_{2})}. \code{Score.bnem} performs a score test of
#' \eqn{H_{0}:\beta_{1}=\beta_{10}}.
#' 
#' @param Y Outcome matrix.
#' @param j Column number of the outcome of interest. By default, \code{j=1}. 
#' @param X List of model matrices, one for each outcome.
#' @param L Logical vector, with as many entires as columns in the target
#'   design matrix, indicating which columns design are fixed under the null.
#' @param b10 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE. 
#' 
#' @importFrom foreach foreach '%do%'
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' @return A numeric vector containing the score statistic, degrees of freedom,
#'   and p-value.
#' @examples 
#' \dontrun{
#' # See rMNR or vignette for data generation
#' # Test b13 = 0, which is FALSE
#' Score.mnr(Y=Y,j=1,X=X,L=c(F,F,T));
#' # Test b24 = 0, which is TRUE
#' Score.mnr(Y=Y,j=2,X=X,L=c(F,F,F,T));
#' # Test b32 = ... = b35 = 0, which is FALSE
#' Score.mnr(Y=Y,j=3,X=X,L=c(F,T,T,T,T));
#' # Test b32 = b34 = 0.1, which is TRUE
#' Score.mnr(Y=Y,j=3,X=X,b10=c(0.1,0.1),L=c(F,T,F,T,F));
#' }

Score.mnr = function(Y,j=1,X,L,b10=NULL,maxit=100,eps=1e-8,report=F){
  # Check input type
  if(!is.matrix(Y)){stop("A numeric matrix is expected for Y.")};
  if(!is.list(X)||(!all(unlist(lapply(X,is.matrix))))){
    stop("A list of numeric model matrices is expected for X.")};
  # Restructure
  k = ncol(Y);
  # Place target outcome first
  perm = c(j,seq(from=1,to=k)[-j]);
  Y = Y[,perm];
  X = X[perm];
  
  # Check test specification
  Z = X[[1]];
  p = ncol(Z);
  if(length(L)!=p){stop("L should have as many entries as columns the design matrix for the outcome of interest.")};
  df = sum(L);
  if(df==0){stop("At least 1 entry of L should be TRUE.")};
  if(df==p){stop("At least 1 entry of L should be FALSE.")};
  
  # Check for missingness
  Miss = sum(is.na(Y))+sum(unlist(lapply(X,function(x){sum(is.na(x))})));
  if(Miss>0){stop("Inputs should contain no missing data.")};
  # Null coefficient
  if(is.null(b10)){b10=rep(0,times=df)};
  
  # Partition target design
  Za = Z[,L,drop=F];
  Zb = Z[,!L,drop=F];
  X[[1]] = Zb;
  # Adjust response for fixed component
  Y[,1] = Y[,1]-fastMMp(Za,b10);
  # Fit null model
  M0 = fit.mnr(Y=Y,X=X,maxit=maxit,eps=eps,report=report);
  # Extract covariance
  L = vcov(M0,type="Outcome",inv=T);
  # Extract residuals
  E = resid(M0);
  # Score vector
  U = fastIP(Za,fastMMp(E,L[1,]));
  # Covariance matrix
  I11 = fastIP(Za,Za)*L[1,1];
  I22 = vcov(M0,type="Information",inv=F);
  # Cross Information
  i = NULL;
  I12 = foreach(i=1:k,.combine=cbind) %do% {
    return(fastIP(Za,X[[i]])*L[1,i])
  }
  # Efficient information
  V = SchurC(I11=I11,I22=I22,I12=I12);
  # Score statistic
  Ts = fastQF(X=U,A=fastInv(V));
  # P value
  p = pchisq(q=Ts,df=df,lower.tail=F);
  # Output
  Out = c("Score"=Ts,"df"=df,"p"=p);
  return(Out);
}