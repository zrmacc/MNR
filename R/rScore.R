# Purpose: Repeated Score test for multivariate outcome regression
# Updated: 180910

########################
# Multivariate case
########################

#' Multivariate Regression Repeated Score Test
#' 
#' Individually tests each column of \code{G} for association with the target 
#' outcome. \code{G} may contain missing elements, although the remaining model
#' matrices may not. \code{rScore.mnr} accelerates association testing by
#' recycling the same null model for each hypothesis test.
#' 
#' @param Y Outcome matrix.
#' @param j Column number of the outcome of interest. By default, \code{j=1}.
#' @param G Numeric matrix of covariates for the target outcome whose regression
#'   coefficients are zero under the null.
#' @param X List of model matrices, one for each outcome.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE.
#' @param parallel Run association testing in parallel? Must register parallel
#'   backend first.
#' 
#' @importFrom foreach foreach '%do%'
#' @importFrom plyr aaply
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' @return A numeric vector of p-values, one for each column of \code{G}.
#' @examples 
#' \dontrun{
#' # See rMNR or vignette for data generation
#' # Genotype matrix
#' G = replicate(2000,rbinom(n=1000,size=2,prob=0.25));
#' storage.mode(G) = "numeric";
#' # Introduce missingness
#' G[sort(sample(length(G),size=0.01*length(G)))] = NA;
#' # Repeated Score test
#' R = rScore.mnr(Y=Y,G=G,X=X,report=T);
#' # Estimated size
#' mean(R<=0.05);
#' }

rScore.mnr = function(Y,j=1,G,X,maxit=100,eps=1e-8,report=F,parallel=F){
  # Check input type
  if(!is.matrix(Y)){stop("A numeric matrix is expected for Y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")}
  if(!is.list(X)||(!all(unlist(lapply(X,is.matrix))))){
    stop("A list of numeric model matrices is expected for X.")};
  # Restructure
  k = ncol(Y);
  # Place target outcome first
  perm = c(j,seq(from=1,to=k)[-j]);
  Y = Y[,perm];
  X = X[perm];
  rm(j);
  # Check for missingness
  Miss = sum(is.na(Y))+sum(unlist(lapply(X,function(x){sum(is.na(x))})));
  if(Miss>0){stop("Inputs other than G should contain no missing data.")};
  # Fit null model
  M0 = fit.mnr(Y=Y,X=X,maxit=maxit,eps=eps,report=report);
  # Extract covariance
  L = vcov(M0,type="Outcome",inv=T);
  # Extract residuals
  E = resid(M0);
  # Working vector
  u = fastMMp(E,L[1,]);
  # Function to calculate score statistics
  aux = function(g){
    # Check for missingness
    key = !is.na(g);
    # If missingness, exclude 
    if(sum(!key)>0){
      u0 = u[key];
      g0 = g[key];
      X0 = lapply(X,function(x){x[key,,drop=F]});
      X1 = lapply(X,function(x){x[!key,,drop=F]});
    } else {
      u0 = u;
      g0 = g;
      X0 = X;
    }
    ## Information matrices
    # Marginal information
    I11 = L[1,1]*fastIP(g0,g0);
    i = NULL;
    I12 = foreach(i=1:k,.combine=cbind) %do% {
      return(fastIP(g0,X0[[i]])*L[1,i])
    };
    # Nusiance information
    I22 = vcov(M0,type="Information");
    if(sum(!key)>0){
      l = NULL;
      # Loss of information
      J22 = foreach(i=1:k,.combine=rbind) %do% {
        Slice = foreach(l=1:k,.combine=cbind) %do% {
          Out = fastIP(X1[[i]],X1[[l]])*L[i,l];
          return(Out);
        };
      };
      I22 = I22-J22;
    };
    # Efficient information
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    # Score vector
    U = as.numeric(fastIP(g0,u0));
    # Score statistic
    Ts = (U^2)/V;
    return(Ts);
  }
  # Score statistics
  S = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # P values
  p = pchisq(q=S,df=1,lower.tail=F);
  return(p);
}