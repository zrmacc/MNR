# Purpose: Fitting procedure for multivariate outcome regression
# Updated: 180314

#' @useDynLib MNR
#' @importFrom Rcpp sourceCpp
NULL

########################
# Bivariate case
########################

#' Fit Bivariate Outcome Model
#'
#' Fit the bivariate outcome regression model.
#'
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome. 
#' @param Z.t Target model matrix.
#' @param Z.s Surrogate model matrix.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report number of iterations?
#'
#' @importFrom methods new
#' @importFrom stats coef resid
#' @importFrom RcppEigen fastLmPure
#' @export

fit.bnr = function(y.t,y.s,Z.t,Z.s,maxit=10,eps=1e-4,report=T){
  # Observations
  n = nrow(Z.s);
  p = ncol(Z.t);
  q = ncol(Z.s);
  # Initialize estimates
  M.t.0 = RcppEigen::fastLmPure(X=Z.t,y=y.t);
  b0 = coef(M.t.0);
  M.s.0 = RcppEigen::fastLmPure(X=Z.s,y=y.s)
  a0 = coef(M.s.0);
  E0 = cbind(resid(M.t.0),resid(M.s.0));
  s0 = fastIP(E0,E0)/n;
  # Objective function
  Q = function(s){-(n/2)*log(fastDet(s))};
  # Initial objective
  q0 = Q(s0);
  # Incomplete projections
  A.t = incP(Z.t);
  A.s = incP(Z.s);
  # Update wrapper
  Update = function(b0,a0,s0){
    Out = updateBVR(yt=y.t,ys=y.s,Zt=Z.t,Zs=Z.s,At=A.t,As=A.s,b0=b0,a0=a0,s0=s0);
  }
  # Update interations
  for(i in 1:maxit){
    # Initial objective
    q0 = Q(s0);
    # Proposal
    U = Update(b0,a0,s0);
    # Final objective
    q1 = Q(U$s1);
    # Update parameters if objective increases sufficiently
    if((q1-q0)>eps){
      s0 = U$s1;
      b0 = U$b1;
      a0 = U$a1;
    } else {break};
  }
  if(report){cat(paste0(i," updates performed before tolerance limit."),"\n")};
  # Information matrices
  Ibb = s0[1,1]*fastIP(Z.t,Z.t);
  Iaa = s0[2,2]*fastIP(Z.s,Z.s);
  Iba = s0[1,2]*fastIP(Z.t,Z.s);
  J = list("Ibb"=Ibb,"Iaa"=Iaa,"Iba"=Iba);
  # Residuals
  eT = as.numeric(y.t-Z.t%*%b0);
  eS = as.numeric(y.s-Z.s%*%a0);
  # Name regression coefficients
  Alpha = as.numeric(a0);
  names(Alpha) = colnames(Z.s);
  Beta = as.numeric(b0);
  names(Beta) = colnames(Z.t);
  # Name covariance matrix
  rownames(s0) = colnames(s0) = c("y.t","y.s");
  # Output
  Out = new(Class="mvr",Coefficients=list("Beta"=Beta,"Alpha"=Alpha),
            Covariance=s0, Information=J, Residuals=list("Target"=eT,"Surrogate"=eS));
  return(Out);
}