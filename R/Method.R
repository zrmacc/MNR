# Purpose: Fitting procedure for multivariate outcome regression
# Updated: 180404

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
#' @param y.t Target outcome vector.
#' @param y.s Surrogate outcome vector. 
#' @param Z.t Target model matrix.
#' @param Z.s Surrogate model matrix.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report number of iterations?
#'
#' @importFrom methods new
#' @importFrom stats coef resid
#' @export

fit.bnr = function(y.t,y.s,Z.t,Z.s,maxit=10,eps=1e-6,report=T){
  # Dimensions
  n = nrow(Z.t);
  p = ncol(Z.t);
  q = ncol(Z.s);
  # Initialize beta
  B = incP(Z.t);
  b0 = fastMvp(B,y.t);
  # Initialize alpha
  A = incP(Z.s);
  a0 = fastMvp(A,y.s);
  # Initialize sigma
  E0 = cbind(y.t-fastMvp(Z.t,b0),y.s-fastMvp(Z.s,a0));
  s0 = fastIP(E0,E0)/n;
  # Objective function
  Q = function(s){-n*log(fastDet(s))};
  # Initial objective
  q00 = q0 = Q(s0);
  # Update wrapper
  Update = function(b,a,s){
    Out = updateBVR(t=y.t,s=y.s,X=Z.t,Xi=Z.s,B=B,A=A,b0=b,a0=a,s0=s);
  }
  # Update interations
  for(i in 1:maxit){
    # Proposal
    U = Update(b0,a0,s0);
    # Final objective
    q1 = Q(U$s1);
    # Update parameters if objective increases sufficiently
    d = q1-q0;
    if(d>eps){
      if(report){cat("Log likelihood increment:",signif(d,2),"\n")};
      q0 = q1;
      s0 = U$s1;
      b0 = U$b1;
      a0 = U$a1;
      E0 = U$E1;
    } else {break};
  }
  if(report){cat(paste0(i-1," updates performed before tolerance limit."),"\n")};
  # Precision matrix
  L = fastInv(s0);
  # Partition precision
  Ltt = L[1,1];
  Lts = L[1,2];
  Lss = L[2,2];
  # Information matrices
  Ibb = Ltt*fastIP(Z.t,Z.t);
  Iaa = Lss*fastIP(Z.s,Z.s);
  Iba = Lts*fastIP(Z.t,Z.s);
  J = list("Ibb"=Ibb,"Iaa"=Iaa,"Iba"=Iba);
  # Residuals
  eT = E0[,1,drop=F];
  colnames(eT) = colnames(y.t);
  eS = E0[,2,drop=F];
  colnames(eS) = colnames(y.s);
  # Name regression coefficients
  Alpha = as.numeric(a0);
  names(Alpha) = colnames(Z.s);
  Beta = as.numeric(b0);
  names(Beta) = colnames(Z.t);
  # Name covariance matrix
  rownames(s0) = colnames(s0) = c("t","s");
  # Output
  Out = new(Class="mnr",Coefficients=list("Beta"=Beta,"Alpha"=Alpha),
            Covariance=s0, Information=J, Residuals=list("Target"=eT,"Surrogate"=eS));
  return(Out);
}

########################
# Multivariate case
########################

#' Fit Multivariate Outcome Model
#'
#' Fit the multivariate outcome regression model.
#'
#' @param y.t Target outcome vector.
#' @param y.s Surrogate outcome matrix. 
#' @param Z.t Target model matrix.
#' @param Z.s Surrogate model matrix.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report number of iterations?
#'
#' @importFrom methods new
#' @importFrom stats coef resid
#' @export

fit.mnr = function(y.t,y.s,Z.t,Z.s,maxit=10,eps=1e-4,report=T){
  # Observations
  n = nrow(Z.t);
  # Surrogates
  m = ncol(y.s);
  # Covariates
  p = ncol(Z.t);
  q = ncol(Z.s);
  # Initialize beta
  B = incP(Z.t);
  b0 = fastMvp(B,y.t);
  # Form surrogate vector
  s = c(fastT(y.s));
  # Initialize alpha
  a0 = alpha0(Z.s,s);
  # Initialize sigma
  E0 = cbind(y.t-fastMvp(Z.t,b0),matrix(s-fastMvp(Z.s,a0),ncol=m,byrow=T));
  s0 = fastIP(E0,E0)/n;
  # Objective function
  Q = function(s){-(n/2)*log(fastDet(s))};
  # Initial objective
  q0 = Q(s0);
  # Update wrapper
  Update = function(b,a,s0){
    Out = updateMNR(t=y.t,s=s,X=Z.t,Xi=Z.s,B=B,b0=b,a0=a,s0=s0);
  }
  # Update interations
  for(i in 1:maxit){
    # Proposal
    U = Update(b0,a0,s0);
    # Final objective
    q1 = Q(U$s1);
    # Update parameters if objective increases sufficiently
    d = q1-q0;
    if(d>eps){
      if(report){cat("Log likelihood increment:",signif(d,2),"\n")};
      q0 = q1;
      s0 = U$s1;
      b0 = U$b1;
      a0 = U$a1;
      E0 = U$E1;
    } else {break};
  }
  if(report){cat(paste0(i-1," updates performed before tolerance limit."),"\n")};
  # Precision matrix
  L = fastInv(s0);
  # Partition
  Ltt = L[1,1];
  Lts = L[1,2:(m+1),drop=F];
  Lss = L[2:(m+1),2:(m+1)];
  # Information matrices
  Ibb = infoMNR(n,Z.t,Z.t,Ltt);
  Iaa = infoMNR(n,Z.s,Z.s,Lss);
  Iba = infoMNR(n,Z.t,Z.s,Lts);
  J = list("Ibb"=Ibb,"Iaa"=Iaa,"Iba"=Iba);
  # Residuals
  eT = E0[,1,drop=F];
  colnames(eT) = colnames(y.t);
  eS = E0[,2:(m+1)];
  colnames(eS) = colnames(y.s);
  # Name regression coefficients
  Alpha = as.numeric(a0);
  names(Alpha) = colnames(Z.s);
  Beta = as.numeric(b0);
  names(Beta) = colnames(Z.t);
  # Name covariance matrix
  rownames(s0) = colnames(s0) = c("t",colnames(y.s));
  # Output
  Out = new(Class="mnr",Coefficients=list("Beta"=Beta,"Alpha"=Alpha),
            Covariance=s0, Information=J, Residuals=list("Target"=eT,"Surrogate"=eS));
  return(Out);
}