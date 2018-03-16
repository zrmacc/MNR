# Purpose: Fitting procedure for multivariate outcome regression
# Updated: 180315

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
#' @importFrom RcppEigen fastLmPure
#' @export

fit.bnr = function(y.t,y.s,Z.t,Z.s,maxit=10,eps=1e-4,report=T){
  # Ensure matrix
  y.t = as.matrix(y.t,nrow=nrow(Z.t));
  y.s = as.matrix(y.s,nrow=nrow(Z.s));
  # Observations
  n = nrow(Z.s);
  p = ncol(Z.t);
  q = ncol(Z.s);
  # Initialize beta
  M.t.0 = RcppEigen::fastLmPure(X=Z.t,y=y.t);
  b0 = coef(M.t.0);
  # Initialize alpha
  M.s.0 = RcppEigen::fastLmPure(X=Z.s,y=y.s)
  a0 = coef(M.s.0);
  # Initialize sigma
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
      E0 = U$E1;
    } else {break};
  }
  if(report){cat(paste0(i," updates performed before tolerance limit."),"\n")};
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
  rownames(s0) = colnames(s0) = c("y.t","y.s");
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
#' @importFrom RcppEigen fastLmPure
#' @export

fit.mnr = function(y.t,y.s,Z.t,Z.s,maxit=10,eps=1e-4,report=T){
  # Ensure matrix
  y.t = as.matrix(y.t,nrow=nrow(Z.t));
  # Observations
  n = nrow(Z.t);
  # Surrogates
  m = ncol(y.s);
  # Covariates
  p = ncol(Z.t);
  q = ncol(Z.s);
  # Initialize beta
  M.t.0 = RcppEigen::fastLmPure(X=Z.t,y=y.t);
  b0 = coef(M.t.0);
  # Form surrogate vector
  S = c(fastT(y.s));
  # Initialize alpha
  M.s.0 = RcppEigen::fastLmPure(X=Z.s,y=S);
  a0 = coef(M.s.0);
  # Initialize sigma
  E0 = cbind(resid(M.t.0),matrix(resid(M.s.0),ncol=2,byrow=T));
  s0 = fastIP(E0,E0)/n;
  # Objective function
  Q = function(s){-(n/2)*log(fastDet(s))};
  # Initial objective
  q0 = Q(s0);
  # Incomplete projections
  A.t = incP(Z.t);
  # Update wrapper
  Update = function(b0,a0,s0){
    Out = updateMNR(yt=y.t,S=S,Zt=Z.t,Zs=Z.s,At=A.t,b0=b0,a0=a0,s0=s0);
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
      E0 = U$E1;
    } else {break};
  }
  if(report){cat(paste0(i," updates performed before tolerance limit."),"\n")};
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
  rownames(s0) = colnames(s0) = c(colnames(y.t),colnames(y.s));
  # Output
  Out = new(Class="mnr",Coefficients=list("Beta"=Beta,"Alpha"=Alpha),
            Covariance=s0, Information=J, Residuals=list("Target"=eT,"Surrogate"=eS));
  return(Out);
}