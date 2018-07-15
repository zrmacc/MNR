# Purpose: Fitting procedure for bivariate normal regression
# Updated: 180714

#' @useDynLib MNR
#' @importFrom Rcpp sourceCpp
NULL

#' Fit Bivariate Outcome Model
#'
#' Fit the bivariate outcome regression model.
#' 
#' The target and surrogate model matrices are expected in numeric format. 
#' Expand factors and interactions in advance. If an intercept is required, 
#' include a vector of constants. The inputs should contain no missing values.
#'
#' @param yt Target outcome vector.
#' @param ys Surrogate outcome vector. 
#' @param Zt Target model matrix.
#' @param Zs Surrogate model matrix.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report fitting progress? Default is FALSE. 
#'
#' @importFrom methods new
#' @importFrom stats coef resid
#' @export
#' @return An object of class \code{mnr} containing the fitted regression 
#'   parameters, the estimated target-surrogate covariance matrix, the 
#'   regression parameter information matrix, and the estimated residuals.

fit.bnr = function(yt,ys,Zt,Zs,maxit=100,eps=1e-8,report=T){
  # Dimensions
  n = nrow(Zt);
  # Objective function
  Q = function(S){-n*log(fastDet(S))};
  
  ## Initialize
  # List to hold initial parameters
  theta0 = list();
  # Initialize beta
  Mb = fitNorm(y=yt,Z=Zt);
  theta0$b = Mb$Beta;
  # Initialize alpha
  Ma = fitNorm(y=ys,Z=Zs);
  theta0$a = Ma$Beta;
  # Initialize sigma
  theta0$E = cbind(Mb$Resid,Ma$Resid);
  theta0$S = fastIP(theta0$E,theta0$E)/n;
  # Pre-computable structures
  B = Mb$Ibb*Mb$Tau;
  A = Ma$Ibb*Ma$Tau;
  # Update wrapper
  Update = function(a,S){
    # Precision
    L = fastInv(S);
    # Update Beta
    b1 = BNRbeta(s=ys,Zt=Zt,B=B,Zs=Zs,b0=theta0$b,a=a,L=L);
    # Update Alpha
    a1 = BNRalpha(t=yt,Zt=Zt,Zs=Zs,A=A,b=b1,a0=theta0$a,L=L);
    # Residual matrix
    E1 = cbind(yt-fastMMp(Zt,b1),ys-fastMMp(Zs,a1));
    # Update covariance
    S1 = fastIP(E1,E1)/n;
    # Output
    Out = list("b"=b1,"a"=a1,"E"=E1,"S"=S1);
    return(Out);
  }
  
  # Initialize objective
  q0 = q1 = Q(theta0$S);
  # Set current parameters to initial values
  theta1 = theta0;

  # Update interations
  for(i in 1:maxit){
    # Proposal
    U = Update(a=theta1$a,S=theta1$S);
    # Final objective
    q2 = Q(U$S);
    # Update parameters if objective increases sufficiently
    delta = q2-q1;
    if(delta>eps){
      if(report){cat("Log likelihood increment:",signif(delta,2),"\n")};
      # Update objective
      q1 = q2;
      # Update parameters
      theta1 = U;
    } else {break};
  }
  
  ## Fitting report
  if(report){
    if(i<maxit){
      cat(paste0(i-1," available data update(s) performed before tolerance limit."),"\n");
    } else {
      cat(paste0(i," update(s) performed without reaching tolerance limit."));
    }
  };
  
  # Precision matrix
  theta1$L = fastInv(theta1$S);
  # Information matrices
  J = list();
  J$Ibb = theta1$L[1,1]*B;
  J$Iba = theta1$L[1,2]*fastIP(Zt,Zs);
  J$Iaa = theta1$L[2,2]*A;
  # Residuals
  eT = theta1$E[,1];
  eS = theta1$E[,2];
  # Regression coefficients
  Alpha = as.numeric(theta1$a);
  names(Alpha) = colnames(Zs);
  Beta = as.numeric(theta1$b);
  names(Beta) = colnames(Zt);
  # Covariance matrix
  Sigma = theta1$S;
  rownames(Sigma) = colnames(Sigma) = c("t","s");
  # Output
  Out = new(Class="mnr",Coefficients=list("Beta"=Beta,"Alpha"=Alpha),
            Covariance=Sigma, Information=J, Residuals=list("Target"=eT,"Surrogate"=eS));
  return(Out);
}