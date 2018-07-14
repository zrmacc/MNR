# Purpose: Fitting procedure for multivariate normal regression
# Updated: 180713

#' Fit Multivariate Outcome Model
#'
#' Fit the multivariate outcome regression model.
#' 
#' The target and surrogate model matrices are expected in numeric format.
#' Expand factors and interactions in advance. If an intercept is required,
#' include a vector of constants.
#'
#' @param yt Target outcome vector.
#' @param Ys Surrogate outcome matrix. 
#' @param Zt Target model matrix.
#' @param Zs Surrogate model matrix.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report number of iterations?
#'
#' @importFrom methods new
#' @importFrom stats coef resid
#' @export

fit.mnr = function(yt,Ys,Zt,Zs,maxit=10,eps=1e-6,report=T){
  # Dimensions
  n = nrow(Zt);
  m = ncol(Ys);
  # Objective function
  Q = function(S){-(n/2)*log(fastDet(S))};
  
  ## Initialize
  # List to hold initial parameters
  theta0 = list();
  # Initialize beta
  Mb = fitNorm(y=yt,Z=Zt);
  theta0$b = Mb$Beta;
  B = Mb$Ibb*Mb$Tau;
  # Vecotrize surrogate
  ys = c(fastT(Ys));
  # Initialize alpha
  theta0$a = alpha0(Zs,ys);
  # Initialize sigma
  theta0$E = cbind(yt-fastMMp(Zt,theta0$b),matrix(ys-fastMMp(Zs,theta0$a),ncol=m,byrow=T));
  theta0$S = fastIP(theta0$E,theta0$E)/n;
  # Initial objective
  q0 = Q(theta0$S);
  # Update wrapper
  Update = function(a,S){
    # Precision
    L = fastInv(S);
    # Update Beta
    b1 = MNRbeta(s=ys,Zt=Zt,B=B,Zs=Zs,b0=theta0$b,a=a,L=L);
    # Update Alpha
    a1 = MNRalpha(t=yt,s=ys,Zt=Zt,Zs=Zs,b=b1,L=L);
    # Residual matrix
    E1 = cbind(yt-fastMMp(Zt,b1),matrix(ys-fastMMp(Zs,a1),ncol=m,byrow=T)); 
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
    delta = q1-q0;
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
  L = fastInv(theta1$S);
  # Partition
  Ltt = L[1,1];
  Lts = L[1,2:(m+1),drop=F];
  Lss = L[2:(m+1),2:(m+1)];
  # Information matrices
  J = list();
  J$Ibb = infoMNR(n,Zt,Zt,Ltt);
  J$Iaa = infoMNR(n,Zs,Zs,Lss);
  J$Iba = infoMNR(n,Zt,Zs,Lts);
  # Residuals
  eT = theta1$E[,1];
  ES = theta1$E[,2:(m+1)];
  colnames(ES) = colnames(Ys);
  # Regression coefficients
  Alpha = as.numeric(theta1$a);
  names(Alpha) = colnames(Zs);
  Beta = as.numeric(theta1$b);
  names(Beta) = colnames(Zt);
  # Covariance matrix
  Sigma = theta1$S;
  rownames(Sigma) = colnames(Sigma) = c("t",colnames(Ys));
  # Output
  Out = new(Class="mnr",Coefficients=list("Beta"=Beta,"Alpha"=Alpha),
            Covariance=Sigma, Information=J, Residuals=list("Target"=eT,"Surrogate"=ES));
  return(Out);
}
