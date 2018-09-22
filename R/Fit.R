# Purpose: Fitting procedure for multivariate normal regression
# Updated: 180910

#' @useDynLib MNR
#' @importFrom Rcpp sourceCpp
NULL

#' Fit Multivariate Outcome Model
#' 
#' Fits a regression model in which a multivariate normal random vector is 
#' observed for each subject. Regression models are specified using a list of
#' numeric matrices, one for each column of \code{Y}. 
#' 
#' @param Y Outcome matrix.
#' @param X List of model matrices, one for each outcome.
#' @param sig Significance level, for confidence intervals. 
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report fitting progress? Default is FALSE.
#'  
#' @importFrom foreach foreach '%do%'    
#' @importFrom methods new
#' @importFrom stats coef pnorm qnorm resid
#' @export
#' @return An object of class \code{mnr} containing the estimated regression 
#'   parameters, covariance matrix, the information matrix for regression
#'   parameters, and the residuals.
#' @examples
#' \dontrun{
#' # See rMNR or vignette for data generation
#' M = fit.mnr(Y=Y,X=X,eps=1e-8);
#' # Coefficients
#' coef(M);
#' # Outcome covariance
#' vcov(M,type="Outcome");
#' # Information matrix
#' vcov(M,type="Information");
#' # Residuals
#' resid(M);
#' }

fit.mnr = function(Y,X,sig=0.05,maxit=10,eps=1e-6,report=T){
  # Check input type
  if(!is.matrix(Y)){stop("A numeric matrix is expected for Y.")};
  if(!is.list(X)||(!all(unlist(lapply(X,is.matrix))))){
    stop("A list of numeric model matrices is expected for X.")};
  # Dimensions
  n = nrow(Y);
  k = ncol(Y);
  m = lapply(X,ncol);
  # Objective function
  Q = function(S){-n*log(det(S))};
  # Reusable matrices
  XtX = lapply(X,function(x){matIP(x,x)});
  XtXi = lapply(XtX,matInv);
  
  ## Initialize
  # List to hold initial parameters
  theta0 = list();
  # Initialize beta
  b0 = list();
  i = NULL;
  for(i in 1:k){
    b0[[i]] = MMP(XtXi[[i]],matIP(X[[i]],Y[,i]));
  }
  theta0$b = b0;
  # Residual matrix
  E0 = array(0,dim=c(n,0));
  for(i in 1:k){
    E0 = cbind(E0,Y[,i]-MMP(X[[i]],b0[[i]]));
  }
  theta0$E = E0;
  # Initialize sigma
  theta0$S = matIP(theta0$E,theta0$E)/n;
  # Precision
  theta0$L = matInv(theta0$S);
  
  # Update wrapper
  Update = function(theta){
    # Initial objective
    q0 = Q(theta$S);
    # Precision
    L = theta$L;
    # Residual matrix
    E = theta$E;
    # Update beta
    b1 = list();
    for(i in 1:k){
      wi = MMP(E[,-c(i),drop=F],L[-c(i),c(i),drop=F]/L[i,i]);
      b1[[i]] = b0[[i]]+MMP(XtXi[[i]],matIP(X[[i]],wi));
      # Update residual
      E[,i] = Y[,i]-MMP(X[[i]],b1[[i]]);
    }
    # Update covariance
    S1 = matIP(E,E)/n;
    L1 = matInv(S1);
    # New objective
    q1 = Q(S1);
    # Increment
    d = q1-q0;
    # Output
    Out = list("b"=b1,"E"=E,"S"=S1,"L"=L1,"d"=d);
    return(Out);
  }
  
  ## Maximzation
  for(i in 1:maxit){
    # Update
    theta1 = Update(theta0);
    # Accept if increment is positive
    if(theta1$d>0){
      theta0 = theta1;
      if(report){cat("Objective increment: ",signif(theta1$d,digits=3),"\n")}
    }
    # Terminate if increment is below tolerance
    if(theta1$d<eps){
      rm(theta1);
      break;
    }
  };
  
  ## Fitting report
  if(report){
    if(i<maxit){
      cat(paste0(i-1," update(s) performed before tolerance limit.\n\n"));
    } else {
      cat(paste0(i," update(s) performed without reaching tolerance limit.\n\n"));
    }
  };
  
  # Precision matrix
  L = matInv(theta0$S);
  ## Information for beta
  j = NULL;
  Ibb = foreach(i=1:k,.combine=rbind) %do% {
    Slice = foreach(j=1:k,.combine=cbind) %do% {
      if(j==i){
        Out = XtX[[i]]*L[i,i];
      } else {
        Out = matIP(X[[i]],X[[j]])*L[i,j];
      };
      return(Out);
    };
  };
  # Inverse 
  Ibbi = matInv(Ibb);
  
  ## Residuals
  E = theta0$E;
  colnames(E) = colnames(Y);
  rownames(E) = rownames(Y);
  
  ## Regression coefficients
  # Overall
  Point = unlist(theta0$b);
  # SE
  SE = sqrt(diag(Ibbi));
  # CIs
  z = qnorm(p=1-(sig/2));
  L = Point-z*SE;
  U = Point+z*SE;
  P = 2*pnorm(q=abs(Point/SE),lower.tail=F);
  # Labeling
  if(is.null(colnames(Y))){colnames(Y)=paste0("y",seq(1:k))};
  Lab = foreach(i=1:k,.combine=rbind) %do% {
    # Outcome
    Outcome = rep(colnames(Y)[i],times=m[[i]]);
    # Coefficient
    Coeff = colnames(X[[i]]);
    if(is.null(Coeff)){Coeff=paste0("x",i,seq(1:m[[i]]))};
    return(data.frame("Outcome"=Outcome,"Coeff"=Coeff));
  }
  Coeff = data.frame(Lab,"Point"=Point,"SE"=SE,"L"=L,"U"=U,"p"=P);
  rownames(Ibb) = colnames(Ibb) = Coeff$Coeff;
  # Covariance matrix
  S = theta0$S;
  rownames(S) = colnames(S) = colnames(Y);
  # Output
  Out = new(Class="mnr",Coefficients=Coeff,Covariance=S,Information=Ibb,Residuals=E);
  return(Out);
}
