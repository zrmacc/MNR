# Purpose: Score tests for multivariate outcome regression
# Updated: 180315

########################
# Bivariate case
########################

#' Bivariate Regression Score Test
#' 
#' Score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the bivariate outcome model. 
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param D.t Target design matrix.
#' @param D.s Surrogate design matrix.
#' @param L Logical vector, with as many entires as columns in the target
#'   design, indicating which columns of the target design are fixed under the
#'   null.
#' @param b0 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' 
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' 
#' @examples
#' y.t = D.bnr$y.t;
#' y.s = D.bnr$y.s;
#' D.t = D.bnr$D.t;
#' D.s = D.bnr$D.s;
#' # Test for an overall effect
#' Score.bnr(y.t,y.s,D.t,D.s,L=c(TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}
#' Score.bnr(y.t,y.s,D.t,D.s,L=c(FALSE,FALSE,TRUE));

Score.bnr = function(y.t,y.s,D.t,D.s,L,b0,maxit=10,eps=1e-6){
  # Input checks
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(D.t)){stop("L should have as many entries as columns in D.t.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Check for missingness
  A = cbind(y.t,y.s,D.t,D.s);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    D.t = D.t[keep,];
    D.s = D.s[keep,];
    };
  # Partition target design
  D.t.reduced = D.t[,!L,drop=F];
  D.t.null = D.t[,L,drop=F];
  # Convert to model matrices
  X1 = model.matrix(~0+.,data=data.frame(D.t.null));
  if(ncol(D.t.reduced)==0){
    n = nrow(D.t.null);
    X2 = matrix(rep(1,n),ncol=1);
  } else {
    X2 = model.matrix(~.,data=data.frame(D.t.reduced));
  };
  Xi = model.matrix(~.,data=data.frame(D.s));
  # Fit null model
  M0 = fit.bnr(y.t=y.t,y.s=y.s,Z.t=X2,Z.s=Xi,maxit=maxit,eps=eps,report=F);
  # Extract precision
  Lambda = vcov(M0,type="Outcome",inv=T);
  # Partition precision
  LTT = Lambda[1,1];
  LTS = Lambda[1,2];
  # Extract residuals
  eT = resid(M0,type="Target") - as.numeric(X1%*%b0);
  eS = resid(M0,type="Surrogate");
  # Score vector
  Score = fastT(X1)%*%(LTT*eT+LTS*eS);
  # Covariance matrix
  I11 = LTT*fastIP(X1,X1);
  I12 = cbind(LTT*fastIP(X1,X2),LTS*fastIP(X1,Xi));
  I22 = vcov(M0,type="Regression",inv=F);
  # Efficient information
  V = SchurC(I11=I11,I22=I22,I12=I12);
  # Score statistic
  Ts = fastQF(x=Score,A=fastInv(V));
  # Degrees of freedom
  df = ncol(X1);
  # P value
  p = pchisq(q=Ts,df=df,lower.tail=F);
  # Output
  Out = c("Score"=Ts,"df"=df,"p"=p);
  return(Out);
}

########################
# Multivariate case
########################

#' Multivariate Regression Score Test
#' 
#' Score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the multivariate outcome model. 
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param D.t Target design matrix.
#' @param D.s Surrogate design matrix.
#' @param L Logical vector, with as many entires as columns in the target
#'   design, indicating which columns of the target design are fixed under the
#'   null.
#' @param b0 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' 
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' 
#' @examples
#' y.t = D.mnr$y.t;
#' y.s = D.mnr$y.s;
#' D.t = D.mnr$D.t;
#' D.s = D.mnr$D.s;
#' # Test for an overall effect
#' Score.mnr(y.t,y.s,D.t,D.s,L=c(TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}
#' Score.mnr(y.t,y.s,D.t,D.s,L=c(FALSE,FALSE,TRUE));

Score.mnr = function(y.t,y.s,D.t,D.s,L,b0,maxit=10,eps=1e-6){
  # Input checks
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(D.t)){stop("L should have as many entries as columns in D.t.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Check for missingness
  A = cbind(y.t,y.s,D.t);
  aux = function(x){sum(is.na(x))>0};
  keep1 = !apply(A,MARGIN=1,FUN=aux);
  keep2 = !unlist(lapply(X=D.s,FUN=aux));
  keep = (keep1&keep2);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    D.t = D.t[keep,];
    D.s = D.s[keep,];
  };
  # Partition target design
  D.t.reduced = D.t[,!L,drop=F];
  D.t.null = D.t[,L,drop=F];
  # Convert to model matrices
  X1 = model.matrix(~0+.,data=data.frame(D.t.null));
  if(ncol(D.t.reduced)==0){
    n = nrow(D.t.null);
    X2 = matrix(rep(1,n),ncol=1);
  } else {
    X2 = model.matrix(~.,data=data.frame(D.t.reduced));
  };
  Xi = do.call(rbind,D.s);
  Xi = model.matrix(~.,data=data.frame(Xi));
  # Fit null model
  M0 = fit.mnr(y.t=y.t,y.s=y.s,Z.t=X2,Z.s=Xi,maxit=maxit,eps=eps,report=F);
  # Extract precision
  Lambda = vcov(M0,type="Outcome",inv=T);
  # Partition precision
  LTT = Lambda[1,1];
  LTS = Lambda[1,2:ncol(Lambda),drop=F];
  # Extract residuals
  eT = resid(M0,type="Target") - as.numeric(X1%*%b0);
  eS = c(fastT(resid(M0,type="Surrogate")));
  # Score vector
  n = length(eT);
  Score = infoMNR(n=n,A=X1,B=eT,L=LTT) + infoMNR(n=n,A=X1,B=eS,L=LTS);
  # Covariance matrix
  I11 = infoMNR(n=n,A=X1,B=X1,L=LTT);
  I12 = cbind(infoMNR(n=n,A=X1,B=X2,L=LTT),infoMNR(n=n,A=X1,B=Xi,L=LTS));
  I22 = vcov(M0,type="Regression",inv=F);
  # Efficient information
  V = SchurC(I11=I11,I22=I22,I12=I12);
  # Score statistic
  Ts = fastQF(x=Score,A=fastInv(V));
  # Degrees of freedom
  df = ncol(X1);
  # P value
  p = pchisq(q=Ts,df=df,lower.tail=F);
  # Output
  Out = c("Score"=Ts,"df"=df,"p"=p);
  return(Out);
}