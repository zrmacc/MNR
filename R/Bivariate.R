# Purpose: Score tests for bivariate outcome regression
# Updated: 180409

#' Bivariate Regression Wald Test
#' 
#' Wald test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the bivariate outcome model. 
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param Z.t Numeric model matrix for the target outcome.
#' @param Z.s Numeric model matrix for the surrogate outcome.
#' @param L Logical vector, with as many entires as columns in the target
#'   model, indicating which columns of the target model are fixed under the
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
#' Z.t = D.bnr$Z.t;
#' Z.s = D.bnr$Z.s;
#' # Test for an overall effect
#' # Note: Null model should contain the intercept
#' Wald.bnr(y.t,y.s,Z.t,Z.s,L=c(FALSE,TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}
#' Wald.bnr(y.t,y.s,Z.t,Z.s,L=c(FALSE,FALSE,FALSE,TRUE));

Wald.bnr = function(y.t,y.s,Z.t,Z.s,L,b0,maxit=10,eps=1e-6){
  # Input checks
  if(is.matrix(y.s)){if(ncol(y.s)>1){stop("For multiple surrogates, use Score.mnr")}};
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(Z.t)){stop("L should have as many entries as columns in Z.t.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(sum(L)==length(L)){stop("At least 1 entry of L should be FALSE.")};
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Check for missingness
  A = cbind(y.t,y.s,Z.t,Z.s);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    Z.t = Z.t[keep,];
    Z.s = Z.s[keep,];
  };
  # Partition target design
  X1 = Z.t[,L,drop=F];
  X2 = Z.t[,!L,drop=F];
  Xi = Z.s;
  # Fit full model
  M0 = fit.bnr(y.t=y.t,y.s=y.s,Z.t=Z.t,Z.s=Xi,maxit=maxit,eps=eps,report=F);
  # Efficient information
  V = effInfo(L,M0);
  # Wald statistic
  b1 = coef(M0,type="Beta")[L];
  Tw = fastQF(x=(b1-b0),A=V);
  # Degrees of freedom
  df = ncol(X1);
  # P value
  p = pchisq(q=Tw,df=df,lower.tail=F);
  # Output
  Out = c("Wald"=Tw,"df"=df,"p"=p);
  return(Out);
}

#' Bivariate Regression Score Test
#' 
#' Score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the bivariate outcome
#' model.
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param Z.t Numeric model matrix for the target outcome.
#' @param Z.s Numeric model matrix for the surrogate outcome.
#' @param L Logical vector, with as many entires as columns in the target 
#'   design, indicating which columns have fixed coefficients under the null.
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
#' Z.t = D.bnr$Z.t;
#' Z.s = D.bnr$Z.s;
#' # Test for an overall effect
#' # Note: Null model should contain the intercept
#' Score.bnr(y.t,y.s,Z.t,Z.s,L=c(FALSE,TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}
#' Score.bnr(y.t,y.s,Z.t,Z.s,L=c(FALSE,FALSE,FALSE,TRUE));

Score.bnr = function(y.t,y.s,Z.t,Z.s,L,b0,maxit=10,eps=1e-6){
  # Input checks
  if(is.matrix(y.s)){if(ncol(y.s)>1){stop("For multiple surrogates, use Score.mnr")}};
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(Z.t)){stop("L should have as many entries as columns in Z.t.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(sum(L)==length(L)){stop("At least 1 entry of L should be FALSE.")};
  # Check for missingness
  A = cbind(y.t,y.s,Z.t,Z.s);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    Z.t = Z.t[keep,];
    Z.s = Z.s[keep,];
  };
  # Null coefficient
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Partition target design
  X1 = Z.t[,L,drop=F];
  X2 = Z.t[,!L,drop=F];
  Xi = Z.s;
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
# Repeated Score Test
########################

#' Bivariate Regression Score Test
#' 
#' Individually tests \eqn{H_{0}:\beta_{j}=0} for each column of \code{X1}, 
#' adjusting for \code{X2}.
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param X1 Numeric matrix of covariates for the target outcome whose
#'   regression coefficients are constrained under the null.
#' @param X2 Numeric matrix of covariates for the target outcome whose
#'   regression coefficients are unconstrained under the null.
#' @param Z.s Numeric model matrix for the surrogate outcome.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#'   
#' @importFrom plyr aaply
#' @importFrom stats model.matrix pchisq resid vcov
#' @export

rScore.bnr = function(y.t,y.s,X1,X2,Z.s,maxit=10,eps=1e-6){
  # Input checks
  if(is.matrix(y.s)){if(ncol(y.s)>1){stop("For multiple surrogates, use Score.mnr")}};
  # Check for missingness
  A = cbind(y.t,y.s,X1,X2,Z.s);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    Z.t = Z.t[keep,];
    Z.s = Z.s[keep,];
  };
  # Fit null model
  M0 = fit.bnr(y.t=y.t,y.s=y.s,Z.t=X2,Z.s=Z.s,maxit=maxit,eps=eps,report=F);
  # Extract precision
  Lambda = vcov(M0,type="Outcome",inv=T);
  # Partition precision
  LTT = Lambda[1,1];
  LTS = Lambda[1,2];
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = resid(M0,type="Surrogate");
  # Fixed information component
  I22 = vcov(M0,type="Regression",inv=F);
  # Function to calculate score statistics
  aux = function(x){
    # Score vector
    Score = fastIP(x,LTT*eT+LTS*eS);
    # Test-dependent information components
    I11 = LTT*fastIP(x,x);
    I12 = cbind(LTT*fastIP(x,X2),LTS*fastIP(x,Z.s));
    # Efficient information
    V = SchurC(I11=I11,I22=I22,I12=I12);
    # Score statistic
    Ts = fastQF(x=Score,A=fastInv(V));
    return(Ts);
  }
  # Score statistics
  S = aaply(.data=X1,.margins=2,.fun=aux);
  # P values
  p = pchisq(q=S,df=1,lower.tail=F);
  return(p);
}