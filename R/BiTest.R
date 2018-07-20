# Purpose: Score tests for bivariate outcome regression
# Updated: 180713

#' Bivariate Regression Score Test
#' 
#' Score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the bivariate outcome 
#' model. The score test is specified using a logical vector \code{L}, with as
#' many entries as columns in the target model matrix \code{Zt}. The values of
#' \code{L} set to \code{T} are fixed at \eqn{\beta_{10}} under the null. The
#' values of \code{L} set to \code{F} are estimated.
#' 
#' @param yt Target outcome.
#' @param ys Surrogate outcome.
#' @param Zt Numeric model matrix for the target outcome.
#' @param Zs Numeric model matrix for the surrogate outcome.
#' @param L Logical vector, with as many entires as columns in the target 
#'   design, indicating which columns have fixed coefficients under the null.
#' @param b10 Value of the regression coefficient for the selected columns under 
#'   the null. Defaults to zero.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE. 
#'   
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value estimated using the chi-square distribution.
#' 
#' @examples
#' yt = D.bnr$yt;
#' ys = D.bnr$ys;
#' Zt = D.bnr$Zt;
#' Zs = D.bnr$Zs;
#' # Test for an overall effect, should reject.
#' Score.bnr(yt,ys,Zt,Zs,L=c(FALSE,TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}, should not reject.
#' Score.bnr(yt,ys,Zt,Zs,L=c(FALSE,FALSE,FALSE,TRUE));

Score.bnr = function(yt,ys,Zt,Zs,L,b10,maxit=100,eps=1e-6,report=F){
  # Input checks
  if(!is.vector(yt)){stop("A numeric vector is expected for yt.")};
  if(!is.vector(ys)){stop("A numeric vector is expected for ys.")};
  if(!is.matrix(Zt)){stop("A numeric matrix is expected for Zt.")};
  if(!is.matrix(Zs)){stop("A numeric matrix is expected for Zs.")};
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  if(length(L)!=ncol(Zt)){stop("L should have as many entries as columns in Zt.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")};
  if(sum(L)==length(L)){stop("At least 1 entry of L should be FALSE.")};
  # Check for missingness
  Miss = sum(is.na(yt))+sum(is.na(ys))+sum(is.na(Zt))+sum(is.na(Zs));
  if(Miss>0){stop("Inputs should contain no missing data.")};
  # Degrees of freedom
  df = sum(L);
  # Null coefficient
  if(missing(b10)){b10=rep(0,times=df)};
  
  # Partition target design
  # Zt1 is fixed under the null.
  # Zt2 is estimated under the null.
  Zt1 = Zt[,L,drop=F];
  Zt2 = Zt[,!L,drop=F];
  # Adjust response for fixed component
  yt = yt-as.numeric(fastMMp(Zt1,b10));
  # Fit null model
  M0 = fit.bnr(yt=yt,ys=ys,Zt=Zt2,Zs=Zs,maxit=maxit,eps=eps,report=report);
  # Extract precision
  Lambda = vcov(M0,type="Outcome",inv=T);
  # Partition precision
  LTT = Lambda[1,1];
  LTS = Lambda[1,2];
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = resid(M0,type="Surrogate");
  # Score vector
  u = LTT*eT+LTS*eS;
  Score = fastIP(Zt1,u);
  # Covariance matrix
  I11 = LTT*fastIP(Zt1,Zt1);
  I12 = cbind(LTT*fastIP(Zt1,Zt2),LTS*fastIP(Zt1,Zs));
  I22 = vcov(M0,type="Regression",inv=F);
  # Efficient information
  V = SchurC(I11=I11,I22=I22,I12=I12);
  # Score statistic
  Ts = fastQF(X=Score,A=fastInv(V));
  # P value
  p = pchisq(q=Ts,df=df,lower.tail=F);
  # Output
  Out = c("Score"=Ts,"df"=df,"p"=p);
  return(Out);
}

########################
# Repeated Score Test
########################

#' Repeated Bivariate Regression Score Test
#' 
#' Individually tests \eqn{H_{0}:\beta=0} for each column of \code{Zt1}, 
#' adjusting for \code{Zt2}. Missing values are permitted in \code{Zt1}, 
#' though not in the other model matrices. 
#' 
#' @param yt Target outcome vector.
#' @param ys Surrogate outcome vector.
#' @param Zt1 Numeric matrix of covariates for the target outcome whose 
#'   regression coefficients are zero under the null.
#' @param Zt2 Numeric matrix of covariates for the target outcome whose 
#'   regression coefficients are unconstrained, i.e. estimated, under the null.
#' @param Zs Numeric model matrix for the surrogate outcome.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE.  
#'   
#' @importFrom plyr aaply
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' @return A numeric vector of p-values, one for each column of \code{Zt1}. The
#'   p-values are estimated based on the chi-square distribution with one degree
#'   of freedom.
#' 
#' @examples
#' yt = D.bnr$yt;
#' ys = D.bnr$ys;
#' Zt2 = D.bnr$Zt;
#' Zs = D.bnr$Zs;
#' Zt1 = G;
#' # Test each column of G for association with the target outcome
#' R = rScore.bnr(yt,ys,Zt1,Zt2,Zs);

rScore.bnr = function(yt,ys,Zt1,Zt2,Zs,maxit=100,eps=1e-6,report=F){
  # Input checks
  if(!is.vector(yt)){stop("A numeric vector is expected for yt.")};
  if(!is.vector(ys)){stop("A numeric vector is expected for ys.")};
  if(!is.matrix(Zt1)){stop("A numeric matrix is expected for Zt1.")};
  if(!is.matrix(Zt2)){stop("A numeric matrix is expected for Zt2.")};
  if(!is.matrix(Zs)){stop("A numeric matrix is expected for Zs.")};
  # Check for missingness
  Miss = sum(is.na(yt))+sum(is.na(ys))+sum(is.na(Zt2))+sum(is.na(Zs));
  if(Miss>0){stop("Inputs other than Zt1 should contain no missing data.")};
  # Fit null model
  M0 = fit.bnr(yt=yt,ys=ys,Zt=Zt2,Zs=Zs,maxit=maxit,eps=eps,report=report);
  # Extract precision
  Lambda = vcov(M0,type="Outcome",inv=T);
  # Partition precision
  LTT = Lambda[1,1];
  LTS = Lambda[1,2];
  LSS = Lambda[2,2];
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = resid(M0,type="Surrogate");
  # Working vector
  z = LTT*eT+LTS*eS;
  # Function to calculate score statistics
  aux = function(x){
    # Check for missingness
    retain = !is.na(x);
    # If missingness, exclude 
    if(sum(!retain)>0){
      x.obs = x[retain];
      z.obs = z[retain];
      Zt2.obs = Zt2[retain,];
      Zs.obs = Zs[retain,];
    } else {
      x.obs = x;
      z.obs = z;
      Zt2.obs = Zt2;
      Zs.obs = Zs;
    }
    # Information matrices
    I11 = LTT*fastIP(x.obs,x.obs);
    I12 = cbind(LTT*fastIP(x.obs,Zt2.obs),LTS*fastIP(x.obs,Zs.obs));
    I22 = rbind(cbind(LTT*fastIP(Zt2.obs,Zt2.obs),LTS*fastIP(Zt2.obs,Zs.obs)),
                cbind(LTS*fastIP(Zs.obs,Zt2.obs),LSS*fastIP(Zs.obs,Zs.obs)));
    # Efficient information
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    # Score vector
    a = as.numeric(fastIP(x.obs,z.obs));
    # Score statistic
    Ts = (a^2)/V;
    return(Ts);
  }
  # Score statistics
  S = aaply(.data=Zt1,.margins=2,.fun=aux);
  # P values
  p = pchisq(q=S,df=1,lower.tail=F);
  return(p);
}