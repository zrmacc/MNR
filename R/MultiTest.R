# Purpose: Score test for multivariate outcome regression
# Updated: 180714

########################
# Construct Xi
########################

#' Construct Overall Surrogate Design
#' 
#' Constructs the overall surrogate model matrix from a list of model matrices for the
#' individual surrogate outcomes.
#' 
#' @param Ls List of \eqn{(k-1)} numeric model matrices for the surrogate
#'   outcomes. See \code{\link{Score.mnr}}.
#' @param parallel If true, assumes a parallel coefficient design, where each 
#'   surrogate outcome has the same regression parameters.
#' @export
#' @return A numeric matrix of dimension \eqn{n(k-1)\times q}, where \eqn{n} is 
#' the number of observations, \eqn{k-1} is the number of surrogate outcomes,
#' and \eqn{q=q_{1}+\cdots+q_{(k-1)}} is the total number of covariates for
#' all surrogate outcomes.   

constrXi = function(Ls,parallel){
  # Number of surrogates
  m = length(Ls);
  flag = F;
  # Check each matrix has n observations
  obs = unlist(lapply(Ls,nrow));
  if(length(unique(obs))>1){
    flag = T;
    stop("Each design matrix must have n observations.");
  };
  n = unique(obs);
  # If parallel, check each matrix has q columns
  cols = unlist(lapply(Ls,ncol));
  if(parallel&(length(unique(cols))>1)){
    flag = T;
    stop("If parallel, each design matrix must have q columns.");
  };
  if(!parallel){
    # Total columns in Zs
    Q = sum(cols);
    # Empty lists for design matrices and column names
    Zs = Names = list();
    # First matri
    Names[[1]] = colnames(Ls[[1]]);
    A = cbind(Ls[[1]],array(0,dim=c(n,Q-cols[1])));
    Zs[[1]] = A;
    # Last matrix
    Names[[m]] = colnames(Ls[[m]]);
    B = cbind(array(0,dim=c(n,Q-cols[m])),Ls[[m]]);
    Zs[[m]] = B;
    # Intervening matrices
    if(m>2){
      for(j in 2:(m-1)){
        # Left sum
        ls = sum(cols[1:(j-1)]);
        # Right sum
        rs = sum(cols[(j+1):m]);
        # Formatted matrix
        Names[[j]] = colnames(Ls[[j]]);
        C = cbind(array(0,dim=c(n,ls)),Ls[[j]],array(0,dim=c(n,rs)));
        Zs[[j]] = C;
      };
    };
  } else {
    Zs = Ls;
    Names = list(colnames(Zs[[1]]));
  };
  # Function to intercalate matrices
  aux = function(x){do.call(rbind,x)[order(sequence(sapply(x,nrow))),]};
  # Final matrix
  Zs = aux(Zs);
  colnames(Zs) = c(unlist(Names));
  return(Zs);
};

########################
# Multivariate case
########################

#' Multivariate Regression Score Test
#' 
#' Score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the multivariate outcome
#' model. The score test is specified using a logical vector \code{L}, with as
#' many entries as columns in the target model matrix. The values of \code{L} 
#' set to \code{T} are fixed at \eqn{\beta_{10}} under the null. The values of 
#' \code{L} set to \code{F} are estimated.
#' 
#' @section Useage notes:
#' \itemize{
#'  \item \code{Ls} is a list of model matrices, one for each surrogate outcome.
#'  Each matrix should the same number of rows, but the number of columns may
#'  vary. The overall design matrix is constructed internally. By default, each
#'  surrogate outcome is allocated its own regression coefficients. If a
#'  parallel design is specified, then all surrogate outcomes share the same set
#'  of regression coefficients. In this case, each model matrix in \code{Ls}
#'  must have the same number of columns as well.
#' }
#' 
#' @param yt Target outcome vector.
#' @param Ys Surrogate outcome matrix.
#' @param Zt Numeric model matrix for the target outcome.
#' @param Ls List of numeric model matrices for the surrogate outcomes. 
#' @param L Logical vector, with as many entires as columns in the target
#'   design, indicating which columns design are fixed under the null.
#' @param b10 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param parallel If true, assumes a parallel coefficient design, where each 
#'   secondary outcome has the same regression parameters.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param REML Apply REML correction to covariance matrix? Default is FALSE. 
#' @param report Report model fitting progress? Default is FALSE. 
#' 
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#' 
#' @examples
#' yt = D.mnr$yt;
#' Ys = D.mnr$Ys;
#' Zt = D.mnr$Zt;
#' Ls = D.mnr$Ls;
#' # Test for an overall effect, should reject.
#' Score.mnr(yt,Ys,Zt,Ls,L=c(FALSE,TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}, should not reject.
#' Score.mnr(yt,Ys,Zt,Ls,L=c(FALSE,FALSE,FALSE,TRUE));

Score.mnr = function(yt,Ys,Zt,Ls,L,b10,parallel=F,maxit=100,eps=1e-8,REML=F,report=F){
  # Check input type
  if(!is.vector(yt)){stop("A numeric vector is expected for yt.")};
  if(!is.matrix(Ys)){stop("A numeric vector is expected for Ys.")};
  if(!is.matrix(Zt)){stop("A numeric matrix is expected for Zt.")};
  if(!all(unlist(lapply(Ls,is.matrix)))){stop("A list of numeric matrices is expected for Ls.")};
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  
  # Check test specification
  if(length(L)!=ncol(Zt)){stop("L should have as many entries as columns in Zt.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")};
  if(sum(L)==length(L)){stop("At least 1 entry of L should be FALSE.")};
  
  # Check for missingness
  Miss = sum(is.na(yt))+sum(is.na(Ys))+sum(is.na(Zt))+sum(is.na(Ls));
  if(Miss>0){stop("Inputs should contain no missing data.")};
  # Degrees of freedom
  df = sum(L);
  # Null coefficient
  if(missing(b10)){b10=rep(0,times=df)};
  
  # Observations
  n = length(yt);
  # Partition target design
  # Zt1 is fixed under the null.
  # Zt2 is estimated under the null.
  Zt1 = Zt[,L,drop=F];
  Zt2 = Zt[,!L,drop=F];
  # Form surrogate design 
  Zs = constrXi(Ls=Ls,parallel=parallel);
  # Regression parameters estimated
  k = ncol(Zs)+ncol(Zt2);
  # Adjust response for fixed component
  yt = yt-as.numeric(fastMMp(Zt1,b10));
  # Fit null model
  M0 = fit.mnr(yt=yt,Ys=Ys,Zt=Zt2,Zs=Zs,maxit=maxit,eps=eps,report=report);
  # Extract covariance
  Sigma = vcov(M0,type="Outcome",inv=F);
  # REML-type adjustment
  if(REML){
    diag(Sigma) = n/(n-k)*diag(Sigma);
  };
  # Precision matrix
  Lambda = fastInv(Sigma);
  # Partition precision
  LTT = Lambda[1,1];
  LTS = Lambda[1,2:ncol(Lambda),drop=F];
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = c(fastT(resid(M0,type="Surrogate")));
  # Score vector
  n = length(eT);
  Score = infoMNR(n=n,A=Zt1,B=eT,L=LTT) + infoMNR(n=n,A=Zt1,B=eS,L=LTS);
  # Covariance matrix
  I11 = infoMNR(n=n,A=Zt1,B=Zt1,L=LTT);
  I12 = cbind(infoMNR(n=n,A=Zt1,B=Zt2,L=LTT),infoMNR(n=n,A=Zt1,B=Zs,L=LTS));
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