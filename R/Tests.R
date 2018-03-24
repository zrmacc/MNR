# Purpose: Score tests for multivariate outcome regression
# Updated: 180319

########################
# Construct Xi
########################

#' Construct Overall Surrogate Design
#' 
#' @param D.s List of \eqn{k-1} design matrices
#' @param parallel If true, assumes a parallel coefficient design, where each 
#'   secondary outcome has the same regression parameters.
#' @export

constrXi = function(D.s,parallel){
  # Number of surrogates
  m = length(D.s);
  flag = F;
  # Check each matrix has n observations
  obs = unlist(lapply(D.s,nrow));
  if(length(unique(obs))>1){
    flag = T;
    stop("Each design matrix must have n observations.");
  };
  n = unique(obs);
  # If parallel, check each matrix has q columns
  cols = unlist(lapply(D.s,ncol));
  if(parallel&(length(unique(cols))>1)){
    flag = T;
    stop("If parallel, each design matrix must have q columns.");
  }
  if(!parallel){
    # Total columns in Z.s
    Q = sum(cols);
    # Empty lists for design matrices and column names
    Z.s = Cn.s = list();
    # First matrix
    A = model.matrix(~0+.,data=data.frame(D.s[[1]]));
    Cn.s[[1]] = colnames(A);
    A = cbind(A,array(0,dim=c(n,Q-cols[1])));
    Z.s[[1]] = A;
    # Last matrix
    B = model.matrix(~0+.,data=data.frame(D.s[[m]]));
    Cn.s[[m]] = colnames(B)
    B = cbind(array(0,dim=c(n,Q-cols[m])),B);
    Z.s[[m]] = B;
    # Intervening matrices
    if(m>2){
      for(j in 2:(m-1)){
        # Left sum
        ls = sum(cols[1:(j-1)]);
        # Right sum
        rs = sum(cols[(j+1):m]);
        # Formatted matrix
        C = model.matrix(~0+.,data=data.frame(D.s[[j]]));
        Cn.s[[j]] = colnames(C);
        C = cbind(array(0,dim=c(n,ls)),C,array(0,dim=c(n,rs)));
        Z.s[[j]] = C;
      }
    }
  } else {
    Z.s = lapply(D.s,FUN=function(x){model.matrix(~0+.,data=data.frame(x))});
    Cn.s = list(colnames(Z.s[[1]]));
  }
  # Function to intercalate matrices
  aux = function(x){do.call(rbind,x)[order(sequence(sapply(x,nrow))),]};
  
  # Final matrix
  Z.s = aux(Z.s);
  Z.s = cbind(1,Z.s);
  colnames(Z.s) = c("(Intercept)",unlist(Cn.s));
  return(Z.s);
}

########################
# Bivariate case
########################

#' Bivariate Regression Wald Test
#' 
#' Wald test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the bivariate outcome model. 
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param Z.t Target model matrix. Include intercept and expand factors.
#' @param Z.s Surrogate model matrix. Include intercept and expand factors.
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
#' Z.t = model.matrix(~.,data=data.frame(D.bnr$D.t));
#' Z.s = model.matrix(~.,data=data.frame(D.bnr$D.s));
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
#' Score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} in the bivariate outcome model. 
#' 
#' @param y.t Target outcome.
#' @param y.s Surrogate outcome.
#' @param Z.t Target model matrix. Include intercept and expand factors.
#' @param Z.s Surrogate model matrix. Include intercept and expand factors.
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
#' Z.t = model.matrix(~.,data=data.frame(D.bnr$D.t));
#' Z.s = model.matrix(~.,data=data.frame(D.bnr$D.s));
#' # Test for an overall effect
#' # Note: Null model should contain the intercept
#' Score.bnr(y.t,y.s,Z.t,Z.s,L=c(FALSE,TRUE,TRUE,TRUE));
#' # Test for the effect of \eqn{\beta_{3}}
#' Score.bnr(y.t,y.s,Z.t,Z.s,L=c(FALSE,FALSE,FALSE,TRUE));

Score.bnr = function(y.t,y.s,Z.t,Z.s,L,b0,maxit=10,eps=1e-6){
  # Input checks
  if(is.matrix(y.s)){if(ncol(y.s)>1){stop("For multiple surrogates, use Score.mnr")}};
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(Z.t)){stop("L should have as many entries as columns in D.t.")};
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
  X2 = Z.t[,!L,drop=F];
  X1 = Z.t[,L,drop=F];
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
#' @param parallel If true, assumes a parallel coefficient design, where each 
#'   secondary outcome has the same regression parameters.
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

Score.mnr = function(y.t,y.s,D.t,D.s,L,b0,parallel=F,maxit=10,eps=1e-6){
  # Input checks
  if(is.vector(y.s)){stop("For a single surrogate, use Score.bnr")};
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(D.t)){stop("L should have as many entries as columns in D.t.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Check for missingness
  A = cbind(y.t,y.s,D.t,do.call(cbind,D.s));
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    D.t = D.t[keep,];
    D.s = lapply(D.s,FUN=function(x){x[keep,]});
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
  # Form Z.s
  Z.s = constrXi(D.s=D.s,parallel=parallel);
  # Fit null model
  M0 = fit.mnr(y.t=y.t,y.s=y.s,Z.t=X2,Z.s=Z.s,maxit=maxit,eps=eps,report=F);
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
  I12 = cbind(infoMNR(n=n,A=X1,B=X2,L=LTT),infoMNR(n=n,A=X1,B=Z.s,L=LTS));
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