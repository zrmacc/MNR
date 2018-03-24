# Purpose: Matrix operations
# Updated: 180324

#' Efficient Information
#' 
#' Calculate efficient information for a subset of components if \eqn{\beta}.
#' 
#' @param L Logical vector indicating the components of interest
#' @param M Fitted model

effInfo = function(L,M){
  # Input checks
  if(class(M)!="mnr"){stop("M is the output of fit.bnr or fit.mnr.")};
  if(!is.logical(L)){stop("L should be a logical vector.")};
  n.b = length(M@Coefficients$Beta);
  if(length(L)!=n.b){stop("L should the same length as beta.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")};
  # Complete case
  if(sum(L)==length(L)){
    Out = SchurC(I11=M@Information$Ibb,I22=M@Information$Iaa,I12=M@Information$Iba);
  } else {
    I11 = M@Information$Ibb[L,L,drop=F];
    I12 = cbind(M@Information$Ibb[L,!L,drop=F],M@Information$Iba[L,,drop=F]);
    I22.a = cbind(M@Information$Ibb[!L,!L,drop=F],M@Information$Iba[!L,,drop=F]);
    I22.b = cbind(fastT(M@Information$Iba[!L,,drop=F]),M@Information$Iaa);
    I22 = rbind(I22.a,I22.b);
    Out = SchurC(I11=I11,I22=I22,I12=I12);
  }
  # Output
  return(Out);
} 