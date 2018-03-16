#' Multivariate Regression Model
#'
#' @slot Coefficients Regression coefficients.
#' @slot Covariance Outcome covariance matrix.
#' @slot Information Information.
#' @slot Residuals Phenotypic residuals.
#' @name mnr-class
#' @rdname mnr-class
#' @exportClass mnr
setClass(Class="mnr",representation=representation(Coefficients="list",Covariance="matrix",Information="list",Residuals="list"));

########################
# Print Method
########################

#' Print for Multivariate Regression Model
#'
#' @param x A \code{mnr} object.
#' @param ... Unused.
#' @export

print.mnr = function(x,...){
  cat("Target Outcome Regression Coefficients:\n");
  print(signif(x@Coefficients$Beta,digits=3));
  cat("\n");
  cat("Surrogate Outcome Regression Coefficients:\n");
  print(signif(x@Coefficients$Alpha,digits=3));
}

########################
# Show Method
########################

#' Show for Multivariate Regression Model
#' @param object A \code{mnr} object.
#' @rdname mnr-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="mnr"),definition=function(object){print.mnr(x=object)});

########################
# Coef Method
########################

#' Convert to proper case
#'
#' @param s A string.

toProper = function(s){
  a = unlist(strsplit(s,split=""));
  a[1] = toupper(a[1]);
  a[2:length(a)] = tolower(a[2:length(a)]);
  return(paste0(a,collapse=""));
}

#' Extract Coefficients from Multivariate Regression Model
#'
#' @param object A \code{mnr} object.
#' @param ... Unused.
#' @param type Either Alpha or Beta.
#' @export

coef.mnr = function(object,...,type="Beta"){
  # Ensure proper case
  type = toProper(type);
  if(! type %in% c("Alpha","Beta")){stop("Select Alpha or Beta.")};
  if(type=="Alpha"){
    return(object@Coefficients[["Alpha"]]);
  } else {
    return(object@Coefficients[["Beta"]]);
  }
}

########################
# Resid Method
########################

#' Extract Residuals from Multivariate Regression Model
#'
#' @param object A \code{mnr} object.
#' @param ... Unused.
#' @param type Either Surrogate or Target.
#' @export

residuals.mnr = function(object,...,type="Target"){
  # Ensure first letter is capitalized
  type = toProper(type);
  if(! type %in% c("Surrogate","Target")){stop("Select Surrogate or Target.")};
  if(type=="Target"){
    return(object@Residuals[["Target"]]);
  } else {
    return(object@Residuals[["Surrogate"]]);
  }
}

########################
# Vcov Method
########################

#' Extract Covariance Matrix from Multivariate Regression Model
#'
#' Returns the either the estimated covariance matrix of the outcome,
#' or the information matrix of the regression coefficients. 
#' 
#' @param object A \code{mnr} object.
#' @param ... Unused.
#' @param type Either Outcome or Regression. Default is Regression.  
#' @param inv Invert information matrix? Default is TRUE.
#' @export

vcov.mnr = function(object,...,type="Regression",inv=T){
  # Ensure proper case
  type = toProper(type);
  if(! type %in% c("Regression","Outcome")){stop("Select Regression or Outcome.")};
  if(type=="Outcome"){
    Out = object@Covariance;
    if(inv){Out = fastInv(Out)};
    return(Out);
  } else {
    # Beta information
    Ibb = object@Information$Ibb;
    # Alpha information
    Iaa = object@Information$Iaa;
    # Cross inforation
    Iba = object@Information$Iba;
    # Binding
    A = cbind(Ibb,Iba);
    B = cbind(fastInv(Iba),Iaa);
    Out = rbind(A,B);
    # Invert
    if(inv){Out = fastInv(Out);}
    return(Out);
  }
};
