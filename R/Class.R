#' Multivariate Regression Model
#'
#' Defines the object class returned by\code{\link{fit.mnr}}. 
#'
#' @slot Coefficients Regression coefficients.
#' @slot Covariance Outcome covariance matrix.
#' @slot Information Information for regression coefficients. 
#' @slot Residuals Outcome residuals. 
#' @name mnr-class
#' @rdname mnr-class
#' @exportClass mnr

setClass(Class="mnr",representation=representation(Coefficients="data.frame",Covariance="matrix",Information="matrix",Residuals="matrix"));

########################
# Print Method
########################

#' Print for Multivariate Regression Model
#' 
#' Print method for objects of class \code{mnr}. 
#'
#' @param x A \code{mnr} object.
#' @param ... Unused.
#' @export

print.mnr = function(x,...){
  Coeff = x@Coefficients;
  aux = function(v){
    if(is.numeric(v)){return(signif(v,digits=3))}
    else{return(v)};
  };
  Coeff[] = lapply(Coeff,aux);
  print(Coeff);
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

#' Extract Coefficients from Multivariate Regression Model
#' 
#' Returns the estimated regression coefficients from a fitted \code{mnr} model. 
#'
#' @param object A \code{mnr} object.
#' @param ... Unused.
#' @param outcome Select an outcome from among 
#' @export
#' @return A data.frame containing estimated coefficients

coef.mnr = function(object,...,outcome=NULL){
  # Coefficient frame
  Coeff = object@Coefficients;
  if(is.null(outcome)){
    return(Coeff);
  } else {
    Choices = unique(Coeff$Outcome);
    if(!(outcome %in% Choices)){
      stop(paste("Select outcome from among:",paste(Choices,collapse=" ")));
    } else {
      Out = Coeff[Coeff$Outcome==outcome,];
      return(Out);
    }
  }
};

########################
# Resid Method
########################

#' Extract Residuals from Multivariate Regression Model
#' 
#' Returns the estimated residuals from a fitted \code{mnr} model.
#'
#' @param object A \code{mnr} object.
#' @param ... Unused.
#' @export
#' @return A numeric matrix of residuals. 

residuals.mnr = function(object,...){
  return(object@Residuals);
}

########################
# Vcov Method
########################

#' Extract Covariance Matrix from Multivariate Regression Model
#'
#' Returns estimated covariance matrices from a fitted \code{mnr} model. Specify
#' "Regression" for the information matrix of the regression parameters. Specify
#' "Outcome" for the outcome covariance matrix. 
#' 
#' @param object A \code{mnr} object.
#' @param ... Unused.
#' @param type Either "Information" or "Outcome". Default is "Information".  
#' @param inv Invert information matrix? Default is FALSE.
#' @export
#' @return A numeric matrix. 

vcov.mnr = function(object,...,type="Information",inv=F){
  if(! type %in% c("Information","Outcome")){stop("Select Information or Outcome.")};
  if(type=="Outcome"){
    Out = object@Covariance;
    if(inv){Out = fastInv(Out)};
    return(Out);
  } else {
    # Beta information
    Out = object@Information;
    # Invert
    if(inv){Out = fastInv(Out);}
    return(Out);
  }
};
