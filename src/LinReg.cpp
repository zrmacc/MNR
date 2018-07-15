// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Univariate OLS model.
//' 
//' Fits the standard OLS model.
//' 
//' @param y Outcome.
//' @param Z Model matrix.
//' @export 
//' 
//' @return A list containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{Tau}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//' 
// [[Rcpp::export]]

SEXP fitNorm(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Z){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = Z.cols();
  // Gram matrix
  const Eigen::MatrixXd ZtZ = Z.transpose()*Z;
  // Estimate beta
  const Eigen::VectorXd b = (ZtZ).llt().solve(Z.transpose()*y);
  // Calculate residuals
  const Eigen::VectorXd eps = (y-Z*b);
  // Scale
  const double qf = (eps.transpose()*eps);
  const double tau = qf/(n-p);
  // Information
  const Eigen::MatrixXd Ibb = ZtZ/tau;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("Tau")=tau,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}