// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for fitting bivariate outcome model
// Updated: 180315
#include <RcppEigen.h>

//' Bivariate Parameter Update
//' 
//' Parameter update for bivariate outcome model.
//' 
//' @param t Target outcome
//' @param s Surrogate outcome
//' @param X Target design
//' @param Xi Surrogate design
//' @param B Target inc. proj.
//' @param A Surrogate inc. proj.
//' @param b0 Current beta
//' @param a0 Current alpha
//' @param s0 Current sigma
//' @export 
// [[Rcpp::export]]

SEXP updateBVR(const Eigen::Map<Eigen::VectorXd> t, const Eigen::Map<Eigen::VectorXd> s,
               const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Xi,
               const Eigen::Map<Eigen::MatrixXd> B, const Eigen::Map<Eigen::MatrixXd> A,
               const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a0,
               const Eigen::Map<Eigen::MatrixXd> s0){
  // Observations
  const int n = X.rows();
  // Invert Sigma
  const Eigen::MatrixXd l0 = s0.completeOrthogonalDecomposition().pseudoInverse();
  // Surrogate residual
  const Eigen::VectorXd es = (s-Xi*a0);
  // Target working vector
  const double ws = l0(0,1)/l0(0,0);
  const Eigen::VectorXd zt = (t+ws*es);
  // Update beta
  const Eigen::VectorXd b1 = B*zt;
  // Target residual
  const Eigen::VectorXd et = (t-X*b1);
  // Surrogate working vector
  const double wt = l0(1,0)/l0(1,1);
  const Eigen::VectorXd zs = (wt*et+s);
  // Update alpha
  const Eigen::VectorXd a1 = A*zs;
  // Residual matrix
  Eigen::MatrixXd E(n,2);
  E.setZero();
  E.col(0) = (t-X*b1);
  E.col(1) = (s-Xi*a1); 
  // Covariance matrix
  const Eigen::MatrixXd ip = E.transpose()*E;
  const Eigen::MatrixXd s1 = (ip/n);
  // Output
  return Rcpp::List::create(Rcpp::Named("b1")=b1,Rcpp::Named("a1")=a1,Rcpp::Named("s1")=s1,Rcpp::Named("E1")=E);
}