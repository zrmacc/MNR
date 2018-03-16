// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for fitting bivariate outcome model
// Updated: 180315
#include <RcppEigen.h>

//' Bivariate Parameter Update
//' 
//' Parameter update for bivariate outcome model.
//' 
//' @param yt Target outcome
//' @param ys Surrogate outcome
//' @param Zt Target design
//' @param Zs Surrogate design
//' @param At Target inc. proj.
//' @param As Surrogate inc. proj.
//' @param b0 Current beta
//' @param a0 Current alpha
//' @param s0 Current sigma
// [[Rcpp::export]]

SEXP updateBVR(const Eigen::Map<Eigen::VectorXd> yt, const Eigen::Map<Eigen::VectorXd> ys,
               const Eigen::Map<Eigen::MatrixXd> Zt, const Eigen::Map<Eigen::MatrixXd> Zs,
               const Eigen::Map<Eigen::MatrixXd> At, const Eigen::Map<Eigen::MatrixXd> As,
               const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a0,
               const Eigen::Map<Eigen::MatrixXd> s0){
  // Observations
  const int n = Zt.rows();
  // Invert Sigma
  const Eigen::MatrixXd l0 = s0.completeOrthogonalDecomposition().pseudoInverse();
  // Surrogate residual
  const Eigen::VectorXd eS = (ys-Zs*a0);
  // Target working vector
  const double w1 = l0(0,1)/l0(0,0);
  const Eigen::VectorXd zT = (yt+w1*eS);
  // Update beta
  const Eigen::VectorXd b1 = At*zT;
  // Target residual
  const Eigen::VectorXd eT = (yt-Zt*b1);
  // Surrogate working vector
  const double w2 = l0(1,0)/l0(1,1);
  const Eigen::VectorXd zS = (w2*eT+ys);
  // Update alpha
  const Eigen::VectorXd a1 = As*zS;
  // Residual matrix
  Eigen::MatrixXd E(n,2);
  E.setZero();
  E.col(0) = (yt-Zt*b1);
  E.col(1) = (ys-Zs*a1); 
  // Covariance matrix
  const Eigen::MatrixXd ip = E.transpose()*E;
  const Eigen::MatrixXd s1 = (ip/n);
  // Output
  return Rcpp::List::create(Rcpp::Named("b1")=b1,Rcpp::Named("a1")=a1,Rcpp::Named("s1")=s1,Rcpp::Named("E1")=E);
}