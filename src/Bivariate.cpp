// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Bivariate Normal, Update for Beta
//' 
//' Parameter update \eqn{\beta} in the bivariate outcome model.
//' 
//' @param s Surrogate outcome
//' @param Zt Target design
//' @param B Target design inner product.
//' @param Zs Surrogate design
//' @param b0 Initial beta
//' @param a Current alpha
//' @param L Current precision.
//' @export 
// [[Rcpp::export]]

SEXP BNRbeta(const Eigen::Map<Eigen::VectorXd> s, const Eigen::Map<Eigen::MatrixXd> Zt,
             const Eigen::Map<Eigen::MatrixXd> B, const Eigen::Map<Eigen::MatrixXd> Zs,
             const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a,
             const Eigen::Map<Eigen::MatrixXd> L){
  // Surrogate residuals
  const Eigen::VectorXd es = (s-Zs*a);
  // Weight
  const double w = L(0,1)/L(0,0);
  // Weighted residuals
  const Eigen::VectorXd wes = w*es;
  // Update
  const Eigen::MatrixXd Ztwes = Zt.transpose()*wes;
  const Eigen::VectorXd b1 = b0 + B.llt().solve(Ztwes);
  return Rcpp::wrap(b1);
}

//' Bivariate Normal, Update for Alpha
//' 
//' Parameter update for \eqn{\alpha} in the bivariate outcome model.
//' 
//' @param t Target outcome
//' @param Zt Target design
//' @param Zs Surrogate design
//' @param A Surrogate design inner product
//' @param b Current beta
//' @param a0 Initial alpha
//' @param L Current precision
//' @export 
// [[Rcpp::export]]

SEXP BNRalpha(const Eigen::Map<Eigen::VectorXd> t, const Eigen::Map<Eigen::MatrixXd> Zt,
              const Eigen::Map<Eigen::MatrixXd> Zs, const Eigen::Map<Eigen::MatrixXd> A,
              const Eigen::Map<Eigen::VectorXd> b, const Eigen::Map<Eigen::VectorXd> a0,
              const Eigen::Map<Eigen::MatrixXd> L){
  // Target residuals
  const Eigen::VectorXd et = (t-Zt*b);
  // Weight
  const double w = L(0,1)/L(1,1);
  // Weighted residuals
  const Eigen::VectorXd wet = w*et;
  // Update
  const Eigen::MatrixXd Zswet = Zs.transpose()*wet;
  const Eigen::VectorXd a1 = a0 + A.llt().solve(Zswet);
  return Rcpp::wrap(a1);
}