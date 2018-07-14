// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for fitting multivariate outcome model
// Updated: 180404
#include <RcppEigen.h>

//' Initialize Alpha
//' 
//' Initialize \eqn{\alpha} for multivariate outcome model.
//' 
//' @param Zs Surrogate design matrix
//' @param s Surrogate outcome
//' @export 
// [[Rcpp::export]]

SEXP alpha0(const Eigen::MatrixXd Zs, const Eigen::VectorXd s){
  const Eigen::MatrixXd A = (Zs.transpose()*Zs);
  const Eigen::VectorXd a = A.llt().solve(Zs.transpose()*s);
  return Rcpp::wrap(a);
}

//' Multivariate Normal, Update for Beta
//' 
//' Parameter update for multivariate outcome model.
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

SEXP MNRbeta(const Eigen::Map<Eigen::VectorXd> s, const Eigen::Map<Eigen::MatrixXd> Zt,
             const Eigen::Map<Eigen::MatrixXd> B, const Eigen::Map<Eigen::MatrixXd> Zs,
             const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a,
             const Eigen::Map<Eigen::MatrixXd> L){
  // Surrogate residuals
  const Eigen::VectorXd es = (s-Zs*a);
  // Dimensions
  const int n = Zt.rows();
  const int m = L.cols()-1;
  // Partition Lambda
  const double Ltt = L(0,0);
  const Eigen::VectorXd Lst = L.block(1,0,m,1);
  // Weight
  const Eigen::RowVectorXd w = (Lst.transpose())/Ltt;
  // Calculate weighted surrogate residuals
  Eigen::VectorXd wes(n);
  wes.setZero();
  for(int i=0; i<n; i++){
    wes(i) = (w*es.segment(i*m,m));
  }
  // Update
  const Eigen::VectorXd b1 = b0 + B.llt().solve(Zt.transpose()*wes);
  return Rcpp::wrap(b1);
}

//' Multivariate Normal, Update for Alpha
//' 
//' Parameter update for multivariate outcome model.
//' 
//' @param t Target outcome
//' @param s Surrogate outcome
//' @param Zt Target design
//' @param Zs Surrogate design
//' @param b Current beta
//' @param a0 Initial alpha
//' @param L Current precision
//' @export 
// [[Rcpp::export]]

SEXP MNRalpha(const Eigen::Map<Eigen::VectorXd> t, const Eigen::Map<Eigen::VectorXd> s,
              const Eigen::Map<Eigen::MatrixXd> Zt, const Eigen::Map<Eigen::MatrixXd> Zs, 
              const Eigen::Map<Eigen::VectorXd> b, const Eigen::Map<Eigen::MatrixXd> L){
  // Target residuals
  const Eigen::VectorXd et = (t-Zt*b);
  // Dimensions
  const int n = Zt.rows();
  const int m = L.cols()-1;
  const int q = Zs.cols();
  // Partition Lambda
  const Eigen::VectorXd Lst = L.block(1,0,m,1);
  const Eigen::MatrixXd Lss = L.block(1,1,m,m);
  // Create matrix A = Zs'Lss Zs;
  Eigen::MatrixXd A(q,q);
  A.setZero();
  for(int i=0; i<n; i++){
    A += (Zs.block(i*m,0,m,q).transpose())*Lss*Zs.block(i*m,0,m,q);
  };
  // Surrogate working vector
  Eigen::VectorXd zs(n*m);
  zs.setZero();
  // Create (Lst*et+Lss*s)
  for(int i=0; i<n; i++){
    zs.segment(i*m,m) = (Lst*et(i))+(Lss*s.segment(i*m,m));
  };
  // Update
  const Eigen::VectorXd a1 = A.llt().solve(Zs.transpose()*zs);
  // Output
  return Rcpp::wrap(a1);
}

//' Information
//'
//' Calculate \eqn{A'(I_{n\times n}\otimes L)B}, the general form of information matrices
//' for the multivariate outcome model. 
//'
//' @param n Observations
//' @param A Left hand matrix
//' @param B Right hand matrix
//' @param L Precision component
//' @export 
// [[Rcpp::export]]

SEXP infoMNR(const int n, const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B, 
             const Eigen::Map<Eigen::MatrixXd> L){
  // Dimensions
  const int a = L.rows();
  const int b = L.cols();
  const int p = A.cols();
  const int q = B.cols();
  // Output structure
  Eigen::MatrixXd Out(p,q);
  Out.setZero();
  // Loop over observations
  for(int i=0; i<n; i++){
    Out += (A.block(i*a,0,a,p).transpose())*L*(B.block(i*b,0,b,q));
  }
  // Output
  return Rcpp::wrap(Out);
}