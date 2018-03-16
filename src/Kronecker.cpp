// [[Rcpp::depends(RcppEigen)]]
// Purpose: Sparse kronecker product
// Updated: 180315
#include <RcppEigen.h>

// Note: R's implementation is faster. 

//' Sparse Kronecker Product
//' 
//' Forms the product \eqn{I_{n\times n}\otimes A}.
//' 
//' @param n Dimension of identity matrix
//' @param A Numeric matrix
// [[Rcpp::export]]

SEXP sKroneckerP(const int n, const Eigen::Map<Eigen::MatrixXd> A){
  // Generate identity matrix
  const Eigen::VectorXd i = Eigen::VectorXd::Constant(n,1);
  const Eigen::MatrixXd I = i.asDiagonal();
  // Output structure
  const int a = A.rows();
  const int b = A.cols();
  Eigen::MatrixXd Out(n*a,n*b);
  Out.setZero();
  // Fill blocks
  for(int i=0; i<n; i++){
    Out.block(i*a,i*b,a,b) = A;
  }
  // Output
  return Rcpp::wrap(Out);
}