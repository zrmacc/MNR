// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Matrix trace
//'
//' Calculates the trace of \eqn{A}.
//'
//' @param A Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP tr(const Eigen::Map<Eigen::MatrixXd> A){
  const double t = A.diagonal().sum();
  return Rcpp::wrap(t);
}

//' Matrix matrix product
//'
//' Calculates \eqn{AB};
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP fastMMp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  const Eigen::MatrixXd C = A*B;
  return Rcpp::wrap(C);
}

//' Matrix Transpose
//'
//' Constructs \eqn{A'} from \eqn{A}.
//'
//' @param A Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP fastT(const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd At = A.transpose();
  return Rcpp::wrap(At);
}

//' Matrix Inner Product
//'
//' Calculates \eqn{A'B}.
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP fastIP(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  const Eigen::MatrixXd AtB = (A.transpose() * B);
  return Rcpp::wrap(AtB);
}

//' Matrix Inverse
//'
//' @param A Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP fastInv(const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd Ai = A.completeOrthogonalDecomposition().pseudoInverse();
  return Rcpp::wrap(Ai);
}

//' Matrix Determinant
//'
//' Calculates \eqn{\det(A)}.
//'
//' @param A Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP fastDet(const Eigen::Map<Eigen::MatrixXd> A){
  const double d = A.determinant();
  return Rcpp::wrap(d);
}

//' Quadratic Form
//' 
//' Calculates \eqn{x'Ax}.
//' 
//' @param x Numeric vector.
//' @param A Numeric matrix.
//' @export
// [[Rcpp::export]]
SEXP fastQF(const Eigen::Map<Eigen::VectorXd> x, const Eigen::Map<Eigen::MatrixXd> A){
  const double q = x.transpose()*A*x;
  return Rcpp::wrap(q);
}

//' Schur complement
//'
//' Calculates the efficient information \eqn{I_{11}-I_{12}I_{22}^{-1}I_{21}};
//'
//' @param I11 Information of target parameter
//' @param I22 Information of nuisance parameter
//' @param I12 Cross information between target and nuisance parameters
//' @export
//'
// [[Rcpp::export]]
SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> I11, const Eigen::Map<Eigen::MatrixXd> I22,
            const Eigen::Map<Eigen::MatrixXd> I12){
  // Kernel matrix
  const Eigen::MatrixXd K = I11 - I12 * I22.llt().solve(I12.transpose());
  return Rcpp::wrap(K);
}