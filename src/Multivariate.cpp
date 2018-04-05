// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for fitting multivariate outcome model
// Updated: 180404
#include <RcppEigen.h>

//' Initialize Alpha
//' 
//' @param Xi Surrogate design matrix
//' @param s Surrogate outcome
//' @export 
// [[Rcpp::export]]

SEXP alpha0(const Eigen::MatrixXd Xi, const Eigen::VectorXd s){
  const Eigen::VectorXd a0 = (Xi.transpose()*Xi).llt().solve(Xi.transpose()*s);
  return Rcpp::wrap(a0);
}

// Target Working Vector
// 
// Construct target working vector, i.e. the outcome using for updating \eqn{\beta}.
// 
// @param t Target outcome
// @param es Surrogate residual vector
// @param Ltt Target precision
// @param Lst Surrogate target precision
Eigen::VectorXd mnrZT(const Eigen::VectorXd t, const Eigen::VectorXd es,
                      const double Ltt, const Eigen::VectorXd Lst){
  // Dimension
  const int n = t.size();
  const int k = Lst.size();
  // Output structure
  Eigen::VectorXd Out(n);
  Out.setZero();
  // Weight
  Eigen::RowVectorXd w = Lst.transpose();
  const Eigen::RowVectorXd W = w/Ltt;
  // Calculate entries
  for(int i=0; i<n; i++){
    Out(i) = (t(i)+W*es.segment(i*k,k));
  }
  return Out;
}

// Surrogate Working Vector
// 
// Constrcut srrogate working vector, i.e. the outcome using for updating \eqn{\alpha}.
// 
// @param t Target residual vector
// @param s Surrogate outcome vector
// @param Lst Surrogate target precision
// @param Lss Surrogate precision
Eigen::VectorXd mnrZS(const Eigen::VectorXd et, const Eigen::VectorXd s,
                      const Eigen::VectorXd Lst, const Eigen::MatrixXd Lss){
  // Dimension
  const int n = et.size();
  const int m = Lst.size();
  // Output structure
  Eigen::VectorXd Out(n*m);
  Out.setZero();
  // Loop over observations
  for(int i=0; i<n; i++){
    Out.segment(i*m,m) = (Lst*et(i))+(Lss*s.segment(i*m,m));
  }
  return Out;
}

// Update Alpha
// 
// Calculate \eqn{\alpha_{1}} for the multivariate outcome model. 
// 
// @param n Sample size
// @param m Surrogate outcomes
// @param zs Surrogate working vector
// @param Xi Surrogate design
// @param Lss Surrogate precision
Eigen::VectorXd updateA(const int n, const int m, const Eigen::VectorXd zs,
                        const Eigen::MatrixXd Xi, const Eigen::MatrixXd Lss){
  // Dimension
  const int q = Xi.cols();
  // Weight matrix
  Eigen::MatrixXd W(q,q);
  W.setZero();
  for(int i=0; i<n; i++){
    W += (Xi.block(i*m,0,m,q).transpose())*Lss*Xi.block(i*m,0,m,q);
  }
  // Update alpha
  const Eigen::VectorXd a1 = W.llt().solve(Xi.transpose()*zs);
  return a1;
}

// Residual Matrix
// 
// Calculate \eqn{E} for the multivariate outcome model
// 
// @param m Surrogate outcomes
// @param et Target residual
// @param s Surrogate outcome
// @param Xi Surrogate design
// @param a1 Alpha
Eigen::MatrixXd resid(const int m, const Eigen::VectorXd et, 
                      const Eigen::VectorXd s, const Eigen::MatrixXd Xi,
                      const Eigen::VectorXd a1){
  // Observations
  const int n = et.size();
  const int q = Xi.cols();
  // Output structure
  Eigen::MatrixXd E(n,m+1);
  E.setZero();
  // Target residual
  E.col(0) = et;
  // Surrogate residuals
  for(int i=0; i<n; i++){
    E.block(i,1,1,m) = (s.segment(i*m,m)-(Xi.block(i*m,0,m,q)*a1)).transpose();
  }
  return E;
}

//' Multivariate Parameter Update
//' 
//' Parameter update for multivariate outcome model.
//' 
//' @param t Target outcome
//' @param s Surrogate outcome vector
//' @param X Target design
//' @param Xi Surrogate design
//' @param B Target inc. proj.
//' @param b0 Current beta
//' @param a0 Current alpha
//' @param s0 Current sigma
//' @export 
// [[Rcpp::export]]

SEXP updateMNR(const Eigen::Map<Eigen::VectorXd> t, const Eigen::Map<Eigen::VectorXd> s,
               const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Xi,
               const Eigen::Map<Eigen::MatrixXd> B,
               const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a0,
               const Eigen::Map<Eigen::MatrixXd> s0){
  // Observations
  const int n = X.rows();
  // Surrogates
  const int m = s0.cols()-1;
  // Invert Sigma
  const Eigen::MatrixXd l0 = s0.completeOrthogonalDecomposition().pseudoInverse();
  // Partition Lambda
  const double Ltt = l0(0,0);
  const Eigen::VectorXd Lst = l0.block(1,0,m,1);
  const Eigen::MatrixXd Lss = l0.block(1,1,m,m);
  // Surrogate residual
  const Eigen::VectorXd es = (s-Xi*a0);
  // Calculate target working vector
  const Eigen::VectorXd zt = mnrZT(t,es,Ltt,Lst);
  // Update beta
  const Eigen::VectorXd b1 = B*zt;
  // Target residual
  const Eigen::VectorXd et = (t-X*b1);
  // Calculate surrogate working vector
  const Eigen::VectorXd zs = mnrZS(et,s,Lst,Lss);
  // Update alpha
  const Eigen::VectorXd a1 = updateA(n,m,zs,Xi,Lss);
  // Residual matrix
  const Eigen::MatrixXd E = resid(m,et,s,Xi,a1);
  // Covariance matrix
  const Eigen::MatrixXd ip = E.transpose()*E;
  const Eigen::MatrixXd s1 = (ip/n);
  // Output
  return Rcpp::List::create(Rcpp::Named("b1")=b1,Rcpp::Named("a1")=a1,Rcpp::Named("s1")=s1,Rcpp::Named("E1")=E);
}

//' Information
//'
//' Calculate \eqn{A'(I_{n\times n}\otimes L)B}, the general form of information matrices
//' for the multivariate outcome model. 
//'
//' @param n Observations
//' @param A Left hand matrix
//' @param B Right hand matrix
//' @param L Precision
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


