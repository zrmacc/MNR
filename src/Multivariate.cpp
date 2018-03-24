// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for fitting multivariate outcome model
// Updated: 180315
#include <RcppEigen.h>

// Target Working Vector
// 
// Construct target working vector, i.e. the outcome using for updating \eqn{\beta}.
// 
// @param yt Target outcome
// @param eS Surrogate residual vector
// @param Ltt Target variance
// @param Lst Surrogate target covariance vector
Eigen::VectorXd mnrZT(const Eigen::VectorXd yt, const Eigen::VectorXd eS,
           const double Ltt, const Eigen::VectorXd Lst){
  // Dimension
  const int n = yt.size();
  const int k = Lst.size();
  // Output structure
  Eigen::VectorXd Out(n);
  Out.setZero();
  // Weight
  Eigen::RowVectorXd w = Lst.transpose();
  const Eigen::RowVectorXd W = w/Ltt;
  // Calculate entries
  for(int i=0; i<n; i++){
    Out(i) = (yt(i)+W*eS.segment(i*k,k));
  }
  return Out;
}

// Surrogate Working Vector
// 
// Constrcut srrogate working vector, i.e. the outcome using for updating \eqn{\alpha}.
// 
// @param eT Target residual vector
// @param S Surrogate outcome vector
// @param Lst Surrogate target covariance vector
// @param Lss Surrogate covariance matrix
Eigen::VectorXd mnrZS(const Eigen::VectorXd eT, const Eigen::VectorXd S,
                      const Eigen::VectorXd Lst, const Eigen::MatrixXd Lss){
  // Dimension
  const int n = eT.size();
  const int m = Lst.size();
  // Output structure
  Eigen::VectorXd Out(n*m);
  Out.setZero();
  // Loop over observations
  for(int i=0; i<n; i++){
    Out.segment(i*m,m) = (Lst*eT(i))+(Lss*S.segment(i*m,m));
  }
  return Out;
}

// Update Alpha
// 
// Calculate \eqn{\alpha_{1}} for the multivariate outcome model. 
// 
// @param n Sample size
// @param m Surrogate outcomes
// @param zS Surrogate working vector
// @param Zs Surrogate design
// @param Lss Surrogate covariance matrix
Eigen::VectorXd updateA(const int n, const int m, const Eigen::VectorXd zS,
                        const Eigen::MatrixXd Zs, const Eigen::MatrixXd Lss){
  // Dimension
  const int q = Zs.cols();
  // Weight matrix
  Eigen::MatrixXd W(q,q);
  W.setZero();
  for(int i=0; i<n; i++){
    W += (Zs.block(i*m,0,m,q).transpose())*Lss*Zs.block(i*m,0,m,q);
  }
  // Update alpha
  const Eigen::VectorXd a1 = W.llt().solve(Zs.transpose()*zS);
  return a1;
}

// Residual Matrix
// 
// Calculate \eqn{E} for the multivariate outcome model
// 
// @param m Surrogate outcomes
// @param eT Target residual
// @param S Surrogate outcome
// @param Zs Surrogate design
// @param a1 Alpha
Eigen::MatrixXd resid(const int m, const Eigen::VectorXd eT, 
                      const Eigen::VectorXd S, const Eigen::MatrixXd Zs,
                      const Eigen::VectorXd a1){
  // Observations
  const int n = eT.size();
  const int q = Zs.cols();
  // Output structure
  Eigen::MatrixXd E(n,m+1);
  E.setZero();
  // Target residual
  E.col(0) = eT;
  // Surrogate residuals
  for(int i=0; i<n; i++){
    E.block(i,1,1,m) = (S.segment(i*m,m)-(Zs.block(i*m,0,m,q)*a1)).transpose();
  }
  return E;
}

//' Multivariate Parameter Update
//' 
//' Parameter update for multivariate outcome model.
//' 
//' @param yt Target outcome
//' @param S Surrogate outcome vector
//' @param Zt Target design
//' @param Zs Surrogate design
//' @param At Target inc. proj.
//' @param b0 Current beta
//' @param a0 Current alpha
//' @param s0 Current sigma
//' @export 
// [[Rcpp::export]]

SEXP updateMNR(const Eigen::Map<Eigen::VectorXd> yt, const Eigen::Map<Eigen::VectorXd> S,
               const Eigen::Map<Eigen::MatrixXd> Zt, const Eigen::Map<Eigen::MatrixXd> Zs,
               const Eigen::Map<Eigen::MatrixXd> At,
               const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a0,
               const Eigen::Map<Eigen::MatrixXd> s0){
  // Observations
  const int n = Zt.rows();
  const int m = s0.cols()-1;
  // Invert Sigma
  const Eigen::MatrixXd l0 = s0.completeOrthogonalDecomposition().pseudoInverse();
  // Partition Lambda
  const double Ltt = l0(0,0);
  const Eigen::VectorXd Lst = l0.block(1,0,m,1);
  const Eigen::MatrixXd Lss = l0.block(1,1,m,m);
  // Surrogate residual
  const Eigen::VectorXd eS = (S-Zs*a0);
  // Calculate target working vector
  const Eigen::VectorXd zT = mnrZT(yt,eS,Ltt,Lst);
  // Update beta
  const Eigen::VectorXd b1 = At*zT;
  // Target residual
  const Eigen::VectorXd eT = (yt-Zt*b1);
  // Calculate surrogate working vector
  const Eigen::VectorXd zS = mnrZS(eT,S,Lst,Lss);
  // Update alpha
  const Eigen::VectorXd a1 = updateA(n,m,zS,Zs,Lss);
  // Residual matrix
  const Eigen::MatrixXd E = resid(m,eT,S,Zs,a1);
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


