// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// updateBVR
SEXP updateBVR(const Eigen::Map<Eigen::VectorXd> yt, const Eigen::Map<Eigen::VectorXd> ys, const Eigen::Map<Eigen::MatrixXd> Zt, const Eigen::Map<Eigen::MatrixXd> Zs, const Eigen::Map<Eigen::MatrixXd> At, const Eigen::Map<Eigen::MatrixXd> As, const Eigen::Map<Eigen::VectorXd> b0, const Eigen::Map<Eigen::VectorXd> a0, const Eigen::Map<Eigen::MatrixXd> s0);
RcppExport SEXP _MNR_updateBVR(SEXP ytSEXP, SEXP ysSEXP, SEXP ZtSEXP, SEXP ZsSEXP, SEXP AtSEXP, SEXP AsSEXP, SEXP b0SEXP, SEXP a0SEXP, SEXP s0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type yt(ytSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Zs(ZsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type At(AtSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type As(AsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type s0(s0SEXP);
    rcpp_result_gen = Rcpp::wrap(updateBVR(yt, ys, Zt, Zs, At, As, b0, a0, s0));
    return rcpp_result_gen;
END_RCPP
}
// fastT
SEXP fastT(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _MNR_fastT(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastT(A));
    return rcpp_result_gen;
END_RCPP
}
// fastIP
SEXP fastIP(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _MNR_fastIP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fastIP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// fastInv
SEXP fastInv(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _MNR_fastInv(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastInv(A));
    return rcpp_result_gen;
END_RCPP
}
// fastDet
SEXP fastDet(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _MNR_fastDet(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastDet(A));
    return rcpp_result_gen;
END_RCPP
}
// fastQF
SEXP fastQF(const Eigen::Map<Eigen::VectorXd> x, const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _MNR_fastQF(SEXP xSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastQF(x, A));
    return rcpp_result_gen;
END_RCPP
}
// incP
SEXP incP(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _MNR_incP(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(incP(A));
    return rcpp_result_gen;
END_RCPP
}
// SchurC
SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> I11, const Eigen::Map<Eigen::MatrixXd> I22, const Eigen::Map<Eigen::MatrixXd> I12);
RcppExport SEXP _MNR_SchurC(SEXP I11SEXP, SEXP I22SEXP, SEXP I12SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type I11(I11SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type I22(I22SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type I12(I12SEXP);
    rcpp_result_gen = Rcpp::wrap(SchurC(I11, I22, I12));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MNR_updateBVR", (DL_FUNC) &_MNR_updateBVR, 9},
    {"_MNR_fastT", (DL_FUNC) &_MNR_fastT, 1},
    {"_MNR_fastIP", (DL_FUNC) &_MNR_fastIP, 2},
    {"_MNR_fastInv", (DL_FUNC) &_MNR_fastInv, 1},
    {"_MNR_fastDet", (DL_FUNC) &_MNR_fastDet, 1},
    {"_MNR_fastQF", (DL_FUNC) &_MNR_fastQF, 2},
    {"_MNR_incP", (DL_FUNC) &_MNR_incP, 1},
    {"_MNR_SchurC", (DL_FUNC) &_MNR_SchurC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MNR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
