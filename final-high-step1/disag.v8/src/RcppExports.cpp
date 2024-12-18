// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Arma_colSums
arma::rowvec Arma_colSums(const arma::mat& x);
RcppExport SEXP _disag_v8_Arma_colSums(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Arma_colSums(x));
    return rcpp_result_gen;
END_RCPP
}
// fitCpp
double fitCpp(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const arma::mat& S5, const arma::mat& S6, const arma::mat& S7, const arma::mat& S8, const arma::mat& S9, const double& sldiagA, const int& it);
RcppExport SEXP _disag_v8_fitCpp(SEXP rhoSEXP, SEXP QtySEXP, SEXP RSEXP, SEXP RtRSEXP, SEXP rSEXP, SEXP NSEXP, SEXP b0SEXP, SEXP S1SEXP, SEXP S2SEXP, SEXP S3SEXP, SEXP S4SEXP, SEXP S5SEXP, SEXP S6SEXP, SEXP S7SEXP, SEXP S8SEXP, SEXP S9SEXP, SEXP sldiagASEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Qty(QtySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type RtR(RtRSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S3(S3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S4(S4SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S5(S5SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S6(S6SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S7(S7SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S8(S8SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S9(S9SEXP);
    Rcpp::traits::input_parameter< const double& >::type sldiagA(sldiagASEXP);
    Rcpp::traits::input_parameter< const int& >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(fitCpp(rho, Qty, R, RtR, r, N, b0, S1, S2, S3, S4, S5, S6, S7, S8, S9, sldiagA, it));
    return rcpp_result_gen;
END_RCPP
}
// hbCpp
arma::vec hbCpp(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const arma::mat& S5, const arma::mat& S6, const arma::mat& S7, const arma::mat& S8, const arma::mat& S9);
RcppExport SEXP _disag_v8_hbCpp(SEXP rhoSEXP, SEXP QtySEXP, SEXP RSEXP, SEXP RtRSEXP, SEXP rSEXP, SEXP NSEXP, SEXP b0SEXP, SEXP S1SEXP, SEXP S2SEXP, SEXP S3SEXP, SEXP S4SEXP, SEXP S5SEXP, SEXP S6SEXP, SEXP S7SEXP, SEXP S8SEXP, SEXP S9SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Qty(QtySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type RtR(RtRSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S3(S3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S4(S4SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S5(S5SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S6(S6SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S7(S7SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S8(S8SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S9(S9SEXP);
    rcpp_result_gen = Rcpp::wrap(hbCpp(rho, Qty, R, RtR, r, N, b0, S1, S2, S3, S4, S5, S6, S7, S8, S9));
    return rcpp_result_gen;
END_RCPP
}
// ll_at_min
double ll_at_min(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const arma::mat& S5, const arma::mat& S6, const arma::mat& S7, const arma::mat& S8, const arma::mat& S9, const double& sldiagA);
RcppExport SEXP _disag_v8_ll_at_min(SEXP rhoSEXP, SEXP QtySEXP, SEXP RSEXP, SEXP RtRSEXP, SEXP rSEXP, SEXP NSEXP, SEXP b0SEXP, SEXP S1SEXP, SEXP S2SEXP, SEXP S3SEXP, SEXP S4SEXP, SEXP S5SEXP, SEXP S6SEXP, SEXP S7SEXP, SEXP S8SEXP, SEXP S9SEXP, SEXP sldiagASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Qty(QtySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type RtR(RtRSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S3(S3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S4(S4SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S5(S5SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S6(S6SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S7(S7SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S8(S8SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S9(S9SEXP);
    Rcpp::traits::input_parameter< const double& >::type sldiagA(sldiagASEXP);
    rcpp_result_gen = Rcpp::wrap(ll_at_min(rho, Qty, R, RtR, r, N, b0, S1, S2, S3, S4, S5, S6, S7, S8, S9, sldiagA));
    return rcpp_result_gen;
END_RCPP
}
// getIndex
arma::vec getIndex(const arma::vec& v, const int b);
RcppExport SEXP _disag_v8_getIndex(SEXP vSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(getIndex(v, b));
    return rcpp_result_gen;
END_RCPP
}
// getRowIndex_NA
arma::uvec getRowIndex_NA(const arma::mat& m);
RcppExport SEXP _disag_v8_getRowIndex_NA(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(getRowIndex_NA(m));
    return rcpp_result_gen;
END_RCPP
}
// tpmm_agg_noblock_fast
Rcpp::List tpmm_agg_noblock_fast(const arma::mat& Xs, const arma::vec& groups, const arma::vec& blocks, const arma::vec& b);
RcppExport SEXP _disag_v8_tpmm_agg_noblock_fast(SEXP XsSEXP, SEXP groupsSEXP, SEXP blocksSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(tpmm_agg_noblock_fast(Xs, groups, blocks, b));
    return rcpp_result_gen;
END_RCPP
}
// tpmm_vec_q_noblock
arma::vec tpmm_vec_q_noblock(const arma::mat& Xs, const arma::mat& b, int type);
RcppExport SEXP _disag_v8_tpmm_vec_q_noblock(SEXP XsSEXP, SEXP bSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(tpmm_vec_q_noblock(Xs, b, type));
    return rcpp_result_gen;
END_RCPP
}
// tpmm_vec_noblock
arma::vec tpmm_vec_noblock(const arma::mat& Xs, const arma::mat& b);
RcppExport SEXP _disag_v8_tpmm_vec_noblock(SEXP XsSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(tpmm_vec_noblock(Xs, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_disag_v8_Arma_colSums", (DL_FUNC) &_disag_v8_Arma_colSums, 1},
    {"_disag_v8_fitCpp", (DL_FUNC) &_disag_v8_fitCpp, 18},
    {"_disag_v8_hbCpp", (DL_FUNC) &_disag_v8_hbCpp, 16},
    {"_disag_v8_ll_at_min", (DL_FUNC) &_disag_v8_ll_at_min, 17},
    {"_disag_v8_getIndex", (DL_FUNC) &_disag_v8_getIndex, 2},
    {"_disag_v8_getRowIndex_NA", (DL_FUNC) &_disag_v8_getRowIndex_NA, 1},
    {"_disag_v8_tpmm_agg_noblock_fast", (DL_FUNC) &_disag_v8_tpmm_agg_noblock_fast, 4},
    {"_disag_v8_tpmm_vec_q_noblock", (DL_FUNC) &_disag_v8_tpmm_vec_q_noblock, 3},
    {"_disag_v8_tpmm_vec_noblock", (DL_FUNC) &_disag_v8_tpmm_vec_noblock, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_disag_v8(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
