// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// l0araC
List l0araC(arma::mat x, arma::vec y, String family, double lam, int maxit, double eps);
RcppExport SEXP _l0ara_l0araC(SEXP xSEXP, SEXP ySEXP, SEXP familySEXP, SEXP lamSEXP, SEXP maxitSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(l0araC(x, y, family, lam, maxit, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_l0ara_l0araC", (DL_FUNC) &_l0ara_l0araC, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_l0ara(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}