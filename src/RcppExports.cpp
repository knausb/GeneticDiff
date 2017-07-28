// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rgt
Rcpp::CharacterMatrix rgt(int nsamp, int nvar, int pphased, Rcpp::IntegerVector pploid, Rcpp::IntegerVector pallele);
RcppExport SEXP GeneticDiff_rgt(SEXP nsampSEXP, SEXP nvarSEXP, SEXP pphasedSEXP, SEXP pploidSEXP, SEXP palleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type pphased(pphasedSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pploid(pploidSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pallele(palleleSEXP);
    rcpp_result_gen = Rcpp::wrap(rgt(nsamp, nvar, pphased, pploid, pallele));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"GeneticDiff_rgt", (DL_FUNC) &GeneticDiff_rgt, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GeneticDiff(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
