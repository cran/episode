// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ode_field
Rcpp::List ode_field(Rcpp::List ode_struct, arma::vec x, Rcpp::Nullable<Rcpp::List> param_, bool differentials);
RcppExport SEXP _episode_ode_field(SEXP ode_structSEXP, SEXP xSEXP, SEXP param_SEXP, SEXP differentialsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ode_struct(ode_structSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type param_(param_SEXP);
    Rcpp::traits::input_parameter< bool >::type differentials(differentialsSEXP);
    rcpp_result_gen = Rcpp::wrap(ode_field(ode_struct, x, param_, differentials));
    return rcpp_result_gen;
END_RCPP
}
// aim_design
Rcpp::List aim_design(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_, arma::mat x, Rcpp::Nullable<Rcpp::List> param_);
RcppExport SEXP _episode_aim_design(SEXP ode_structSEXP, SEXP opt_structSEXP, SEXP sc_SEXP, SEXP xSEXP, SEXP param_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ode_struct(ode_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type opt_struct(opt_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sc_(sc_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type param_(param_SEXP);
    rcpp_result_gen = Rcpp::wrap(aim_design(ode_struct, opt_struct, sc_, x, param_));
    return rcpp_result_gen;
END_RCPP
}
// numint
arma::mat numint(arma::vec time, arma::mat x, std::string type, arma::mat A, arma::mat B);
RcppExport SEXP _episode_numint(SEXP timeSEXP, SEXP xSEXP, SEXP typeSEXP, SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(numint(time, x, type, A, B));
    return rcpp_result_gen;
END_RCPP
}
// bull
Rcpp::List bull(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_, Rcpp::List params_, arma::cube x0s, bool trace);
RcppExport SEXP _episode_bull(SEXP ode_structSEXP, SEXP opt_structSEXP, SEXP sc_SEXP, SEXP params_SEXP, SEXP x0sSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ode_struct(ode_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type opt_struct(opt_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sc_(sc_SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params_(params_SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type x0s(x0sSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(bull(ode_struct, opt_struct, sc_, params_, x0s, trace));
    return rcpp_result_gen;
END_RCPP
}
// cattle
Rcpp::List cattle(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_, Rcpp::Nullable<Rcpp::List> param_, arma::mat x, bool trace);
RcppExport SEXP _episode_cattle(SEXP ode_structSEXP, SEXP opt_structSEXP, SEXP sc_SEXP, SEXP param_SEXP, SEXP xSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ode_struct(ode_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type opt_struct(opt_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sc_(sc_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type param_(param_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(cattle(ode_struct, opt_struct, sc_, param_, x, trace));
    return rcpp_result_gen;
END_RCPP
}
// bronc
Rcpp::List bronc(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_, Rcpp::List params_, arma::cube x0s, arma::uvec indices, bool adjust_lambda, bool adjust_scales, bool adjust_weights, bool trace);
RcppExport SEXP _episode_bronc(SEXP ode_structSEXP, SEXP opt_structSEXP, SEXP sc_SEXP, SEXP params_SEXP, SEXP x0sSEXP, SEXP indicesSEXP, SEXP adjust_lambdaSEXP, SEXP adjust_scalesSEXP, SEXP adjust_weightsSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ode_struct(ode_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type opt_struct(opt_structSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sc_(sc_SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params_(params_SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type x0s(x0sSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust_lambda(adjust_lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust_scales(adjust_scalesSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust_weights(adjust_weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(bronc(ode_struct, opt_struct, sc_, params_, x0s, indices, adjust_lambda, adjust_scales, adjust_weights, trace));
    return rcpp_result_gen;
END_RCPP
}
// ode_solve
Rcpp::List ode_solve(Rcpp::List ode_struct, arma::mat x0, Rcpp::List param_, arma::vec time, bool sensitivity, bool approx);
RcppExport SEXP _episode_ode_solve(SEXP ode_structSEXP, SEXP x0SEXP, SEXP param_SEXP, SEXP timeSEXP, SEXP sensitivitySEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ode_struct(ode_structSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param_(param_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< bool >::type sensitivity(sensitivitySEXP);
    Rcpp::traits::input_parameter< bool >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(ode_solve(ode_struct, x0, param_, time, sensitivity, approx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_episode_ode_field", (DL_FUNC) &_episode_ode_field, 4},
    {"_episode_aim_design", (DL_FUNC) &_episode_aim_design, 5},
    {"_episode_numint", (DL_FUNC) &_episode_numint, 5},
    {"_episode_bull", (DL_FUNC) &_episode_bull, 6},
    {"_episode_cattle", (DL_FUNC) &_episode_cattle, 6},
    {"_episode_bronc", (DL_FUNC) &_episode_bronc, 10},
    {"_episode_ode_solve", (DL_FUNC) &_episode_ode_solve, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_episode(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}