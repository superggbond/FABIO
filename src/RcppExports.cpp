// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SampleZ
void SampleZ(const arma::vec& y, const arma::vec& z_hat, arma::vec& z);
RcppExport SEXP _FABIO_SampleZ(SEXP ySEXP, SEXP z_hatSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_hat(z_hatSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    SampleZ(y, z_hat, z);
    return R_NilValue;
END_RCPP
}
// CalcPosterior
double CalcPosterior(List cHyp, const double ng_test, const double yty, const arma::mat& Xgamma, const int s_size, const arma::mat& XtX, const arma::vec& Xty, arma::vec& Xb, arma::vec& beta);
RcppExport SEXP _FABIO_CalcPosterior(SEXP cHypSEXP, SEXP ng_testSEXP, SEXP ytySEXP, SEXP XgammaSEXP, SEXP s_sizeSEXP, SEXP XtXSEXP, SEXP XtySEXP, SEXP XbSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cHyp(cHypSEXP);
    Rcpp::traits::input_parameter< const double >::type ng_test(ng_testSEXP);
    Rcpp::traits::input_parameter< const double >::type yty(ytySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xgamma(XgammaSEXP);
    Rcpp::traits::input_parameter< const int >::type s_size(s_sizeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcPosterior(cHyp, ng_test, yty, Xgamma, s_size, XtX, Xty, Xb, beta));
    return rcpp_result_gen;
END_RCPP
}
// CalcXtX
void CalcXtX(const arma::mat& X, const arma::vec& y, const int s_size, arma::mat& XtX, arma::vec& Xty);
RcppExport SEXP _FABIO_CalcXtX(SEXP XSEXP, SEXP ySEXP, SEXP s_sizeSEXP, SEXP XtXSEXP, SEXP XtySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int >::type s_size(s_sizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xty(XtySEXP);
    CalcXtX(X, y, s_size, XtX, Xty);
    return R_NilValue;
END_RCPP
}
// CalcCC_PVEnZ
void CalcCC_PVEnZ(const arma::vec& Xb, const double ni_test, List& cHyp, arma::vec& z_hat);
RcppExport SEXP _FABIO_CalcCC_PVEnZ(SEXP XbSEXP, SEXP ni_testSEXP, SEXP cHypSEXP, SEXP z_hatSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< const double >::type ni_test(ni_testSEXP);
    Rcpp::traits::input_parameter< List& >::type cHyp(cHypSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z_hat(z_hatSEXP);
    CalcCC_PVEnZ(Xb, ni_test, cHyp, z_hat);
    return R_NilValue;
END_RCPP
}
// ProposeH
List ProposeH(const List& cHyp_old, int rep, const double h_max, const double h_min, const double h_scale);
RcppExport SEXP _FABIO_ProposeH(SEXP cHyp_oldSEXP, SEXP repSEXP, SEXP h_maxSEXP, SEXP h_minSEXP, SEXP h_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type cHyp_old(cHyp_oldSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    Rcpp::traits::input_parameter< const double >::type h_max(h_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< const double >::type h_scale(h_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(ProposeH(cHyp_old, rep, h_max, h_min, h_scale));
    return rcpp_result_gen;
END_RCPP
}
// ProposePi
double ProposePi(const List& cHyp_old, List& cHyp_new, int rep, double logp_max, double logp_min, double logp_scale, double beta_a, double beta_b);
RcppExport SEXP _FABIO_ProposePi(SEXP cHyp_oldSEXP, SEXP cHyp_newSEXP, SEXP repSEXP, SEXP logp_maxSEXP, SEXP logp_minSEXP, SEXP logp_scaleSEXP, SEXP beta_aSEXP, SEXP beta_bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type cHyp_old(cHyp_oldSEXP);
    Rcpp::traits::input_parameter< List& >::type cHyp_new(cHyp_newSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    Rcpp::traits::input_parameter< double >::type logp_max(logp_maxSEXP);
    Rcpp::traits::input_parameter< double >::type logp_min(logp_minSEXP);
    Rcpp::traits::input_parameter< double >::type logp_scale(logp_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type beta_a(beta_aSEXP);
    Rcpp::traits::input_parameter< double >::type beta_b(beta_bSEXP);
    rcpp_result_gen = Rcpp::wrap(ProposePi(cHyp_old, cHyp_new, rep, logp_max, logp_min, logp_scale, beta_a, beta_b));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_iter
List mcmc_iter(const int total_step, const int w_step, const int r_pace, const int w_pace, const int n_mh, const int ng_test, const int ni_test, const double h_max, const double h_min, const double h_scale, const int g_max, const int g_min, const double logp_max, const double logp_min, const double logp_scale, const double beta_a, const double beta_b, const arma::vec& y, arma::vec& z_hat, arma::vec& z, arma::vec& rank_old, arma::vec& beta_old, arma::vec& beta_new, arma::vec& Xtz_old, arma::vec& Xtz_new, arma::vec& Xb_old, arma::vec& Xb_new, const arma::vec& p_gamma, const arma::ivec pos_vec, List& cHyp_old, const arma::mat& X, arma::mat& Xgamma_old, arma::mat& Xgamma_new, arma::mat& XtX_old, arma::mat& XtX_new, arma::mat& Result_hyp, arma::mat& Result_gamma, arma::mat& beta_g, bool display_progress);
RcppExport SEXP _FABIO_mcmc_iter(SEXP total_stepSEXP, SEXP w_stepSEXP, SEXP r_paceSEXP, SEXP w_paceSEXP, SEXP n_mhSEXP, SEXP ng_testSEXP, SEXP ni_testSEXP, SEXP h_maxSEXP, SEXP h_minSEXP, SEXP h_scaleSEXP, SEXP g_maxSEXP, SEXP g_minSEXP, SEXP logp_maxSEXP, SEXP logp_minSEXP, SEXP logp_scaleSEXP, SEXP beta_aSEXP, SEXP beta_bSEXP, SEXP ySEXP, SEXP z_hatSEXP, SEXP zSEXP, SEXP rank_oldSEXP, SEXP beta_oldSEXP, SEXP beta_newSEXP, SEXP Xtz_oldSEXP, SEXP Xtz_newSEXP, SEXP Xb_oldSEXP, SEXP Xb_newSEXP, SEXP p_gammaSEXP, SEXP pos_vecSEXP, SEXP cHyp_oldSEXP, SEXP XSEXP, SEXP Xgamma_oldSEXP, SEXP Xgamma_newSEXP, SEXP XtX_oldSEXP, SEXP XtX_newSEXP, SEXP Result_hypSEXP, SEXP Result_gammaSEXP, SEXP beta_gSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type total_step(total_stepSEXP);
    Rcpp::traits::input_parameter< const int >::type w_step(w_stepSEXP);
    Rcpp::traits::input_parameter< const int >::type r_pace(r_paceSEXP);
    Rcpp::traits::input_parameter< const int >::type w_pace(w_paceSEXP);
    Rcpp::traits::input_parameter< const int >::type n_mh(n_mhSEXP);
    Rcpp::traits::input_parameter< const int >::type ng_test(ng_testSEXP);
    Rcpp::traits::input_parameter< const int >::type ni_test(ni_testSEXP);
    Rcpp::traits::input_parameter< const double >::type h_max(h_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< const double >::type h_scale(h_scaleSEXP);
    Rcpp::traits::input_parameter< const int >::type g_max(g_maxSEXP);
    Rcpp::traits::input_parameter< const int >::type g_min(g_minSEXP);
    Rcpp::traits::input_parameter< const double >::type logp_max(logp_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type logp_min(logp_minSEXP);
    Rcpp::traits::input_parameter< const double >::type logp_scale(logp_scaleSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_a(beta_aSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_b(beta_bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z_hat(z_hatSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rank_old(rank_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_new(beta_newSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xtz_old(Xtz_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xtz_new(Xtz_newSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xb_old(Xb_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xb_new(Xb_newSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p_gamma(p_gammaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec >::type pos_vec(pos_vecSEXP);
    Rcpp::traits::input_parameter< List& >::type cHyp_old(cHyp_oldSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xgamma_old(Xgamma_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xgamma_new(Xgamma_newSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XtX_old(XtX_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XtX_new(XtX_newSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Result_hyp(Result_hypSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Result_gamma(Result_gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_g(beta_gSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_iter(total_step, w_step, r_pace, w_pace, n_mh, ng_test, ni_test, h_max, h_min, h_scale, g_max, g_min, logp_max, logp_min, logp_scale, beta_a, beta_b, y, z_hat, z, rank_old, beta_old, beta_new, Xtz_old, Xtz_new, Xb_old, Xb_new, p_gamma, pos_vec, cHyp_old, X, Xgamma_old, Xgamma_new, XtX_old, XtX_new, Result_hyp, Result_gamma, beta_g, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FABIO_SampleZ", (DL_FUNC) &_FABIO_SampleZ, 3},
    {"_FABIO_CalcPosterior", (DL_FUNC) &_FABIO_CalcPosterior, 9},
    {"_FABIO_CalcXtX", (DL_FUNC) &_FABIO_CalcXtX, 5},
    {"_FABIO_CalcCC_PVEnZ", (DL_FUNC) &_FABIO_CalcCC_PVEnZ, 4},
    {"_FABIO_ProposeH", (DL_FUNC) &_FABIO_ProposeH, 5},
    {"_FABIO_ProposePi", (DL_FUNC) &_FABIO_ProposePi, 8},
    {"_FABIO_mcmc_iter", (DL_FUNC) &_FABIO_mcmc_iter, 39},
    {NULL, NULL, 0}
};

RcppExport void R_init_FABIO(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
