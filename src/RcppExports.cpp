// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BIC_logit_NDH
double BIC_logit_NDH(arma::sp_mat A, Rcpp::List object);
RcppExport SEXP _JANE_BIC_logit_NDH(SEXP ASEXP, SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(BIC_logit_NDH(A, object));
    return rcpp_result_gen;
END_RCPP
}
// BIC_logit_RS
double BIC_logit_RS(arma::sp_mat A, Rcpp::List object);
RcppExport SEXP _JANE_BIC_logit_RS(SEXP ASEXP, SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(BIC_logit_RS(A, object));
    return rcpp_result_gen;
END_RCPP
}
// BIC_logit_RSR
double BIC_logit_RSR(arma::sp_mat A, Rcpp::List object);
RcppExport SEXP _JANE_BIC_logit_RSR(SEXP ASEXP, SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(BIC_logit_RSR(A, object));
    return rcpp_result_gen;
END_RCPP
}
// BIC_ICL_MBC
Rcpp::List BIC_ICL_MBC(Rcpp::List object);
RcppExport SEXP _JANE_BIC_ICL_MBC(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(BIC_ICL_MBC(object));
    return rcpp_result_gen;
END_RCPP
}
// log_like_C
double log_like_C(arma::colvec par, arma::mat X, arma::colvec y);
RcppExport SEXP _JANE_log_like_C(SEXP parSEXP, SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(log_like_C(par, X, y));
    return rcpp_result_gen;
END_RCPP
}
// gradient_C
arma::colvec gradient_C(arma::colvec par, arma::mat X, arma::colvec y);
RcppExport SEXP _JANE_gradient_C(SEXP parSEXP, SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(gradient_C(par, X, y));
    return rcpp_result_gen;
END_RCPP
}
// compute_dist
void compute_dist(arma::mat U, arma::mat& distances, std::string model, arma::mat X, arma::mat indices, bool downsampling);
RcppExport SEXP _JANE_compute_dist(SEXP USEXP, SEXP distancesSEXP, SEXP modelSEXP, SEXP XSEXP, SEXP indicesSEXP, SEXP downsamplingSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< bool >::type downsampling(downsamplingSEXP);
    compute_dist(U, distances, model, X, indices, downsampling);
    return R_NilValue;
END_RCPP
}
// log_Q
double log_Q(arma::sp_mat A, arma::mat U, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, arma::colvec p, arma::rowvec a, double b, double c, arma::mat G, arma::colvec nu, double e, double f, void * X, void * n_control, void * model);
RcppExport SEXP _JANE_log_Q(SEXP ASEXP, SEXP USEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP prob_matrixSEXP, SEXP betaSEXP, SEXP pSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP GSEXP, SEXP nuSEXP, SEXP eSEXP, SEXP fSEXP, SEXP XSEXP, SEXP n_controlSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    Rcpp::traits::input_parameter< void * >::type X(XSEXP);
    Rcpp::traits::input_parameter< void * >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< void * >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Q(A, U, mus, omegas, prob_matrix, beta, p, a, b, c, G, nu, e, f, X, n_control, model));
    return rcpp_result_gen;
END_RCPP
}
// log_Q_RE
double log_Q_RE(arma::sp_mat A, arma::mat U, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, arma::colvec p, arma::rowvec a, double b, double c, arma::mat G, arma::colvec nu, arma::colvec e, arma::mat f, arma::mat X, Rcpp::String model, void * n_control);
RcppExport SEXP _JANE_log_Q_RE(SEXP ASEXP, SEXP USEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP prob_matrixSEXP, SEXP betaSEXP, SEXP pSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP GSEXP, SEXP nuSEXP, SEXP eSEXP, SEXP fSEXP, SEXP XSEXP, SEXP modelSEXP, SEXP n_controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type model(modelSEXP);
    Rcpp::traits::input_parameter< void * >::type n_control(n_controlSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Q_RE(A, U, mus, omegas, prob_matrix, beta, p, a, b, c, G, nu, e, f, X, model, n_control));
    return rcpp_result_gen;
END_RCPP
}
// draw_A_NDH_c
arma::sp_mat draw_A_NDH_c(arma::mat U, double beta0);
RcppExport SEXP _JANE_draw_A_NDH_c(SEXP USEXP, SEXP beta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type beta0(beta0SEXP);
    rcpp_result_gen = Rcpp::wrap(draw_A_NDH_c(U, beta0));
    return rcpp_result_gen;
END_RCPP
}
// draw_A_RS_c
arma::sp_mat draw_A_RS_c(arma::mat U, double beta0, arma::colvec s);
RcppExport SEXP _JANE_draw_A_RS_c(SEXP USEXP, SEXP beta0SEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_A_RS_c(U, beta0, s));
    return rcpp_result_gen;
END_RCPP
}
// draw_A_RSR_c
arma::sp_mat draw_A_RSR_c(arma::mat U, double beta0, arma::colvec s, arma::colvec r);
RcppExport SEXP _JANE_draw_A_RSR_c(SEXP USEXP, SEXP beta0SEXP, SEXP sSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_A_RSR_c(U, beta0, s, r));
    return rcpp_result_gen;
END_RCPP
}
// update_U
void update_U(arma::mat& U, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, void * X, void * n_control, void * model);
RcppExport SEXP _JANE_update_U(SEXP USEXP, SEXP ASEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP prob_matrixSEXP, SEXP betaSEXP, SEXP XSEXP, SEXP n_controlSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< void * >::type X(XSEXP);
    Rcpp::traits::input_parameter< void * >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< void * >::type model(modelSEXP);
    update_U(U, A, mus, omegas, prob_matrix, beta, X, n_control, model);
    return R_NilValue;
END_RCPP
}
// update_U_CC
void update_U_CC(arma::mat& U, double n_control, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, void * X, void * model);
RcppExport SEXP _JANE_update_U_CC(SEXP USEXP, SEXP n_controlSEXP, SEXP ASEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP prob_matrixSEXP, SEXP betaSEXP, SEXP XSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< void * >::type X(XSEXP);
    Rcpp::traits::input_parameter< void * >::type model(modelSEXP);
    update_U_CC(U, n_control, A, mus, omegas, prob_matrix, beta, X, model);
    return R_NilValue;
END_RCPP
}
// update_U_RE
void update_U_RE(arma::mat& U, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, arma::mat X, Rcpp::String model, void * n_control);
RcppExport SEXP _JANE_update_U_RE(SEXP USEXP, SEXP ASEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP prob_matrixSEXP, SEXP betaSEXP, SEXP XSEXP, SEXP modelSEXP, SEXP n_controlSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type model(modelSEXP);
    Rcpp::traits::input_parameter< void * >::type n_control(n_controlSEXP);
    update_U_RE(U, A, mus, omegas, prob_matrix, beta, X, model, n_control);
    return R_NilValue;
END_RCPP
}
// update_U_RE_CC
void update_U_RE_CC(arma::mat& U, double n_control, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, arma::mat X, Rcpp::String model);
RcppExport SEXP _JANE_update_U_RE_CC(SEXP USEXP, SEXP n_controlSEXP, SEXP ASEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP prob_matrixSEXP, SEXP betaSEXP, SEXP XSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type model(modelSEXP);
    update_U_RE_CC(U, n_control, A, mus, omegas, prob_matrix, beta, X, model);
    return R_NilValue;
END_RCPP
}
// update_beta
void update_beta(arma::colvec& beta, arma::sp_mat A, arma::mat U, double f, double e, void * X, void * n_control, void * model);
RcppExport SEXP _JANE_update_beta(SEXP betaSEXP, SEXP ASEXP, SEXP USEXP, SEXP fSEXP, SEXP eSEXP, SEXP XSEXP, SEXP n_controlSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< void * >::type X(XSEXP);
    Rcpp::traits::input_parameter< void * >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< void * >::type model(modelSEXP);
    update_beta(beta, A, U, f, e, X, n_control, model);
    return R_NilValue;
END_RCPP
}
// update_beta_CC
void update_beta_CC(arma::colvec& beta, arma::sp_mat A, double n_control, arma::mat U, double f, double e, void * X, void * model);
RcppExport SEXP _JANE_update_beta_CC(SEXP betaSEXP, SEXP ASEXP, SEXP n_controlSEXP, SEXP USEXP, SEXP fSEXP, SEXP eSEXP, SEXP XSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< void * >::type X(XSEXP);
    Rcpp::traits::input_parameter< void * >::type model(modelSEXP);
    update_beta_CC(beta, A, n_control, U, f, e, X, model);
    return R_NilValue;
END_RCPP
}
// update_beta_RE
void update_beta_RE(arma::colvec& beta, arma::sp_mat A, arma::mat U, arma::mat f, arma::colvec e, arma::mat X, Rcpp::String model, void * n_control);
RcppExport SEXP _JANE_update_beta_RE(SEXP betaSEXP, SEXP ASEXP, SEXP USEXP, SEXP fSEXP, SEXP eSEXP, SEXP XSEXP, SEXP modelSEXP, SEXP n_controlSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type model(modelSEXP);
    Rcpp::traits::input_parameter< void * >::type n_control(n_controlSEXP);
    update_beta_RE(beta, A, U, f, e, X, model, n_control);
    return R_NilValue;
END_RCPP
}
// update_beta_RE_CC
void update_beta_RE_CC(arma::colvec& beta, arma::sp_mat A, double n_control, arma::mat U, arma::mat f, arma::colvec e, arma::mat X, Rcpp::String model);
RcppExport SEXP _JANE_update_beta_RE_CC(SEXP betaSEXP, SEXP ASEXP, SEXP n_controlSEXP, SEXP USEXP, SEXP fSEXP, SEXP eSEXP, SEXP XSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type n_control(n_controlSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type model(modelSEXP);
    update_beta_RE_CC(beta, A, n_control, U, f, e, X, model);
    return R_NilValue;
END_RCPP
}
// update_mus_omegas
void update_mus_omegas(arma::mat prob_matrix, arma::mat U, double b, arma::rowvec a, double c, arma::mat G, arma::mat& mus, arma::cube& omegas);
RcppExport SEXP _JANE_update_mus_omegas(SEXP prob_matrixSEXP, SEXP USEXP, SEXP bSEXP, SEXP aSEXP, SEXP cSEXP, SEXP GSEXP, SEXP musSEXP, SEXP omegasSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type omegas(omegasSEXP);
    update_mus_omegas(prob_matrix, U, b, a, c, G, mus, omegas);
    return R_NilValue;
END_RCPP
}
// update_p
void update_p(arma::mat prob_matrix, arma::colvec& p, arma::colvec nu);
RcppExport SEXP _JANE_update_p(SEXP prob_matrixSEXP, SEXP pSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nu(nuSEXP);
    update_p(prob_matrix, p, nu);
    return R_NilValue;
END_RCPP
}
// update_prob_matrix_DA
void update_prob_matrix_DA(arma::mat& prob_matrix, arma::mat mus, arma::cube omegas, arma::colvec p, arma::mat U, double temp_beta);
RcppExport SEXP _JANE_update_prob_matrix_DA(SEXP prob_matrixSEXP, SEXP musSEXP, SEXP omegasSEXP, SEXP pSEXP, SEXP USEXP, SEXP temp_betaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type prob_matrix(prob_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type temp_beta(temp_betaSEXP);
    update_prob_matrix_DA(prob_matrix, mus, omegas, p, U, temp_beta);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_JANE_BIC_logit_NDH", (DL_FUNC) &_JANE_BIC_logit_NDH, 2},
    {"_JANE_BIC_logit_RS", (DL_FUNC) &_JANE_BIC_logit_RS, 2},
    {"_JANE_BIC_logit_RSR", (DL_FUNC) &_JANE_BIC_logit_RSR, 2},
    {"_JANE_BIC_ICL_MBC", (DL_FUNC) &_JANE_BIC_ICL_MBC, 1},
    {"_JANE_log_like_C", (DL_FUNC) &_JANE_log_like_C, 3},
    {"_JANE_gradient_C", (DL_FUNC) &_JANE_gradient_C, 3},
    {"_JANE_compute_dist", (DL_FUNC) &_JANE_compute_dist, 6},
    {"_JANE_log_Q", (DL_FUNC) &_JANE_log_Q, 17},
    {"_JANE_log_Q_RE", (DL_FUNC) &_JANE_log_Q_RE, 17},
    {"_JANE_draw_A_NDH_c", (DL_FUNC) &_JANE_draw_A_NDH_c, 2},
    {"_JANE_draw_A_RS_c", (DL_FUNC) &_JANE_draw_A_RS_c, 3},
    {"_JANE_draw_A_RSR_c", (DL_FUNC) &_JANE_draw_A_RSR_c, 4},
    {"_JANE_update_U", (DL_FUNC) &_JANE_update_U, 9},
    {"_JANE_update_U_CC", (DL_FUNC) &_JANE_update_U_CC, 9},
    {"_JANE_update_U_RE", (DL_FUNC) &_JANE_update_U_RE, 9},
    {"_JANE_update_U_RE_CC", (DL_FUNC) &_JANE_update_U_RE_CC, 9},
    {"_JANE_update_beta", (DL_FUNC) &_JANE_update_beta, 8},
    {"_JANE_update_beta_CC", (DL_FUNC) &_JANE_update_beta_CC, 8},
    {"_JANE_update_beta_RE", (DL_FUNC) &_JANE_update_beta_RE, 8},
    {"_JANE_update_beta_RE_CC", (DL_FUNC) &_JANE_update_beta_RE_CC, 8},
    {"_JANE_update_mus_omegas", (DL_FUNC) &_JANE_update_mus_omegas, 8},
    {"_JANE_update_p", (DL_FUNC) &_JANE_update_p, 3},
    {"_JANE_update_prob_matrix_DA", (DL_FUNC) &_JANE_update_prob_matrix_DA, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_JANE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
