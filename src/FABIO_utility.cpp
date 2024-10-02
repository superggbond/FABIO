#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;
// using namespace arma;

// [[Rcpp::export]]
void SampleZ(const arma::vec& y, const arma::vec& z_hat, arma::vec& z){
  double d1, d2, z_rand = 0.0;
  for (int i = 0; i < z.size(); ++i){
    d1 = y[i];
    d2 = z_hat[i];
    // y is centered for case control studies
    if (d1 <= 0){
      // Control, right truncated
      z_rand = d2 + rnorm(1, 0, 1)[0];
      while (z_rand > 0){z_rand = d2 + rnorm(1, 0, 1)[0];}
    } else {
      z_rand = d2 + rnorm(1, 0, 1)[0];
      while (z_rand < 0){z_rand = d2 + rnorm(1, 0, 1)[0];}
    }
    z[i] = z_rand;
  }
  return;
}

// [[Rcpp::export]]
double CalcPosterior(List cHyp, const double ng_test, const double yty,
                    const arma::mat& Xgamma, const int s_size, const arma::mat& XtX,
                    const arma::vec& Xty, arma::vec& Xb, arma::vec& beta) {
  double sigma_a2 = (double)cHyp["h"] / ((1.0 - (double)cHyp["h"])
                                           * std::exp((double)cHyp["logp"]) * (double)ng_test);
  double logpost = 0.0;
  double P_yy = yty;

  const arma::mat Xgamma_sub = Xgamma.cols(0,(s_size-1));
  const arma::mat XtX_sub = XtX.submat(0,0,(s_size-1),(s_size-1));
  const arma::vec Xty_sub = Xty.subvec(0,(s_size-1));

  // Calculate Omega
  arma::mat Omega = sigma_a2 * XtX_sub;
  arma::mat M_temp = arma::eye(s_size,s_size);
  Omega = Omega + M_temp;

  // Calculate beta_hat
  double epsilon = 1e-8;
  arma::mat Omega_L;
  bool chol_success = chol(Omega_L, Omega, "lower");

  while (!chol_success) {
    Omega.diag() += epsilon;
    chol_success = chol(Omega_L, Omega, "lower");
    epsilon *= 10;
  }

  // mat Omega_L = chol(Omega, "lower");
  arma::vec logdet_O = Omega_L.diag();
  logdet_O = log(logdet_O);
  double logdet_O_sum = 2.0 * sum(logdet_O);
  arma::vec beta_hat = solve(Omega, Xty_sub);
  beta_hat = sigma_a2 * beta_hat;

  arma::mat d_temp = Xty_sub.t() * beta_hat;
  double d = d_temp(0,0);
  P_yy = P_yy - d;

  // tau = 1 for probit plink
  double tau = 1.0;

  // sample beta
  for (int i=0; i < s_size; ++i) {
    beta(i) = rnorm(1, 0, 1)[0];
  }
  arma::vec beta_sub = beta.subvec(0,(s_size-1));
  arma::mat Omega_new = diagmat(Omega);
  beta_sub = solve(Omega_new, M_temp) * beta_sub;

  // This computes inv(L^T(Omega)) %*% beta
  beta_sub = beta_sub * std::sqrt(sigma_a2 / tau);
  beta_sub = beta_sub + beta_hat;
  Xb = Xgamma_sub * beta_sub;

  logpost = -0.5 * logdet_O_sum;
  logpost = logpost - 0.5 * P_yy;
  logpost = logpost + ((double)cHyp["n_gamma"] - 1.0) * (double)cHyp["logp"]
  + ((double)ng_test - (double)cHyp["n_gamma"]) * log(1.0 - std::exp((double)cHyp["logp"]));

  // update beta
  for (int i=0; i < s_size; ++i) {
    beta[i] = beta_sub[i];
  }

  return logpost;
}

// [[Rcpp::export]]
void CalcXtX(const arma::mat& X, const arma::vec& y, const int s_size, arma::mat& XtX, arma::vec& Xty) {
  const arma::mat X_sub = X.cols(0,(s_size-1));
  XtX.submat(0,0,(s_size-1),(s_size-1)) = X_sub.t() * X_sub;
  Xty.subvec(0,(s_size-1)) = X_sub.t() * y;
  return;
}

// [[Rcpp::export]]
void CalcCC_PVEnZ(const arma::vec& Xb, const double ni_test, List& cHyp, arma::vec& z_hat) {
  arma::mat d = Xb.t() * Xb;
  double temp = d(0)/ni_test;
  temp /= temp + 1.0;
  cHyp["pve"] = temp;
  cHyp["pge"] = 1.0;
  z_hat = Xb;
  return;
}

// [[Rcpp::export]]
List ProposeH(const List& cHyp_old, int rep, const double h_max, const double h_min, const double h_scale) {
  List cHyp_new = clone(cHyp_old);
  double h = cHyp_old["h"];
  double d_h = (h_max - h_min) * h_scale;

  for (int i=0; i<rep; ++i){
    double ran = runif(1, 0, 1)[0];
    h += (ran - 0.5) * d_h;
    if (h < h_min) {h = 2 * h_min - h;}
    if (h > h_max) {h = 2 * h_max - h;}
  }
  cHyp_new["h"] = h;

  return cHyp_new;
}

// [[Rcpp::export]]
double ProposePi(const List& cHyp_old, List& cHyp_new, int rep, double logp_max, double logp_min, double logp_scale, double beta_a, double beta_b) {
  double logp_old = cHyp_old["logp"];
  double logp_new;
  double log_ratio = 0.0;
  double d_logp, temp;
  temp = (logp_max - logp_min) * logp_scale;
  d_logp = std::min(0.1, temp);
  
  if (beta_a==0 & beta_b==0){
    for (int i=0; i<rep; ++i){
      double ran = runif(1, 0, 1)[0];
      logp_new = logp_old + (ran - 0.5) * d_logp;
      if (logp_new < logp_min) {logp_new = 2 * logp_min - logp_new;}
      if (logp_new > logp_max) {logp_new = 2 * logp_max - logp_new;}
      log_ratio = log_ratio + logp_new - logp_old;
      logp_old = logp_new;
    }
  } else {
    double ran = rbeta(1,beta_a,beta_b)[0];
    logp_new = log(ran);
    if (logp_new < logp_min) {logp_new = logp_min;}
    if (logp_new > logp_max) {logp_new = logp_max;}
    log_ratio = log_ratio + logp_new - logp_old;
    logp_old = logp_new;
  }
  cHyp_new["logp"] = logp_new;
  return log_ratio;
}

// [[Rcpp::export]]
List mcmc_iter(const int total_step, const int w_step, const int r_pace, const int w_pace, const int n_mh, const int ng_test, const int ni_test,
               const double h_max, const double h_min, const double h_scale, const int g_max, const int g_min,
               const double logp_max, const double logp_min, const double logp_scale, const double beta_a, const double beta_b,
               const arma::vec& y, arma::vec& z_hat, arma::vec& z, arma::vec& rank_old, arma::vec& beta_old, arma::vec& beta_new,
               arma::vec& Xtz_old, arma::vec& Xtz_new, arma::vec& Xb_old, arma::vec& Xb_new, const arma::vec& p_gamma,
               const arma::ivec pos_vec, List& cHyp_old, const arma::mat& X, arma::mat& Xgamma_old, arma::mat& Xgamma_new, arma::mat& XtX_old, arma::mat& XtX_new,
               arma::mat& Result_hyp, arma::mat& Result_gamma, arma::mat& beta_g, bool display_progress=true){
  Function ProposeGamma("ProposeGamma");
  Function SetXgamma("SetXgamma");

  int rep;
  List cHyp_new;
  List ProposeGamma_results;
  List SetXgamma_results;
  arma::vec rank_new;
  arma::ivec pos_vec_selected;
  arma::uvec indices;
  double logPost_old;
  double logPost_new;
  double logMHratio;
  int accept;
  int n_accept = 0;
  int w = 0;
  int w_col;

  Progress p(total_step, display_progress);
  for (int t=0; t < total_step; ++t){
    p.increment();
    SampleZ(y, z_hat, z);
    double mean_z = mean(z);
    z = z - mean_z;
    double ztz = as_scalar(z.t() * z);

    // first proposal
    arma::mat Xold_sub = Xgamma_old.submat(0,0,(ni_test-1),(rank_old.size()-1));
    arma::vec Xtz_old_new = Xold_sub.t() * z;
    for (int i=0; i < rank_old.size(); ++i){Xtz_old(i) = Xtz_old_new(i);}
    logPost_old = CalcPosterior(cHyp_old,ng_test,ztz,Xgamma_old,rank_old.size(),XtX_old,Xtz_old,Xb_old,beta_old);

    // M-H steps
    for (int i=0; i < n_mh; ++i){
      if (runif(1, 0, 1)[0] < 0.33){
        rep = 1 + floor(runif(1, 0, 20)[0]);
      } else {
        rep = 1;
      }

      logMHratio = 0.0;
      // update h
      cHyp_new = ProposeH(cHyp_old, rep, h_max, h_min, h_scale);

      // update n_gamma
      ProposeGamma_results = ProposeGamma(rank_old, p_gamma, cHyp_old, cHyp_new, rep, g_max, g_min, ng_test);
      cHyp_new = as<List>(ProposeGamma_results["cHyp_new"]);
      rank_new = as<arma::vec>(ProposeGamma_results["rank_new"]);
      logMHratio = logMHratio + (double)ProposeGamma_results["logp"];

      // update pi
      logMHratio = logMHratio + ProposePi(cHyp_old, cHyp_new, rep, logp_max, logp_min, logp_scale, beta_a, beta_b);

      if ((int)cHyp_new["n_gamma"] <= 20 || (int)cHyp_old["n_gamma"] <= 20){
        indices = arma::uvec(rank_new.size());
        for (int j=0; j<rank_new.size(); ++j){indices[j]=rank_new[j]-1;}
        pos_vec_selected = pos_vec(indices);
        indices = arma::uvec(pos_vec_selected.size());
        for (int j=0; j<pos_vec_selected.size(); ++j){indices[j]=pos_vec_selected[j];}
        Xgamma_new.cols(0,(rank_new.size()-1)) = X.cols(indices);
        CalcXtX(Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
      } else {
        SetXgamma_results = SetXgamma(X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new, Xgamma_new, XtX_new, Xtz_new);
        Xgamma_new = as<arma::mat>(SetXgamma_results["X_new"]);
        XtX_new = as<arma::mat>(SetXgamma_results["XtX_new"]);
        Xtz_new = as<arma::vec>(SetXgamma_results["Xty_new"]);
      }
      logPost_new = CalcPosterior(cHyp_new,ng_test,ztz,Xgamma_new,rank_new.size(),XtX_new,Xtz_new,Xb_new,beta_new);
    }
    logMHratio += logPost_new - logPost_old;

    if (logMHratio > 0.0 || log(runif(1,0,1)[0]) < logMHratio) {
      accept = 1;
      n_accept += 1;
    } else {
      accept = 0;
    }

    if (accept == 1){
      logPost_old = logPost_new;
      cHyp_old = clone(cHyp_new);
      for (int j=0; j<Xb_old.size(); ++j){Xb_old[j]=Xb_new[j];}

      rank_old.set_size(rank_new.size());
      for (int j=0; j<rank_old.size(); ++j){rank_old[j]=rank_new[j];}
      Xgamma_old.cols(0,(rank_new.size()-1)) = Xgamma_new.cols(0,(rank_new.size()-1));
      XtX_old.submat(0,0,(rank_new.size()-1),(rank_new.size()-1)) = XtX_new.submat(0,0,(rank_new.size()-1),(rank_new.size()-1));
      Xtz_old.subvec(0,(rank_new.size()-1)) = Xtz_new.subvec(0,(rank_new.size()-1));
      beta_old.subvec(0,(rank_new.size()-1)) = beta_new.subvec(0,(rank_new.size()-1));
    } else {
      cHyp_new = clone(cHyp_old);
    }

    // Calculate z_hat, and pve.
    CalcCC_PVEnZ(Xb_old, ni_test, cHyp_old, z_hat);

    // Sample mu and update z_hat.
    z = z - z_hat;
    double mean_z_new = mean(z);
    z = z - mean_z_new;
    mean_z = mean_z + mean_z_new;
    double sd = sqrt(1.0 / ni_test);
    double new_value = R::rnorm(0, sd);
    mean_z = mean_z + new_value;
    z_hat = z_hat + mean_z;

    // Save data.
    if (t >= w_step) {
      if (t % r_pace == 0) {
        w_col = w % w_pace;

        Result_hyp(w_col,0) = cHyp_old["h"];
        Result_hyp(w_col,1) = cHyp_old["pve"];
        Result_hyp(w_col,2) = cHyp_old["pge"];
        Result_hyp(w_col,3) = cHyp_old["logp"];
        Result_hyp(w_col,4) = cHyp_old["n_gamma"];

        Result_gamma.row(w_col).zeros();
        for (int i=0; i<rank_old.size(); ++i) {
          int rank_ind = rank_old(i) - 1;
          int pos_g = pos_vec(rank_ind);
          Result_gamma(w_col, i) = pos_g + 1;

          beta_g(pos_g,0) = beta_g(pos_g,0) + beta_old(i);
          beta_g(pos_g,1) = beta_g(pos_g,1) + 1;
        }
        w++;
      }
    }
  }
  return List::create(_["w"]=w,
                      _["Result_hyp"]=Result_hyp,
                      _["Result_gamma"]=Result_gamma,
                      _["beta_g"]=beta_g);
}
