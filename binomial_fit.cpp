#include <TMB.hpp>
#include <fenv.h>
using namespace density;
using namespace Eigen;
using namespace R_inla;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Number of observations.
  DATA_INTEGER(n);
  // Binomial response.
  DATA_MATRIX(y);
  // Number of species.
  DATA_INTEGER(n_species);
  // Number of trials.
  DATA_IVECTOR(n_trials);
  // Number of coefficients.
  DATA_INTEGER(n_betas);
  // Design matrix.
  DATA_MATRIX(mat);
  // SSTs.
  DATA_VECTOR(ssts);
  // Centred SSTs.
  DATA_VECTOR(ssts_centred);
  // Centred cosfiltered SSTs.
  DATA_VECTOR(ssts_cf_centred);
  // Julian month in radians for each observation.
  DATA_VECTOR(jmonth_rad);
  // Julian month in radians for each month.
  DATA_VECTOR(month_jmonth_rad);
  // Visitation probabilities.
  DATA_VECTOR(v);
  // Month ID for each observation.
  DATA_IVECTOR(month_id);
  // Mesh node ID for each observation.
  DATA_IVECTOR(mesh_id);
  // Number of months.
  DATA_INTEGER(n_months);
  // Number of mesh nodes.
  DATA_INTEGER(n_meshnodes);
  // Month temperature.
  DATA_VECTOR(month_temp_centred);
  // Cosfiltered month temperature.
  DATA_VECTOR(cosfilt_temp_centred);
  // Month-based design matrix.
  DATA_MATRIX(mat_pred);
  // Something something SPDE (needs demystifying).
  DATA_STRUCT(spde,spde_t);
  // Indicators for latent variables.
  DATA_INTEGER(fit_psi);
  DATA_INTEGER(fit_omega);
  DATA_INTEGER(fit_epsilon);
  DATA_INTEGER(fit_int);
  DATA_INTEGER(fit_cf);
  // Vector of coefficients.
  PARAMETER_MATRIX(betas);
  // Parameters for the AR(1) psi process.
  PARAMETER_VECTOR(link_phi_psi);
  PARAMETER_VECTOR(log_sigma_psi);
  // Parameters for the omega spatial process.
  PARAMETER_VECTOR(log_kappa_omega);
  PARAMETER_VECTOR(log_sigma_omega);
  // Parameters for the AR(1) part of the epsilon process.
  PARAMETER_VECTOR(link_phi_epsilon);
  PARAMETER_VECTOR(log_sigma_epsilon);
  // Parameter for the spatial GMRF of the epsilon process.
  PARAMETER_VECTOR(log_kappa_epsilon);
  // Parameters for the SST-interaction spatial field.
  PARAMETER_VECTOR(log_kappa_u_int);
  PARAMETER_VECTOR(log_tau_u_int);
  // Parameters for the cosfilt-interaction spatial field.
  PARAMETER_VECTOR(log_kappa_u_cf);
  PARAMETER_VECTOR(log_tau_u_cf);
  // Horizontal shift for seasonal interaction term.
  PARAMETER_VECTOR(link_gamma);
  // Vector of temporal process random variables.
  PARAMETER_MATRIX(psi_t_all);
  // Array of omega spatial process random variables.
  PARAMETER_MATRIX(omega_s_all);
  // Array of spatiotemporal field random variables.
  PARAMETER_ARRAY(epsilon_st_all);
  // Vector of SST-interaction spatial field.
  PARAMETER_MATRIX(u_int_all);
  // Vector of cosfilt-interaction spatial field.
  PARAMETER_MATRIX(u_cf_all);
  // Transforming parameters.
  vector<Type> phi_psi = 2*exp(link_phi_psi)/(1 + exp(link_phi_psi)) - 1;
  vector<Type> sigma_psi = exp(log_sigma_psi);
  vector<Type> kappa_omega = exp(log_kappa_omega);
  vector<Type> sigma_omega = exp(log_sigma_omega);
  vector<Type> phi_epsilon = 2*exp(link_phi_epsilon)/(1 + exp(link_phi_epsilon)) - 1;
  vector<Type> sigma_epsilon = exp(log_sigma_epsilon);
  vector<Type> kappa_epsilon = exp(log_kappa_epsilon);
  vector<Type> kappa_u_int = exp(log_kappa_u_int);
  vector<Type> tau_u_int = exp(log_tau_u_int);
  vector<Type> kappa_u_cf = exp(log_kappa_u_cf);
  vector<Type> tau_u_cf = exp(log_tau_u_cf);
  vector<Type> gamma = 3.141593*exp(link_gamma)/(1 + exp(link_gamma));
  ADREPORT(phi_psi);
  ADREPORT(sigma_psi);
  ADREPORT(kappa_omega);
  ADREPORT(sigma_omega);
  ADREPORT(phi_epsilon);
  ADREPORT(sigma_epsilon);
  ADREPORT(kappa_epsilon);
  ADREPORT(kappa_u_int);
  ADREPORT(tau_u_int);
  ADREPORT(kappa_u_cf);
  ADREPORT(tau_u_cf);
  ADREPORT(gamma);
  vector<Type> f_all(n_species);
  array<Type> d_full_logit(n_species,n_meshnodes,n_months);
  for (int s = 0; s < n_species; s++){
    // Extracting species-specific stuff.
    vector<Type> y_s(n);
    y_s = y.col(s);
    vector<Type> betas_s(n_betas);
    betas_s = betas.row(s);
    // Filling latent variables.
    vector<Type> psi_t(n_months);
    vector<Type> omega_s(n_meshnodes);
    array<Type> epsilon_st(n_meshnodes, n_months);
    vector<Type> u_int(n_meshnodes);
    vector<Type> u_cf(n_meshnodes);
    for (int i = 0; i < n_meshnodes; i++){
      omega_s(i) = omega_s_all(s, i);
      u_int(i) = u_int_all(s, i);
      u_cf(i) = u_cf_all(s, i);
      for (int j = 0; j < n_months; j++){
	if (i == 0){
	  psi_t(j) = psi_t_all(s, j);
	}
	epsilon_st(i, j) = epsilon_st_all(s, i, j);
      }
    }
    // Calculating fitted probabilities.
    vector<Type> d_fixed_logit(n);
    vector<Type> d2(n);
    vector<Type> d_fixed_logit_pred(n_months);
    vector<Type> p(n);
    d_fixed_logit = mat*betas_s;
    for (int i = 0; i < n; i++){
      d2(i) = d_fixed_logit(i) + psi_t(month_id(i)) + omega_s(mesh_id(i)) + epsilon_st(mesh_id(i), month_id(i));
      // Adding contribution from u_int.
      if (fit_int == 1){
	d2(i) += ssts_centred(i)*u_int(mesh_id(i))/tau_u_int(s);
      } else if (fit_int == 2){
	d2(i) += cos(jmonth_rad(i) - gamma(s))*u_int(mesh_id(i))/tau_u_int(s);
      }
      // Adding contribution from u_cf.
      if (fit_cf == 1){
	d2(i) += ssts_cf_centred(i)*u_cf(mesh_id(i))/tau_u_cf(s);
      }
    }
    d_fixed_logit_pred = mat_pred*betas_s;
    for (int i = 0; i < n_meshnodes; i++){
      for (int j = 0; j < n_months; j++){
	d_full_logit(s,i,j) = d_fixed_logit_pred(j) + psi_t(i) + omega_s(i) + epsilon_st(i, j);
	// Adding contribution from u_int.
	if (fit_int == 1){
	  d_full_logit(s,i,j) += month_temp_centred(j)*u_int(i)/tau_u_int(s);
	} else if (fit_int == 2){
	  d_full_logit(s,i,j) += cos(month_jmonth_rad(j) - gamma(s))*u_int(i)/tau_u_int(s);
	}
	// Adding contribution from u_cf.
	if (fit_cf == 1){
	  d_full_logit(s,i,j) += cosfilt_temp_centred(j)*u_cf(i)/tau_u_cf(s);
	}
      }
    }
    p = v*exp(d2)/(1 + exp(d2));    
    Type dummy_y;
    Type dummy_n;
    f_all(s) = 0;
    // Component due to f(y | u).
    for (int i = 0; i < n; i++){
      dummy_y = y_s(i);
      dummy_n = n_trials(i);
      f_all(s) -= dbinom(dummy_y, dummy_n, p(i), true);
    }
    // Component due to temopral psi field.
    if (fit_psi == 1){
      f_all(s) += SCALE(AR1(phi_psi(s)),sigma_psi(s))(psi_t);
    }
    // Component due to spatial omega field.
    if (fit_omega == 1){
      SparseMatrix<Type> Q_omega = Q_spde(spde,kappa_omega(s));
      f_all(s) += SCALE(GMRF(Q_omega), sigma_omega(s))(omega_s);
    }
    // Component due to spatiotemporal epsilon field.
    if (fit_epsilon == 1){
      SparseMatrix<Type> Q_epsilon = Q_spde(spde,kappa_epsilon(s));
      f_all(s) += SEPARABLE(SCALE(AR1(phi_epsilon(s)),sigma_epsilon(s)), GMRF(Q_epsilon))(epsilon_st);
    }
    // Component due to SST-interaction spatial field.
    if (fit_int > 0){
      SparseMatrix<Type> Q_int = Q_spde(spde,kappa_u_int(s));
      f_all(s) += GMRF(Q_int)(u_int);
    }
    // Component due to cosfilt-interaction spatial field.
    if (fit_cf == 1){
      SparseMatrix<Type> Q_cf = Q_spde(spde,kappa_u_cf(s));
      f_all(s) += GMRF(Q_cf)(u_cf);
    }
  }
  REPORT(f_all);
  REPORT(d_full_logit);
  // Returning negative of the joint density.
  Type f = sum(f_all);
  return f;
}
