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
  // Month-based design matrix.
  DATA_MATRIX(mat_pred);
  // Something something SPDE (needs demystifying).
  DATA_STRUCT(spde,spde_t);
  // Indicators for spatial fields.
  DATA_INTEGER(fit_st);
  DATA_INTEGER(fit_int);
  // Vector of coefficients.
  PARAMETER_MATRIX(betas);
  // Parameters for the AR(1) process.
  PARAMETER_VECTOR(link_phi);
  PARAMETER_VECTOR(log_sigma_u_t);
  // Parameter for the spatial GMRF.
  PARAMETER_VECTOR(log_kappa_u_s);
  // Parameters for the SST-interaction spatial field.
  PARAMETER_VECTOR(log_kappa_u_int);
  PARAMETER_VECTOR(log_tau_u_int);
  // Horizontal shift for seasonal interaction term.
  PARAMETER_VECTOR(link_gamma);
  // Array of spatiotemporal field random variables.
  PARAMETER_ARRAY(u_st_all);
  // Vector of SST-interaction spatial field.
  PARAMETER_MATRIX(u_int_all);
  // Transforming parameters.
  vector<Type> phi = 2*exp(link_phi)/(1 + exp(link_phi)) - 1;
  vector<Type> sigma_u_t = exp(log_sigma_u_t);
  vector<Type> kappa_u_s = exp(log_kappa_u_s);
  vector<Type> kappa_u_int = exp(log_kappa_u_int);
  vector<Type> tau_u_int = exp(log_tau_u_int);
  vector<Type> gamma = 3.141593*exp(link_gamma)/(1 + exp(link_gamma));
  ADREPORT(phi);
  ADREPORT(sigma_u_t);
  ADREPORT(kappa_u_s);
  ADREPORT(kappa_u_int);
  ADREPORT(tau_u_int);
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
    array<Type> u_st(n_meshnodes, n_months);
    vector<Type> u_int(n_meshnodes);
    for (int i = 0; i < n_meshnodes; i++){
      u_int(i) = u_int_all(s, i);
      for (int j = 0; j < n_months; j++){
	u_st(i, j) = u_st_all(s, i, j);
      }
    }
    // Calculating fitted probabilities.
    vector<Type> d_fixed_logit(n);
    vector<Type> d2(n);
    vector<Type> d_fixed_logit_pred(n_months);
    vector<Type> p(n);
    d_fixed_logit = mat*betas_s;
    for (int i = 0; i < n; i++){
      d2(i) = d_fixed_logit(i) + u_st(mesh_id(i), month_id(i));
      // Addint contribution from u_int.
      if (fit_int == 1){
	d2(i) += ssts_centred(i)*u_int(mesh_id(i))/tau_u_int(s);
      } else if (fit_int == 2){
	d2(i) += cos(jmonth_rad(i) - gamma(s))*u_int(mesh_id(i))/tau_u_int(s);
      }
    }
    d_fixed_logit_pred = mat_pred*betas_s;
    for (int i = 0; i < n_meshnodes; i++){
      for (int j = 0; j < n_months; j++){
	d_full_logit(s,i,j) = d_fixed_logit_pred(j) + u_st(i, j);
	// Adding contribution from u_int.
	if (fit_int == 1){
	  d_full_logit(s,i,j) += month_temp_centred(j)*u_int(i)/tau_u_int(s);
	} else if (fit_int == 2){
	  d_full_logit(s,i,j) += cos(month_jmonth_rad(j) - gamma(s))*u_int(i)/tau_u_int(s);
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
    // Component due to spatiotemporal field.
    if (fit_st == 1){
      SparseMatrix<Type> Q = Q_spde(spde,kappa_u_s(s));
      f_all(s) += SEPARABLE(SCALE(AR1(phi(s)),sigma_u_t(s)), GMRF(Q))(u_st);
    }
    // Component due to SST-interaction spatial field.
    if (fit_int > 0){
      SparseMatrix<Type> Q_int = Q_spde(spde,kappa_u_int(s));
      f_all(s) += GMRF(Q_int)(u_int);
    }
  }
  REPORT(f_all);
  REPORT(d_full_logit);
  // Returning negative of the joint density.
  Type f = sum(f_all);
  return f;
}
