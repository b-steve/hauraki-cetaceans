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
  DATA_INTEGER(fit_epsilon);
  DATA_INTEGER(fit_omega);
  DATA_INTEGER(fit_int);
  // Number of latent factor fields.
  DATA_INTEGER(n_factors);
  // Vector of coefficients.
  PARAMETER_MATRIX(betas);
  // Parameters for the latent factor fields.
  PARAMETER_VECTOR(phi);
  PARAMETER_VECTOR(link_rho);
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(log_kappa_epsilon);
  PARAMETER_VECTOR(log_kappa_omega);
  PARAMETER_VECTOR(log_lambda);
  // Parameters for the SST-interaction spatial field.
  PARAMETER_VECTOR(log_kappa_u_int);
  PARAMETER_VECTOR(log_tau_u_int);
  // Horizontal shift for seasonal interaction term.
  PARAMETER_VECTOR(link_gamma);
  // Values in loading matrix.
  PARAMETER_VECTOR(L_val);
  // Array of epsilon-field values.
  PARAMETER_ARRAY(epsilon_input);
  // Array of omega-field values.
  PARAMETER_MATRIX(omega_input);
  // Vector of SST-interaction spatial field.
  PARAMETER_MATRIX(u_int_all);
  // Transforming parameters.
  vector<Type> rho = 2*exp(link_rho)/(1 + exp(link_rho)) - 1;
  vector<Type> kappa_epsilon = exp(log_kappa_epsilon);
  vector<Type> kappa_omega = exp(log_kappa_omega);
  vector<Type> lambda = exp(log_lambda);
  vector<Type> kappa_u_int = exp(log_kappa_u_int);
  vector<Type> tau_u_int = exp(log_tau_u_int);
  vector<Type> gamma = 3.141593*exp(link_gamma)/(1 + exp(link_gamma));
  ADREPORT(rho);
  ADREPORT(kappa_epsilon);
  ADREPORT(kappa_omega);
  ADREPORT(lambda);
  ADREPORT(kappa_u_int);
  ADREPORT(tau_u_int);
  ADREPORT(gamma);
  vector<Type> f_species(n_species);
  vector<Type> f_factors(n_factors);
  array<Type> d_full_logit(n_species,n_meshnodes,n_months);
  // Assembling the loadings matrix.
  matrix<Type> L_mat(n_species, n_factors);
  int k = 0;
  for (int i = 0; i < n_factors; i++){
    for (int j = 0; j < n_species; j++){
      if (j >= i){
	L_mat(j, i) = L_val(k);
	k++;
      } else {
	L_mat(j, i) = 0.0;
      }
    }
  }
  // Can rotate the L_mat matrix as per supplementary materials of
  // Thorson et al. (2015), but yet to implement this yet.

  // Calculating tau for omega and epsilon fields.
  vector<Type> log_tau_epsilon(n_factors);
  vector<Type> log_tau_omega(n_factors);
  for (int k = 0; k < n_factors; k++){
    log_tau_omega(k) = 0.5*(log(1.0 + lambda(k)) - log(4.0*M_PI*exp(2.0*log_kappa_omega(k))*pow(1 - rho(k), 2)));
    log_tau_epsilon(k) = log_tau_omega(k) + log(1 - rho(k)) + log_kappa_omega(k) - log_kappa_epsilon(k) - 0.5*log(lambda(k)*(1 - pow(rho(k), 2)));
  }
  vector<Type> tau_epsilon = exp(log_tau_epsilon);
  vector<Type> tau_omega = exp(log_tau_omega);
  // Assembling the latent factors.
  array<Type> epsilon(n_factors, n_meshnodes, n_months);
  matrix<Type> omega(n_factors, n_meshnodes);
  array<Type> psi_st(n_factors, n_meshnodes, n_months);
  for (int k = 0; k < n_factors; k++){
    for (int i = 0; i < n_meshnodes; i++){
      omega(k, i) = omega_input(k, i)/tau_omega(k);
      for (int j = 0; j < n_months; j++){
	// Constructing the latent epsilon field from inputs.
	epsilon(k, i, j) = epsilon_input(k, i, j)/tau_epsilon(k);
	// Don't really know what is going on here, but I think it is
	// some clever construction of the latent factor fields that
	// helps with computational efficiency. Something or other
	// about Gompertz processes?
	psi_st(k, i, j) = phi(k)*pow(rho(k), j) + // Temporal autocorrelation term.
	  epsilon(k, i, j) + // Spatiotemporal residual term.
	  alpha(k)/(1 - rho(k)) + // Unknown term I don't understand (diff in height for each field?)
	  omega(k, i)/(1 - rho(k)); // Spatial variance term.
      }
    }
  }
  for (int s = 0; s < n_species; s++){
    // Extracting species-specific stuff.
    vector<Type> y_s(n);
    y_s = y.col(s);
    vector<Type> betas_s(n_betas);
    betas_s = betas.row(s);
    // Filling latent variables for species s.
    array<Type> u_st(n_meshnodes, n_months);
    vector<Type> u_int(n_meshnodes);
    for (int i = 0; i < n_meshnodes; i++){
      u_int(i) = u_int_all(s, i);
      for (int j = 0; j < n_months; j++){
	u_st(i, j) = 0.0;
	for (int k = 0; k < n_factors; k++){
	  u_st(i, j) += psi_st(k, i, j)*L_mat(s, k);
	}
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
      // Adding contribution from u_int.
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
    f_species(s) = 0;
    // Component due to f(y | u).
    for (int i = 0; i < n; i++){
      dummy_y = y_s(i);
      dummy_n = n_trials(i);
      f_species(s) -= dbinom(dummy_y, dummy_n, p(i), true);
    }
    // Component due to interaction spatial field.
    if (fit_int > 0){
      SparseMatrix<Type> Q_int = Q_spde(spde, kappa_u_int(s));
      f_species(s) += GMRF(Q_int)(u_int);
    }
  }
  // Component due to spatiotemporal factor fields.
  if (fit_epsilon + fit_omega > 0){
    array<Type> epsilon_tmp(n_meshnodes, n_months);
    SparseMatrix<Type> Q;
    for (int k = 0; k < n_factors; k++){
      f_factors(k) = 0.0;
      // For the omega field.
      if (fit_omega == 1){
	Q = Q_spde(spde, kappa_omega(k));
	f_factors(k) += GMRF(Q)(omega_input.row(k));
      }
      // For the epsilon field.
      if (fit_epsilon){
	for (int i = 0; i < n_meshnodes; i++){
	  for (int j = 0; j < n_months; j++){
	    epsilon_tmp(i, j) = epsilon_input(k, i, j);
	  }
	}
	SparseMatrix<Type> Q = Q_spde(spde, kappa_epsilon(k));
	f_factors(k) += SEPARABLE(AR1(rho(k)), GMRF(Q))(epsilon_tmp);
      }
    }
  }
  REPORT(d_full_logit);
  // Returning negative of the joint density.
  Type f = sum(f_species) + sum(f_factors);
  return f;
}
