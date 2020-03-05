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
  DATA_IVECTOR(y);
  // Number of trials.
  DATA_IVECTOR(n_trials);
  // Design matrix.
  DATA_MATRIX(mat);
  // SSTs.
  DATA_VECTOR(ssts);
  // Centred SSTs.
  DATA_VECTOR(ssts_centred);
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
  PARAMETER_VECTOR(betas);
  // Parameters for the AR(1) process.
  PARAMETER(link_phi);
  PARAMETER(log_sigma_u_t);
  // Parameter for the spatial GMRF.
  PARAMETER(log_kappa_u_s);
  // Parameters for the SST-interaction spatial field.
  PARAMETER(log_kappa_u_int);
  PARAMETER(log_tau_u_int);
  // Array of spatiotemporal field random variables.
  PARAMETER_ARRAY(u_st);
  // Vector of SST-interaction spatial field.
  PARAMETER_VECTOR(u_int);
  // Transforming parameters.
  Type phi = 2*exp(link_phi)/(1 + exp(link_phi)) - 1;
  Type sigma_u_t = exp(log_sigma_u_t);
  Type kappa_u_s = exp(log_kappa_u_s);
  Type kappa_u_int = exp(log_kappa_u_int);
  Type tau_u_int = exp(log_tau_u_int);
  ADREPORT(phi);
  ADREPORT(sigma_u_t);
  ADREPORT(kappa_u_s);
  ADREPORT(kappa_u_int);
  // Calculating fitted probabilities.
  vector<Type> d_fixed_logit(n);
  vector<Type> d2(n);
  vector<Type> d_fixed_logit_pred(n_months);
  matrix<Type> d_full_logit(n_meshnodes,n_months);
  vector<Type> p(n);
  d_fixed_logit = mat*betas;
  for (int i = 0; i < n; i++){
    d2(i) = d_fixed_logit(i) + u_st(mesh_id(i), month_id(i)) + ssts_centred(i)*u_int(mesh_id(i))/tau_u_int;
  }
  d_fixed_logit_pred = mat_pred*betas;
  for (int i = 0; i < n_meshnodes; i++){
    for (int j = 0; j < n_months; j++){
      d_full_logit(i,j) = d_fixed_logit_pred(j) + u_st(i, j) + month_temp_centred(j)*u_int(i)/tau_u_int;
    }
  }
  REPORT(d_full_logit);
  p = v*exp(d2)/(1 + exp(d2));    
  Type f = 0;
  Type dummy_y;
  Type dummy_n;
  // Component due to f(y | u).
  for (int i = 0; i < n; i++){
    dummy_y = y(i);
    dummy_n = n_trials(i);
    f -= dbinom(dummy_y, dummy_n, p(i), true);
  }
  // Component due to spatiotemporal field.
  if (fit_st == 1){
    SparseMatrix<Type> Q = Q_spde(spde,kappa_u_s);
    f += SEPARABLE(SCALE(AR1(phi),sigma_u_t), GMRF(Q))(u_st);
  }
  // Component due to SST-interaction spatial field.
  if (fit_int == 1){
    SparseMatrix<Type> Q_int = Q_spde(spde,kappa_u_int);
    f += GMRF(Q_int)(u_int);
  }
  // Returning negative of the joint density.
  return f;
}
