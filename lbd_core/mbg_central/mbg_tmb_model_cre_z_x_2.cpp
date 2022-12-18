// ///////////////////////////////////////////////////
// August 2020
// Template file for space-time-Z GPR model.
// Used for fitting IHME Geospatial MBG models
// ///////////////////////////////////////////////////

// ///////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. In our nomenclature, Z is a third interaction (ie age) which defaults to AR1
// 3. Requires same space mesh for all time-Z points
// 4. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 5. ref <<<< FILEPATH REDACTED >>>>
//        <<<< FILEPATH REDACTED >>>>
// ///////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;


// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = exp(2. * logkappa);
  Type kappa4 = kappa2*kappa2;
  Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Robust Inverse Logit that sets min and max values to avoid numerical instability
template<class Type>
Type invlogit_robust(Type x){
  if (x < -20.723){
    x = -20.723; // corresponds to p=1e-9
  } else if ( x > 20.723 ){
    x = 20.723;  // cooresponds to p=1-1e-9
  }
  return 1 / (1 + exp( -1.0 * x ));
}


// AR funtion from neal m
template<class Type>
SparseMatrix<Type> ar_Q(int N, Type rho, Type sigma) {
  SparseMatrix<Type> Q(N,N);
  Q.insert(0,0) = (1.) / pow(sigma, 2.);
  for (size_t n = 1; n < N; n++) {
    Q.insert(n,n) = (1. + pow(rho, 2.)) / pow(sigma, 2.);
    Q.insert(n-1,n) = (-1. * rho) / pow(sigma, 2.);
    Q.insert(n,n-1) = (-1. * rho) / pow(sigma, 2.);
  }
  Q.coeffRef(N-1,N-1) = (1.) / pow(sigma, 2.);
  return Q;
}

// Corresponding list object on the C++ side
template<class Type>
struct option_list {
  int use_priors;
  int adreport_off;
  int pixel_random;
  int country_random;
  int NID_random;
  int int_gp_1_effect;
  int int_gp_2_effect;
  int main_space_effect;
  int main_time_effect;
  int main_age_effect;
  int main_sex_effect;
  int sz_gp_effect;
  int sx_gp_effect;
  int tx_gp_effect;
  int zx_gp_effect;
  int cre_z_gp_effect;
  int cre_x_gp_effect;
  int use_covs;
  int age_in_int_1_gp;
  int sex_in_int_1_gp;
  int sex_in_int_2_gp;
  int use_anc_corrections;
  int use_error_iid_re;
  
  // Way easier to read these in as vectors of integers and then index the first
  option_list(SEXP x){
    use_priors = asVector<int>(getListElement(x,"use_priors"))[0];
    adreport_off = asVector<int>(getListElement(x,"adreport_off"))[0];
    pixel_random = asVector<int>(getListElement(x,"pixel_random"))[0];
    country_random = asVector<int>(getListElement(x,"country_random"))[0];
    NID_random = asVector<int>(getListElement(x,"NID_random"))[0];
    int_gp_1_effect = asVector<int>(getListElement(x,"int_gp_1_effect"))[0];
    int_gp_2_effect = asVector<int>(getListElement(x,"int_gp_2_effect"))[0];
    main_space_effect = asVector<int>(getListElement(x,"main_space_effect"))[0];
    main_time_effect = asVector<int>(getListElement(x,"main_time_effect"))[0];
    main_age_effect = asVector<int>(getListElement(x,"main_age_effect"))[0];
    main_sex_effect = asVector<int>(getListElement(x,"main_sex_effect"))[0];
    sz_gp_effect = asVector<int>(getListElement(x,"sz_gp_effect"))[0];
    sx_gp_effect = asVector<int>(getListElement(x,"sx_gp_effect"))[0];
    tx_gp_effect = asVector<int>(getListElement(x,"tx_gp_effect"))[0];
    zx_gp_effect = asVector<int>(getListElement(x,"zx_gp_effect"))[0];
    cre_z_gp_effect = asVector<int>(getListElement(x,"cre_z_gp_effect"))[0];
    cre_x_gp_effect = asVector<int>(getListElement(x,"cre_x_gp_effect"))[0];
    age_in_int_1_gp = asVector<int>(getListElement(x,"age_in_int_1_gp"))[0];
    sex_in_int_1_gp = asVector<int>(getListElement(x,"sex_in_int_1_gp"))[0];
    sex_in_int_2_gp = asVector<int>(getListElement(x,"sex_in_int_2_gp"))[0];
    use_covs = asVector<int>(getListElement(x,"use_covs"))[0];
    use_anc_corrections = asVector<int>(getListElement(x,"use_anc_corrections"))[0];
    use_error_iid_re = asVector<int>(getListElement(x,"use_error_iid_re"))[0];
  }
};

// how to read in an object that is used as a prior for standard devs of res
template<class Type>
struct prior_type_sigma {
  std::string name;
  Type par1;
  Type par2;
  
  prior_type_sigma(SEXP x){
    name = CHAR(STRING_ELT(getListElement(x,"type"), 0));
    par1 = asVector<float>(getListElement(x,"par1"))[0];
    par2 = asVector<float>(getListElement(x,"par2"))[0];
  }
};

// how to read in an object that is used as a prior for Matern hyperparameters
template<class Type>
struct prior_type_matern {
  std::string name;
  Type par1a; //mean logtau / rho0
  Type par1b; //prec logtau / alpha_rho
  Type par2a; //mean logkappa / sigma0
  Type par2b; //prec logkappa / alpha_sigma
  
  prior_type_matern(SEXP x){
    name = CHAR(STRING_ELT(getListElement(x,"type"), 0));
    par1a = asVector<float>(getListElement(x,"par1"))[0];
    par1b = asVector<float>(getListElement(x,"par1"))[1];
    par2a = asVector<float>(getListElement(x,"par2"))[0];
    par2b = asVector<float>(getListElement(x,"par2"))[1];
  }
};

// evaluate a prior for sigma using the read in object
template<class Type>
Type eval_prior_sigma(prior_type_sigma<Type> prior, Type log_sigma){
  Type penalty;
  // transform log sigma to log tau space to match INLA prior specification
  // NOTE: log tau --> x=log sigma -->
  //    transform: log tau = -2x --> Jacobian: |J| = 2
  // https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html#sec:priors
  Type tau = pow(exp(log_sigma), -2.);
  Type logtau = log(tau);
  
  if(prior.name == "pc.prec") {
    Type lambda = - log(prior.par2) / prior.par1;
    penalty = -lambda * exp(-logtau/Type(2.0)) - logtau/Type(2.0);
  } 
  else if(prior.name == "normal") {
    // prior.par2 needs to be in precision space such as in INLA
    // https://inla.r-inla-download.org/r-inla.org/doc/prior/gaussian.pdf
    penalty = dnorm(logtau, prior.par1, pow(prior.par2, -.5), true);
  }
  else { //loggamma
    // prior.par2 gamma function for TMB uses shape and scale
    // https://kaskr.github.io/adcomp/group__R__style__distribution.html#gab0e2205710a698ad6a0ed39e0652c9a3
    // INLA uses shape and rate so we need to transforma
    // https://inla.r-inla-download.org/r-inla.org/doc/prior/prior-loggamma.pdf
    penalty = dlgamma(logtau, prior.par1, 1./prior.par2, true);
  }
  
  return penalty;
}

// evaluate priors for matern using the read in object
template<class Type>
Type eval_prior_matern(prior_type_matern<Type> prior, Type logtau, Type logkappa){
  Type penalty;
  
  if(prior.name == "pc") {
    Type d = 2.;
    Type lambda1 = -log(prior.par1b) * pow(prior.par1a, d/2.);
    Type lambda2 = -log(prior.par2b) / prior.par2a;
    Type range   = sqrt(8.0) / exp(logkappa);
    Type sigma   = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
    
    penalty = (-d/2. - 1.) * log(range) - lambda1 * pow(range, -d/2.) - lambda2 * sigma;
    // Note: (rho, sigma) --> (x=log kappa, y=log tau) -->
    //  transforms: rho = sqrt(8)/e^x & sigma = 1/(sqrt(4pi)*e^x*e^y)
    //  --> Jacobian: |J| propto e^(-y -2x)
    Type jacobian = - logtau - 2.0*logkappa;
    penalty += jacobian;
  } 
  else { //normal
    penalty = dnorm(logtau, prior.par1a, 1/sqrt(prior.par1b), true) +
      dnorm(logkappa, prior.par2a, 1/sqrt(prior.par2b), true);
  }
  
  return penalty;
}


// Constrain alpha values
template<class Type>
vector<Type> constrain_pars(vector<Type> alpha, vector<int> constraints){
  int K = alpha.size();
  vector<Type> alpha_c(K);
  
  for(int k = 0; k < K; k++){
    if(constraints[k] == 1){
      alpha_c[k] = exp(alpha[k]);
    }
    if(constraints[k] == -1){
      alpha_c[k] = -1. * exp(alpha[k]);
    }
    if(constraints[k] == 0){
      alpha_c[k] = alpha[k];
    }
  }
  
  return alpha_c;
}

// objective function (ie the likelihood function for the model), returns the evaluated negative log likelihood
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // ////////////////////////////////////////////////////////////////////////////
  // INPUTS
  // ////////////////////////////////////////////////////////////////////////////
  DATA_INTEGER(flag); // flag=0 => only prior
  
  // Data
  DATA_VECTOR(y_i);          // obs successes per binomial experiment at point i (aka cluster)
  DATA_VECTOR(n_i);          // trials per cluster
  
  DATA_INTEGER(num_z);       // number of Z groups
  DATA_INTEGER(num_x);       // number of sexes
  
  DATA_IVECTOR(c_re_i);      // country identifiers
  DATA_IVECTOR(nid_re_i);    // NID identifiers
  DATA_IVECTOR(pixel_re_k);  // pixel identifiers
  DATA_IVECTOR(site_iid_i);  // anc identifiers
  DATA_IVECTOR(error_iid_i);  // country-year identifiers
  DATA_IVECTOR(w_i);         // weights for observations
  DATA_MATRIX(X_kj);          // covariate design matrix
  DATA_IVECTOR(fconstraints);// constraints of fixed effects
  
  // instructions for likelihood asessment
  DATA_VECTOR(lik_gaussian_i); // data likelihood for each row
  DATA_VECTOR(lik_binomial_i); // data likelihood for each row
  DATA_VECTOR(sd_i);           // crossalked standard deviation (set to zero if non-existant)
  DATA_IVECTOR(stacker_col_id); // indicator for whether column is sum-to-1 (1) or not (0)
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0_int);    // used to make gmrf precision--interaction
  DATA_SPARSE_MATRIX(M1_int);    // used to make gmrf precision--interaction
  DATA_SPARSE_MATRIX(M2_int);    // used to make gmrf precision--interaction
  DATA_SPARSE_MATRIX(M0_s);    // used to make gmrf precision--space main effect
  DATA_SPARSE_MATRIX(M1_s);    // used to make gmrf precision--space main effect
  DATA_SPARSE_MATRIX(M2_s);    // used to make gmrf precision--space main effect
  DATA_SPARSE_MATRIX(Aproj_int_1); // used to project from space/time/age mesh to data locations
  DATA_IVECTOR(Aproj_int_2); // used to project from time/age/sex mesh to data locations
  DATA_SPARSE_MATRIX(Aproj_s); // used to project from spatial mesh to data locations
  DATA_IVECTOR(Aproj_t); // used to project from time mesh to data locations
  DATA_IVECTOR(Aproj_z); // used to project from age mesh to data locations
  DATA_IVECTOR(Aproj_x); // used to project from sex mesh to data locations  
  DATA_SPARSE_MATRIX(Aproj_sz); // used to project from sex mesh to data locations  
  DATA_SPARSE_MATRIX(Aproj_sx); // used to project from sex mesh to data locations  
  DATA_IVECTOR(Aproj_tx); // used to project from sex mesh to data locations  
  DATA_IVECTOR(Aproj_zx); // used to project from sex mesh to data locations
  DATA_IVECTOR(Aproj_cre_z); // used to project from sex mesh to data locations    
  DATA_IVECTOR(Aproj_cre_x); // used to project from sex mesh to data locations    
  DATA_IVECTOR(ID_k);        // id values (max value of num_i) for k dispersed aggregation values
  DATA_VECTOR(aggweight_k);  // weight instructions for each dispersed aggregation
  DATA_VECTOR(frr);  // fertility rate ratios for ANC data (non-ANC set to 1)
  DATA_IVECTOR(anc_01_i);  // 0/1 vector specifying type as anc (1) or not (0) for later multiplying by site correction
  
  // Options
  // Boolean vector of options to be used to select different models/modelling 
  //   options: see above
  DATA_STRUCT(options, option_list);  
  
  // Prior specifications
  DATA_STRUCT(prior_log_pixelre_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_cre_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_nidre_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_site_iid_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_error_iid_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_t_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_z_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_x_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_int2_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_cre_z_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_cre_x_sigma, prior_type_sigma);
  
  DATA_STRUCT(prior_matern_int, prior_type_matern);
  DATA_STRUCT(prior_matern_s, prior_type_matern);
  //priors for AR1 effects
  DATA_SCALAR(prior_trho_mean);
  DATA_SCALAR(prior_trho_sd);
  DATA_SCALAR(prior_trho_me_mean);
  DATA_SCALAR(prior_trho_me_sd);
  DATA_SCALAR(prior_zrho_mean);
  DATA_SCALAR(prior_zrho_sd);
  DATA_SCALAR(prior_xrho_mean);
  DATA_SCALAR(prior_xrho_sd);
  DATA_SCALAR(prior_cre_zrho_mean);
  DATA_SCALAR(prior_cre_zrho_sd);
  DATA_SCALAR(prior_cre_xrho_mean);
  DATA_SCALAR(prior_cre_xrho_sd);
  
  // Parameters
  PARAMETER_VECTOR(alpha_j);   // fixed effect coefs, including intercept as first index
  PARAMETER(site_fe);           // fixed effect coefs for site FE
  PARAMETER(logtau);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(trho);             // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho);             // Z autocorrelation parameter for AR1, natural scale
  PARAMETER(xrho);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(logtau_me);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa_me);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(logtau_sz);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa_sz);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(logtau_sx);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa_sx);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(trho_me);             // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho_me);             // Z autocorrelation parameter for AR1, natural scale
  PARAMETER(xrho_me);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(trho_int2);             // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho_int2);             // Z autocorrelation parameter for AR1, natural scale
  PARAMETER(xrho_int2);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(trho_tx);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(xrho_tx);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho_zx);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(xrho_zx);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho_sz);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(xrho_sx);             // Sex autocorrelation parameter for AR1, natural scale
  PARAMETER(log_pixelre_sigma); // log of the standard deviation of the normal error nugget term
  PARAMETER(log_cre_sigma);    // log of the standard deviation of the country random effect
  PARAMETER(log_nidre_sigma);  // log of the standard deviation of the nid random effect
  PARAMETER(log_site_iid_sigma);  // log of the standard deviation of the nid random effect
  PARAMETER(log_error_iid_sigma);  // log of the standard deviation of the nid random effect
  PARAMETER(log_gauss_sigma);  // log of sigma for any gaussian observations
  PARAMETER(log_t_sigma); // log of the standard deviation of the normal error time effect
  PARAMETER(log_z_sigma); // log of the standard deviation of the normal error age effect
  PARAMETER(log_x_sigma); // log of the standard deviation of the normal error sex effect
  PARAMETER(log_int2_sigma); // log of the standard deviation of the normal error 2nd interacting effect
  PARAMETER(log_tx_sigma); // log of the standard deviation of the normal error 2nd interacting effect
  PARAMETER(log_zx_sigma); // log of the standard deviation of the normal error 2nd interacting effect
  PARAMETER(log_cre_z_sigma); // log of the standard deviation of the normal error 2nd interacting effect
  PARAMETER(cre_zrho); // log of the standard deviation of the normal error 2nd interacting effect
  PARAMETER(log_cre_x_sigma); // log of the standard deviation of the normal error 2nd interacting effect
  PARAMETER(cre_xrho); // log of the standard deviation of the normal error 2nd interacting effect
  
  // Random effects
  PARAMETER_ARRAY(Epsilon_int_1);  // Random effects for each int_1 mesh location. Should be 3D array of dimensions num_s by num_t by num_z by num x
  PARAMETER_ARRAY(Epsilon_int_2);  // Random effects for each int_2 mesh location. Should be 3D array of dimensions by num_t by num_z by num_x
  PARAMETER_VECTOR(Epsilon_s);  // Random effects for each S mesh location. Should be 1D array of dimensions num_s
  PARAMETER_VECTOR(Epsilon_t);  // Random effects for each T mesh location. Should be 1D array of dimensions num_t
  PARAMETER_VECTOR(Epsilon_z);  // Random effects for each Z mesh location. Should be 1D array of dimensions num_z
  PARAMETER_VECTOR(Epsilon_x);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_ARRAY(Epsilon_sz);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_ARRAY(Epsilon_sx);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_ARRAY(Epsilon_tx);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_ARRAY(Epsilon_zx);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_MATRIX(Epsilon_cre_z);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_MATRIX(Epsilon_cre_x);  // Random effects for each Sex mesh location. Should be 1D array of dimensions num_sex
  PARAMETER_VECTOR(pixel_re);       // Random effects of the nugget
  PARAMETER_VECTOR(cntry_re);    // Random effects values for country intercept
  PARAMETER_VECTOR(nid_re);      // Random effects values for nid intercept
  PARAMETER_VECTOR(site_iid_re);      // Random effects values for nid intercept
  PARAMETER_VECTOR(error_iid_re);      // Random effects values for country year intercept
  
  printf("Epsilon_int_1 size: %ld \n", Epsilon_int_1.size());
  
  // ////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  // ////////////////////////////////////////////////////////////////////////////
  
  // Define the joint-negative log-likelihood as a parallel_accumulator
  // this allows us to add or subtract numbers to the object in parallel
  // parallel_accumulator<Type> jnll(this);
  Type jnll = 0;
  
  // print parallel info
  // max_parallel_regions = omp_get_max_threads();
  // printf("This is thread %ld\n", max_parallel_regions);
  max_parallel_regions = omp_get_max_threads();
  
  // Make spatial precision matrix--interaction
  SparseMatrix<Type> Q_ss_int   = spde_Q(logkappa, logtau, M0_int, M1_int, M2_int);
  printf("Q_ss_int size: %ld \n", Q_ss_int.size());
  
  
  // Make spatial precision matrix--sz
  SparseMatrix<Type> Q_ss_sz   = spde_Q(logkappa_sz, logtau_sz, M0_int, M1_int, M2_int);
  printf("Q_ss_sz size: %ld \n", Q_ss_sz.size());
  
  
  // Make spatial precision matrix--sx
  SparseMatrix<Type> Q_ss_sx   = spde_Q(logkappa_sx, logtau_sx, M0_int, M1_int, M2_int);
  printf("Q_ss_sx size: %ld \n", Q_ss_sx.size());
  
  
  // Make spatial precision matrix--main effect
  SparseMatrix<Type> Q_ss_s   = spde_Q(logkappa_me, logtau_me, M0_s, M1_s, M2_s);
  printf("Q_ss_s size: %ld \n", Q_ss_s.size());
  
  int num_i = y_i.size();           // Number of observations
  int num_j = X_kj.cols();          // Number of covariates
  int num_s = Epsilon_int_1.dim(0);   // Number of mesh nodes
  int num_t = Epsilon_int_1.dim(1);   // Number of time periods
  int num_k = ID_k.size();          // Number of disaggregated points
  int num_cre = cntry_re.size();
  
  // Make transformations of some of our parameters
  Type range         = sqrt(8.0) / exp(logkappa);
  Type spde_sigma    = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
  Type trho_trans    = (exp(trho) - 1) / (exp(trho) + 1);
  Type zrho_trans    = (exp(zrho) - 1) / (exp(zrho) + 1);
  Type xrho_trans    = (exp(xrho) - 1) / (exp(xrho) + 1);
  Type spde_me_sigma    = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau_me) * exp(2.0 * logkappa_me));
  Type spde_sz_sigma    = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau_sz) * exp(2.0 * logkappa_sz));
  Type spde_sx_sigma    = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau_sx) * exp(2.0 * logkappa_sx));
  Type trho_me_trans    = (exp(trho_me) - 1) / (exp(trho_me) + 1);
  Type zrho_me_trans    = (exp(zrho_me) - 1) / (exp(zrho_me) + 1);
  Type xrho_me_trans    = (exp(xrho_me) - 1) / (exp(xrho_me) + 1);
  Type trho_int2_trans    = (exp(trho_int2) - 1) / (exp(trho_int2) + 1);
  Type zrho_int2_trans    = (exp(zrho_int2) - 1) / (exp(zrho_int2) + 1);
  Type xrho_int2_trans    = (exp(xrho_int2) - 1) / (exp(xrho_int2) + 1);
  
  Type trho_tx_trans    = (exp(trho_tx) - 1) / (exp(trho_tx) + 1);
  Type xrho_tx_trans    = (exp(xrho_tx) - 1) / (exp(xrho_tx) + 1);
  Type zrho_zx_trans    = (exp(zrho_zx) - 1) / (exp(zrho_zx) + 1);
  Type xrho_zx_trans    = (exp(xrho_zx) - 1) / (exp(xrho_zx) + 1);
  Type zrho_sz_trans    = (exp(zrho_sz) - 1) / (exp(zrho_sz) + 1);
  Type xrho_sx_trans    = (exp(xrho_sx) - 1) / (exp(xrho_sx) + 1);
  Type pixelre_sigma = exp(log_pixelre_sigma);
  Type cre_sigma     = exp(log_cre_sigma);
  Type site_iid_sigma     = exp(log_site_iid_sigma);
  Type error_iid_sigma     = exp(log_error_iid_sigma);
  Type nidre_sigma   = exp(log_nidre_sigma);
  Type gauss_sigma   = exp(log_gauss_sigma);
  Type t_sigma       = exp(log_t_sigma);
  Type z_sigma       = exp(log_z_sigma);
  Type x_sigma       = exp(log_x_sigma);
  Type int2_sigma       = exp(log_int2_sigma);
  Type tx_sigma       = exp(log_tx_sigma);
  Type zx_sigma       = exp(log_zx_sigma);
  
   Type cre_z_sigma       = exp(log_cre_z_sigma);
   Type cre_zrho_trans    = (exp(cre_zrho) - 1) / (exp(cre_zrho) + 1); //TRANSOFRM from -inf, inf to -1, 1.. // log((1.1 + zrho) / (1.1 - zrho));

   Type cre_x_sigma       = exp(log_cre_x_sigma);
   Type cre_xrho_trans    = (exp(cre_xrho) - 1) / (exp(cre_xrho) + 1); //TRANSOFRM from -inf, inf to -1, 1.. // log((1.1 + zrho) / (1.1 - zrho));
  
  // Define objects for derived values
  vector<Type> fe_k(num_k);                         // fixed effects portion of linear predictor
  vector<Type> eta_k(num_k);                        // estimated linear predictor for each point k
  vector<Type> prob_i(num_i);                       // probability for each data observation
  vector<Type> prob_i_frr_corrected(num_i);                       // probability for each data observation (after anc frr corrections, or else same as prob_i)
  vector<Type> prob_i_pre_error(num_i);             // a version of prob i to add the observation error term to
  vector<Type> prob_k(num_k);                       // probability for each disaggregated data observation
  vector<Type> prob_anc_k(num_k);                       // probability for each disaggregated data observation converted to anc probability
  prob_i.setZero();
  vector<Type> projepsilon_int_1_k(num_k);                // value of gmrf 1 at data points
  vector<Type> projepsilon_sz_k(num_k);                // value of gmrf sz at data points
  vector<Type> projepsilon_sx_k(num_k);                // value of gmrf sx at data points
  vector<Type> projepsilon_s_k(num_k);                // value of gmrf at data points
  vector<Type> epsilon_int_1(num_s * num_t);         // Epsilon_int_1 unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_int_2(num_t * num_z);         // Epsilon_int_2 unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_sz(num_s * num_z);         // Epsilon_sz unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_sx(num_s * num_x);         // Epsilon_sx unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_tx(num_t * num_x);         // Epsilon_tx unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_zx(num_z * num_x);         // Epsilon_zx unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_cre_z(num_cre * num_z);         // Epsilon_cre_z unlisted into a vector for easier matrix multiplication
  vector<Type> epsilon_cre_x(num_cre * num_x);         // Epsilon_cre_x unlisted into a vector for easier matrix multiplication
  
  //Resize epsilon_int_1 depending on structure of interacting GP
  if(options.age_in_int_1_gp == 1) {
    epsilon_int_1.resize(num_s * num_t * num_z);
  }
  
  if(options.sex_in_int_1_gp == 1) {
    epsilon_int_1.resize(num_s * num_t * num_z * num_x);
  }
  
  //Resize epsilon_int_2 depending on structure of interacting GP
  if(options.sex_in_int_2_gp == 1) {
    epsilon_int_2.resize(num_t * num_z * num_x);
  }
  
  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: S, ST, SZ, STZ, and STZX
  if (options.int_gp_1_effect == 1) {
    if (num_t == 1 & num_z == 1)  {
      printf("GP FOR SPACE  ONLY \n");
      PARALLEL_REGION jnll += GMRF(Q_ss_int,false)(epsilon_int_1);
    } else if(num_t > 1 & (num_z == 1 | options.age_in_int_1_gp == 0)) {
      printf("GP FOR SPACE-TIME \n");
      PARALLEL_REGION jnll += SEPARABLE(AR1(trho_trans),GMRF(Q_ss_int,false))(Epsilon_int_1);
    } else if (num_t == 1 & num_z > 1) {
      printf("GP FOR SPACE-Z \n");
      PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_trans),GMRF(Q_ss_int,false))(Epsilon_int_1);
    } else if (num_t > 1 & (num_z > 1 & options.age_in_int_1_gp == 1) & (num_x == 1 | options.sex_in_int_1_gp == 0)) {
      printf("GP FOR SPACE-TIME-Z \n");
      PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_trans),SEPARABLE(AR1(trho_trans),GMRF(Q_ss_int,false)))(Epsilon_int_1);
    } else if (num_t > 1 & (num_z > 1 & options.age_in_int_1_gp == 1) & (num_x > 1 & options.sex_in_int_1_gp == 1)) {
      printf("GP FOR SPACE-TIME-Z-SEX \n");
      PARALLEL_REGION jnll += SEPARABLE(AR1(xrho_trans),SEPARABLE(AR1(zrho_trans),SEPARABLE(AR1(trho_trans),GMRF(Q_ss_int,false))))(Epsilon_int_1);
    } 
  }
  
  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: SZ
  if (options.sz_gp_effect == 1) {
    printf("GP FOR SPACE-Z \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_sz_trans),GMRF(Q_ss_sz,false))(Epsilon_sz);
  }
  
  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: SZ
  if (options.sx_gp_effect == 1) {
    printf("GP FOR SPACE-X \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(xrho_sx_trans),GMRF(Q_ss_sx,false))(Epsilon_sx);
  }
  
  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: TZ, and TZX
  if (options.int_gp_2_effect == 1) {
    if (num_x == 1 | options.sex_in_int_2_gp == 0) {
      printf("GP FOR TIME-Z \n");
      PARALLEL_REGION jnll += SCALE(SEPARABLE(AR1(zrho_int2_trans),AR1(trho_int2_trans)), int2_sigma)(Epsilon_int_2);
    } else if (num_x > 1 & options.sex_in_int_2_gp == 1) {
      printf("GP FOR TIME-Z-SEX \n");
      PARALLEL_REGION jnll += SCALE(SEPARABLE(AR1(xrho_int2_trans),SEPARABLE(AR1(zrho_int2_trans),AR1(trho_int2_trans))), int2_sigma)(Epsilon_int_2);
    } 
  }
  
  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: TZ, and TZX
  if (options.tx_gp_effect == 1) {
    printf("GP FOR T-X \n");
    PARALLEL_REGION jnll += SCALE(SEPARABLE(AR1(trho_tx_trans),AR1(xrho_tx_trans)), tx_sigma)(Epsilon_tx);
  }
  
  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: TZ, and TZX
  if (options.zx_gp_effect == 1) {
    printf("GP FOR Z-X \n");
    PARALLEL_REGION jnll += SCALE(SEPARABLE(AR1(zrho_zx_trans),AR1(xrho_zx_trans)), zx_sigma)(Epsilon_zx);
  }
  
  // pixel random intercept
  if(options.pixel_random == 1 ){
    printf("Pixel RE \n");
    for(int i=0; i<pixel_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(pixel_re(i), Type(0.0), pixelre_sigma, true);
    }
  }
  
  // country random intercept
  if(options.country_random == 1 ){
    printf("Country RE \n");
    for(int i=0; i<cntry_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(cntry_re(i), Type(0.0), cre_sigma, true);
    }
  }
  
  
  // country specific z intercept
  if(options.cre_z_gp_effect == 1 ){
    printf("Country-specific Z \n");
    for(int i=0; i<num_cre; i++){
        PARALLEL_REGION jnll += SCALE(AR1(cre_zrho_trans), cre_z_sigma)(Epsilon_cre_z.row(i));
    }
  }
  
  if(options.cre_x_gp_effect == 1 ){
    printf("Country-specific X \n");
    for(int i=0; i<num_cre; i++){
      PARALLEL_REGION jnll += SCALE(AR1(cre_xrho_trans), cre_x_sigma)(Epsilon_cre_x.row(i));
    }
  }
  
  // nid random intercept
  if(options.NID_random == 1 ){
    printf("NID RE \n");
    for(int i=0; i<nid_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(nid_re(i), Type(0.0), nidre_sigma, true);
    }
  }
  
  // Latent field/Random effect contribution to likelihood.
  if (options.main_space_effect == 1)  {
    printf("GP FOR SPACE  ONLY \n");
    PARALLEL_REGION jnll += GMRF(Q_ss_s,false)(Epsilon_s);
  }   
  if (options.main_time_effect == 1)  {
    printf("GP FOR TIME  ONLY \n");
    PARALLEL_REGION jnll += SCALE(AR1(trho_me_trans), t_sigma)(Epsilon_t);
  }  
  if (options.main_age_effect == 1)  {
    printf("GP FOR AGE  ONLY \n");
    PARALLEL_REGION jnll += SCALE(AR1(zrho_me_trans), z_sigma)(Epsilon_z);
  }
  if (options.main_sex_effect == 1)  {
    printf("GP FOR SEX  ONLY \n");
    PARALLEL_REGION jnll += SCALE(AR1(xrho_me_trans), x_sigma)(Epsilon_x);
  }  
  if (options.use_anc_corrections == 1)  {
    //printf("GP FOR ANC SITE ONLY \n");
    printf("IID RE FOR ANC SITE \n");
    for(int i=0; i<site_iid_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(site_iid_re(i), Type(0.0), site_iid_sigma, true);
    }
  } 
  if (options.use_error_iid_re == 1)  {
    //printf("GP FOR ANC SITE ONLY \n");
    printf("IID RE FOR COUNTRY-YEAR \n");
    for(int i=0; i<error_iid_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(error_iid_re(i), Type(0.0), error_iid_sigma, true);
    }
  } 
  //####################################################################################
  
  // Transform int_1 GMRFs and make vector form
  if (options.int_gp_1_effect == 1)  {
    printf("Transform int_1 GMRF \n");
    for(int s = 0; s < num_s; s++){
      for(int t = 0; t < num_t; t++){
        if(num_z == 1 | options.age_in_int_1_gp == 0) {
          epsilon_int_1[(s + num_s * t )] = Epsilon_int_1(s,t);
        } else {
          if(num_x == 1 | options.sex_in_int_1_gp == 0){
            for(int z = 0; z < num_z; z++){
              epsilon_int_1[(s + num_s * t + num_s * num_t * z)] = Epsilon_int_1(s,t,z);
            }
          } else {
            for(int z = 0; z < num_z; z++){
              for(int x = 0; x < num_x; x++){
                epsilon_int_1[(s + num_s * t + num_s * num_t * z + num_s * num_t * num_z * x)] = Epsilon_int_1(s,t,z,x);
              }
            }
          }
        }
      }
    }
  }
  
  // Transform sz GMRFs and make vector form
  if (options.sz_gp_effect == 1)  {
    printf("Transform SZ GMRF \n");
    for(int s = 0; s < num_s; s++){
      for(int z = 0; z < num_z; z++){
        epsilon_sz[(s + num_s * z )] = Epsilon_sz(s,z);
      }
    }
  }
  
  // Transform sx GMRFs and make vector form
  if (options.sx_gp_effect == 1)  {
    printf("Transform SX GMRF \n");
    for(int s = 0; s < num_s; s++){
      for(int x = 0; x < num_x; x++){
        epsilon_sx[(s + num_s * x )] = Epsilon_sx(s,x);
      }
    }
  }
  
  // Transform int_2 GMRFs and make vector form
  if (options.int_gp_2_effect == 1)  {
    printf("Transform int_2 GMRF \n");
    for(int t = 0; t < num_t; t++){
      if(num_x == 1 | options.sex_in_int_2_gp == 0){
        for(int z = 0; z < num_z; z++){
          epsilon_int_2[(t + num_t * z)] = Epsilon_int_2(t,z);
        }
      } else {
        for(int z = 0; z < num_z; z++){
          for(int x = 0; x < num_x; x++){
            epsilon_int_2[(t + num_t * z + num_t * num_z * x)] = Epsilon_int_2(t,z,x);
          }
        }
      }
    }
  }
  
  // Transform tx GMRFs and make vector form
  if (options.tx_gp_effect == 1)  {
    printf("Transform tx GMRF \n");
    for(int t = 0; t < num_t; t++){
      for(int x = 0; x < num_x; x++){
        epsilon_tx[(t + num_t * x )] = Epsilon_tx(t,x);
      }
    }
  }
  
  // Transform zx GMRFs and make vector form
  if (options.zx_gp_effect == 1)  {
    printf("Transform zx GMRF \n");
    for(int z = 0; z < num_z; z++){
      for(int x = 0; x < num_x; x++){
        epsilon_zx[(z + num_z * x )] = Epsilon_zx(z,x);
      }
    }
  }  
  
  //#######################################################
  // ***********************************
  // Transform cre_z GMRFs and make vector form
  if (options.cre_z_gp_effect == 1)  {
    printf("Transform cre z GMRF \n");
    for(int c = 0; c < num_cre; c++){
      for(int z = 0; z < num_z; z++){
              epsilon_cre_z[(c + num_cre * z )] = Epsilon_cre_z(c,z);
      }
    }
  }
  
  if (options.cre_x_gp_effect == 1)  {
    printf("Transform cre x GMRF \n");
    for(int c = 0; c < num_cre; c++){
      for(int x = 0; x < num_x; x++){
        epsilon_cre_x[(c + num_cre * x )] = Epsilon_cre_x(c,x);
      }
    }
  }
  //**********************************************************
  //#############################################################
  
  // Project from mesh points to data points in order to eval likelihood at each data point, for interaction and main effects
  if (options.int_gp_1_effect == 1)  {
    printf("Project Epsilon int_1 \n");
    projepsilon_int_1_k = Aproj_int_1 * epsilon_int_1.matrix();
  }
  
  if (options.main_space_effect == 1)  {
    printf("Project Epsilon s \n");
    projepsilon_s_k = Aproj_s * Epsilon_s;
  }
  
  if (options.sz_gp_effect == 1)  {
    printf("Project Epsilon sz \n");
    projepsilon_sz_k = Aproj_sz * epsilon_sz.matrix();
  }
  
  if (options.sx_gp_effect == 1)  {
    printf("Project Epsilon sx \n");
    projepsilon_sx_k = Aproj_sx * epsilon_sx.matrix();
  }
  
  // Return un-normalized density on request
  if (flag == 0) return jnll;
  
  // Prior contribution to likelihood. Values are defaulted (for now). Only run if options[0]==1
  if(options.use_priors == 1) {
    if (options.int_gp_1_effect == 1 ) {
      PARALLEL_REGION jnll -= eval_prior_matern(prior_matern_int, logtau, logkappa);    
      PARALLEL_REGION jnll -= dnorm(trho,  Type(prior_trho_mean), Type(prior_trho_sd), true);   //1.2^2   
      if(options.age_in_int_1_gp==1) {
        PARALLEL_REGION jnll -= dnorm(zrho,  Type(prior_zrho_mean), Type(prior_zrho_sd), true);
      }
      if(options.sex_in_int_1_gp==1) {
        PARALLEL_REGION jnll -= dnorm(xrho,  Type(prior_xrho_mean), Type(prior_xrho_sd), true);
      }
    }
    
    if(options.int_gp_2_effect == 1) {
      PARALLEL_REGION jnll -= dnorm(trho_int2,  Type(prior_trho_mean), Type(prior_trho_sd), true);   //1.2^2
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_int2_sigma, log_int2_sigma);
      PARALLEL_REGION jnll -= dnorm(zrho_int2,  Type(prior_zrho_mean), Type(prior_zrho_sd), true);    
      if(options.sex_in_int_2_gp==1) {
        PARALLEL_REGION jnll -= dnorm(xrho_int2,  Type(prior_xrho_mean), Type(prior_xrho_sd), true);
      }  
    }
    
    if(options.sz_gp_effect == 1) {
      PARALLEL_REGION jnll -= dnorm(zrho_sz,  Type(prior_zrho_mean), Type(prior_zrho_sd), true);   //1.2^2
      PARALLEL_REGION jnll -= eval_prior_matern(prior_matern_int, logtau_sz, logkappa_sz); 
    } 
    
    if(options.sx_gp_effect == 1) {
      PARALLEL_REGION jnll -= dnorm(xrho_sx,  Type(prior_xrho_mean), Type(prior_xrho_sd), true);   //1.2^2
      PARALLEL_REGION jnll -= eval_prior_matern(prior_matern_int, logtau_sx, logkappa_sx); 
    } 
    
    if(options.tx_gp_effect == 1) {
      PARALLEL_REGION jnll -= dnorm(trho_tx,  Type(prior_trho_mean), Type(prior_trho_sd), true);   //1.2^2
      PARALLEL_REGION jnll -= dnorm(xrho_tx,  Type(prior_xrho_mean), Type(prior_xrho_sd), true); 
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_int2_sigma, log_tx_sigma);
    }  
    
    if(options.zx_gp_effect == 1) {
      PARALLEL_REGION jnll -= dnorm(zrho_zx,  Type(prior_zrho_mean), Type(prior_zrho_sd), true);   //1.2^2
      PARALLEL_REGION jnll -= dnorm(xrho_zx,  Type(prior_xrho_mean), Type(prior_xrho_sd), true); 
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_int2_sigma, log_zx_sigma);
    } 
    
    if (options.main_space_effect == 1) {
      PARALLEL_REGION jnll -= eval_prior_matern(prior_matern_s, logtau_me, logkappa_me);
    }
    
    if (options.main_time_effect == 1)  {
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_t_sigma, log_t_sigma);
      PARALLEL_REGION jnll -= dnorm(trho_me,  Type(prior_trho_me_mean), Type(prior_trho_me_sd), true);   //1.2^2
    }
    
    
    if (options.main_age_effect == 1)  {
      PARALLEL_REGION jnll -= dnorm(zrho_me,  Type(prior_zrho_mean), Type(prior_zrho_sd), true);
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_z_sigma, log_z_sigma);
    }
    
    if (options.main_sex_effect == 1)  {
      PARALLEL_REGION jnll -= dnorm(xrho_me,  Type(prior_xrho_mean), Type(prior_xrho_sd), true);
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_x_sigma, log_x_sigma);
    }
    
    if (options.cre_z_gp_effect == 1)  {
          PARALLEL_REGION jnll -= dnorm(cre_zrho,  Type(prior_cre_zrho_mean), Type(prior_cre_zrho_sd), true);
          PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_cre_z_sigma, log_cre_z_sigma);   
    }
    
    if (options.cre_x_gp_effect == 1)  {
        PARALLEL_REGION jnll -= dnorm(cre_xrho,  Type(prior_cre_xrho_mean), Type(prior_cre_xrho_sd), true);
        PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_cre_x_sigma, log_cre_x_sigma);  
    }
    
    
    if (options.use_anc_corrections == 1) {
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_site_iid_sigma, log_site_iid_sigma);
      PARALLEL_REGION jnll -= dnorm(site_fe, Type(0.0), Type(3.0), true);
    }
    
    if (options.use_error_iid_re == 1) {
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_error_iid_sigma, log_error_iid_sigma);
    }
    
    for(int j = 0; j < alpha_j.size(); j++){
      if(stacker_col_id(j) == 0) {
        // N(0,3) prior for fixed effects.
        PARALLEL_REGION jnll -= dnorm(alpha_j(j), Type(0.0), Type(3.0), true);
      } else {
        // Dirichlet(1,...,1) prior for stackers
        // done by using gamma and then transforming
        PARALLEL_REGION jnll -= dlgamma(alpha_j(j), Type(1.0), Type(1.0), true);
      }
    }
    
    // if using nugget (option in 3rd index)
    if(options.pixel_random == 1){
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_pixelre_sigma, log_pixelre_sigma);
    }
    // if using country (option in 4th index)
    if(options.country_random == 1 ){
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_cre_sigma, log_cre_sigma);
    }
    // if using nid (option in 5th index)
    if(options.NID_random == 1 ){
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_nidre_sigma, log_nidre_sigma);
    }
    
    
    
  }
  
  // evaluate fixed effects for beta_j values
  vector<Type> beta_j = alpha_j;
  Type norm_factor = 0;
  for(int j = 0; j < alpha_j.size(); j++) {
    if(stacker_col_id(j) == 1) {
      norm_factor += exp(alpha_j(j));
    }
  }
  for(int j = 0; j < alpha_j.size(); j++) {
    if(stacker_col_id(j) == 1) {
      beta_j(j) = exp(alpha_j(j))/norm_factor;
    }
  }
  vector<Type> cbeta_j = constrain_pars(beta_j, fconstraints);
  fe_k = X_kj * cbeta_j.matrix();
  
  // Likelihood contribution from each datapoint i
  // involves iterating over IDs to allow for aggregation over space, age, etc
  printf("Data likelihood \n");
  int i = 0;
  
  for(int k = 0; k < num_k; k++) {
    i = ID_k(k);
    // eta_k(k) = projepsilon_k(k);
    eta_k(k) = fe_k(k) + cntry_re(c_re_i(i)) + nid_re(nid_re_i(i)) + pixel_re(pixel_re_k(k));
    
    if (options.int_gp_1_effect == 1)  {
      eta_k(k) += projepsilon_int_1_k(k);
    }
    
    if (options.sz_gp_effect == 1)  {
      eta_k(k) += projepsilon_sz_k(k);
    }
    
    if (options.int_gp_2_effect == 1)  {
      eta_k(k) += epsilon_int_2(Aproj_int_2(k));
    }
    
    if (options.tx_gp_effect == 1)  {
      eta_k(k) += epsilon_tx(Aproj_tx(k));
    }
    
    if (options.zx_gp_effect == 1)  {
      eta_k(k) += epsilon_zx(Aproj_zx(k));
    }
    
    if(options.main_space_effect==1){
      eta_k(k) += projepsilon_s_k(k);
    }
    
    if(options.main_time_effect==1){
      eta_k(k) += Epsilon_t(Aproj_t(k));
    }
    if(options.main_age_effect==1){
      eta_k(k) += Epsilon_z(Aproj_z(k));
    }
    
    if(options.main_sex_effect==1){
      eta_k(k) += Epsilon_x(Aproj_x(k));
    }
    
    if(options.cre_z_gp_effect==1){
          eta_k(k) += epsilon_cre_z(Aproj_cre_z(k));
    }
    
    if(options.cre_x_gp_effect==1){
      eta_k(k) += epsilon_cre_x(Aproj_cre_x(k));
    }
    
    
    if(lik_binomial_i(i) == 1){
      prob_k(k) = invlogit_robust(eta_k(k));
      
      
      if(options.use_anc_corrections==1){
        //Convert HIV prevalence to ANC prevalence. when frr(i) == 1 (i.e. not ANC data), these will be equivalent.
        prob_anc_k(k) = (prob_k(k) * frr(k)) / (prob_k(k) * frr(k) + 1 - prob_k(k));
        
        prob_i_frr_corrected(i) += prob_anc_k(k) * aggweight_k(k);
        
        prob_i_pre_error(i) = invlogit_robust(logit(prob_i_frr_corrected(i)) + (site_iid_re(site_iid_i(i))+site_fe) * anc_01_i(i));
      } else{
        prob_i_pre_error(i) += prob_k(k) * aggweight_k(k);
      }
      
      if(options.use_error_iid_re==1){
        prob_i(i) = invlogit_robust(logit(prob_i_pre_error(i)) + error_iid_re(error_iid_i(i)));
      } else {
        prob_i(i) = prob_i_pre_error(i);
      }
      
    }

    
    if(lik_gaussian_i(i) == 1){
      prob_i(i) += eta_k(k) * aggweight_k(k);
    }
    
    // Once we have aggregated all relevant probabilities, evaluate the density
    // occurs if the next ID corresponds to a different observation or at the end of ID_k
    if(k != num_k - 1) {
      if(i != ID_k(k+1)) {
        if(lik_binomial_i(i) == 1){
          PARALLEL_REGION jnll -= dbinom(y_i(i), n_i(i), prob_i(i), true) * w_i(i);
        }
        if(lik_gaussian_i(i) == 1){
          // this includes any crosswalked sd (sd_i), other variance (gauss_sigma), and is scaled by sample size (n_i)
          PARALLEL_REGION jnll -= dnorm( y_i(i), prob_i(i),  sqrt( (pow(sd_i(i),2) + (1/n_i(i) * pow(gauss_sigma,2)) ) ), true ) * w_i(i);
        }
      }
    } else { // at the end of all IDS/observations
      if(lik_binomial_i(i) == 1){
        PARALLEL_REGION jnll -= dbinom(y_i(i), n_i(i), prob_i(i), true) * w_i(i);
      }
      if(lik_gaussian_i(i) == 1){
        // this includes any crosswalked sd (sd_i), other variance (gauss_sigma), and is scaled by sample size (n_i)
        PARALLEL_REGION jnll -= dnorm( y_i(i), prob_i(i),  sqrt( (pow(sd_i(i),2) + (1/n_i(i) * pow(gauss_sigma,2)) ) ), true ) * w_i(i);
      }
    }
    
  }
  
  // Report estimates
  if(options.adreport_off == 0){
    ADREPORT(alpha_j);
    ADREPORT(beta_j);
  }
  
  return jnll;
}