// ///////////////////////////////////////////////////
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
  int nugget;
  int country_random;
  int NID_random;
  int useGP;
  // Way easier to read these in as vectors of integers and then index the first
  option_list(SEXP x){
    use_priors = asVector<int>(getListElement(x,"use_priors"))[0];
    adreport_off = asVector<int>(getListElement(x,"adreport_off"))[0];
    nugget = asVector<int>(getListElement(x,"nugget"))[0];
    country_random = asVector<int>(getListElement(x,"country_random"))[0];
    NID_random = asVector<int>(getListElement(x,"NID_random"))[0];
    useGP = asVector<int>(getListElement(x,"useGP"))[0];
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

  // Indices
  DATA_INTEGER(num_i);       // number of datapts in space-time-Z (aka STZ)
  DATA_INTEGER(num_s);       // number of mesh pts in space mesh
  DATA_INTEGER(num_t);       // number of time periods
  DATA_INTEGER(num_z);       // number of Z groups

  // Data (each, excpect for X_ij is a vector of length num_i)
  DATA_VECTOR(y_i);          // obs successes per binomial experiment at point i (aka cluster)
  DATA_VECTOR(n_i);          // trials per cluster
  DATA_IVECTOR(t_i);         // time period of the data point
  DATA_IVECTOR(c_re_i);      // country identifiers
  DATA_IVECTOR(nid_re_i);    // NID identifiers
  DATA_IVECTOR(w_i);         // weights for observations
  DATA_MATRIX(X_ij);         // covariate design matrix (num_i by number of fixed effects matrix)
  DATA_IVECTOR(fconstraints); // constraints of fixed effects

  // instructions for likelihood asessment
  DATA_VECTOR(lik_gaussian_i); // data likelihood for each row
  DATA_VECTOR(lik_binomial_i); // data likelihood for each row
  DATA_VECTOR(sd_i);           // crossalked standard deviation (set to zero if non-existant)

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations

  // Boolean vector of options to be used to select different models/modelling
  //   options: see above
  DATA_STRUCT(options, option_list);

  // Prior specifications
  DATA_STRUCT(prior_log_nugget_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_cre_sigma, prior_type_sigma);
  DATA_STRUCT(prior_log_nidre_sigma, prior_type_sigma);
  DATA_STRUCT(prior_matern, prior_type_matern);

  // Parameters
  PARAMETER_VECTOR(alpha_j);   // fixed effect coefs, including intercept as first index
  PARAMETER(logtau);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(trho);             // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho);             // Z autocorrelation parameter for AR1, natural scale
  PARAMETER(log_nugget_sigma); // log of the standard deviation of the normal error nugget term
  PARAMETER(log_cre_sigma);    // log of the standard deviation of the country random effect (later do as vec if using Random SLOPE TODO)
  PARAMETER(log_nidre_sigma);  // log of the standard deviation of the nid random effect
  PARAMETER(log_gauss_sigma);  // log of sigma for any gaussian observations

  // Random effects
  PARAMETER_ARRAY(Epsilon_stz);  // Random effects for each STZ mesh location. Should be 3D array of dimensions num_s by num_t by num_z
  PARAMETER_VECTOR(nug_i);       // Random effects of the nugget
  PARAMETER_VECTOR(cntry_re);    // Random effects values for country intercept
  PARAMETER_VECTOR(nid_re);      // Random effects values for nid intercept

  printf("Epsilon_stz size: %ld \n", Epsilon_stz.size());

  // ////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  // ////////////////////////////////////////////////////////////////////////////

  // Define the joint-negative log-likelihood as a parallel_accumulator
  // this allows us to add or subtract numbers to the object in parallel
  // parallel_accumulator<Type> jnll(this);
  Type jnll = 0;

  // print parallel info
  max_parallel_regions = omp_get_max_threads();

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss   = spde_Q(logkappa, logtau, M0, M1, M2);
  printf("Q_ss size: %ld \n", Q_ss.size());

  // Make transformations of some of our parameters
  Type range         = sqrt(8.0) / exp(logkappa);
  Type sigma         = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
  Type trho_trans    = (exp(trho) - 1) / (exp(trho) + 1); // TRANSOFRM from -inf, inf to -1, 1.. //log((1.1 + trho) / (1.1 - trho));
  Type zrho_trans    = (exp(zrho) - 1) / (exp(zrho) + 1); //TRANSOFRM from -inf, inf to -1, 1.. // log((1.1 + zrho) / (1.1 - zrho));
  Type nugget_sigma  = exp(log_nugget_sigma);
  Type cre_sigma     = exp(log_cre_sigma);
  Type nidre_sigma   = exp(log_nidre_sigma);
  Type gauss_sigma   = exp(log_gauss_sigma);


  // Define objects for derived values
  vector<Type> fe_i(num_i);                         // main effect X_ij %*% t(alpha_j)
  vector<Type> epsilon_stz(num_s * num_t * num_z);  // Epsilon_stz unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);                // value of gmrf at data points
  vector<Type> prob_i(num_i);                       // Logit estimated prob for each point i

  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: S, ST, SZ, and STZ
  if (num_t == 1 & num_z == 1)  {
    printf("GP FOR SPACE  ONLY \n");
    PARALLEL_REGION jnll += GMRF(Q_ss,false)(epsilon_stz);
  } else if(num_t > 1 & num_z == 1) {
    printf("GP FOR SPACE-TIME \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(trho_trans),GMRF(Q_ss,false))(Epsilon_stz);
  } else if (num_t == 1 & num_z > 1) {
    printf("GP FOR SPACE-Z \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_trans),GMRF(Q_ss,false))(Epsilon_stz);
  } else if (num_t > 1 & num_z > 1) {
    printf("GP FOR SPACE-TIME-Z \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_trans),SEPARABLE(AR1(trho_trans),GMRF(Q_ss,false)))(Epsilon_stz);
  }

  // nugget contribution to the likelihood
  if(options.nugget == 1){
    printf("Nugget \n");
    for (int i = 0; i < num_i; i++){
      // binomial models with sd_i, the additional variance gets put into the nugget.
      // for gaussian models nuggets seem unidentifiable so, the sd_i is in the data likelihood
      PARALLEL_REGION jnll -= dnorm(nug_i(i), Type(0.0), sqrt( pow(sd_i(i),2) + pow(nugget_sigma,2) ), true);
    }
  }

  // country random intercept
  if(options.country_random == 1){
    printf("Country RE \n");
    for(int i=0; i<cntry_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(cntry_re(i), Type(0.0), cre_sigma, true);
    }
  }

  // nid random intercept
  if(options.NID_random == 1){
    printf("NID RE \n");
    for(int i=0; i<nid_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(nid_re(i), Type(0.0), nidre_sigma, true);
    }
  }

  // Transform GMRFs and make vector form
  printf("Transform GMRF \n");
  for(int s = 0; s < num_s; s++){
    for(int t = 0; t < num_t; t++){
      if(num_z == 1) {
        epsilon_stz[(s + num_s * t )] = Epsilon_stz(s,t);
      } else {
        for(int z = 0; z < num_z; z++){
          epsilon_stz[(s + num_s * t + num_s * num_t * z)] = Epsilon_stz(s,t,z);
        }
      }
    }
  }

  // Project from mesh points to data points in order to eval likelihood at each data point
  printf("Project Epsilon \n");
  projepsilon_i = Aproj * epsilon_stz.matrix();

  // evaluate fixed effects for alpha_j values
  vector<Type> calpha_j = constrain_pars(alpha_j, fconstraints);
  fe_i = X_ij * calpha_j.matrix();


  // Return un-normalized density on request
  if (flag == 0) return jnll;

  // Prior contribution to likelihood. Values are defaulted.
  // Only run if options.use_priors==1
  if(options.use_priors == 1) {
    PARALLEL_REGION jnll -= eval_prior_matern(prior_matern, logtau, logkappa);
    if(num_t > 1) {
      // N(0,2.58^2) prior on log((1+rho)/(1-rho))
      PARALLEL_REGION jnll -= dnorm(trho, Type(0.0), Type(2.58), true);
    }
    if(num_z > 1) {
      // N(0,2.58^2) prior on log((1+rho)/(1-rho))
      PARALLEL_REGION jnll -= dnorm(zrho, Type(0.0), Type(2.58), true);
    }
    for(int j = 0; j < alpha_j.size(); j++){
      // N(0,3) prior for fixed effects.
      PARALLEL_REGION jnll -= dnorm(alpha_j(j), Type(0.0), Type(3.0), true);
    }
    // if using nugget (option in 3rd index)
    if(options.nugget == 1){
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_nugget_sigma, log_nugget_sigma);
    }
    // if using country (option in 4th index)
    if(options.country_random == 1){
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_cre_sigma, log_cre_sigma);
    }
    // if using nid (option in 5th index)
    if(options.NID_random == 1){
      PARALLEL_REGION jnll -= eval_prior_sigma(prior_log_nidre_sigma, log_nidre_sigma);
    }
  }

  // Likelihood contribution from each datapoint i
  printf("Data likelihood \n");
  for (int i = 0; i < num_i; i++){

    // mean model
    if(options.useGP==1){
      prob_i(i) = fe_i(i) + projepsilon_i(i) + nug_i(i) + cntry_re(c_re_i(i)) + nid_re(nid_re_i(i));
    } else {
      prob_i(i) = fe_i(i) +  nug_i(i) + cntry_re(c_re_i(i)) + nid_re(nid_re_i(i));
    }

    if(!isNA(y_i(i))){

      if(lik_binomial_i(i) == 1){
        PARALLEL_REGION jnll -= dbinom( y_i(i), n_i(i), invlogit_robust(prob_i(i)), true ) * w_i(i);
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
  }

  return jnll;
}
