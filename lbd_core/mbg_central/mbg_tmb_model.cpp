// ///////////////////////////////////////////////////
// Roy Burstein and Aaron Osgood-Zimmerman
// August 2017
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
// 5. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
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
  int GEO_random;
  int useGP;
  // Way easier to read these in as vectors of integers and then index the first
  option_list(SEXP x){
    use_priors = asVector<int>(getListElement(x,"use_priors"))[0];
    adreport_off = asVector<int>(getListElement(x,"adreport_off"))[0];
    nugget = asVector<int>(getListElement(x,"nugget"))[0];
    country_random = asVector<int>(getListElement(x,"country_random"))[0];
    NID_random = asVector<int>(getListElement(x,"NID_random"))[0];
    GEO_random = asVector<int>(getListElement(x,"GEO_random"))[0];
    useGP = asVector<int>(getListElement(x,"useGP"))[0];
  }
};

// how to read in an object that is used as a prior
template<class Type>
struct prior_type {
    int logNormal;
    Type par1;
    Type par2;

    prior_type(SEXP x){
        logNormal = asVector<int>(getListElement(x,"logNormal"))[0];
        par1 = asVector<float>(getListElement(x,"par1"))[0];
        par2 = asVector<float>(getListElement(x,"par2"))[0];
    }
};

// evaluate a prior using the read in object
template<class Type>
Type eval_prior(prior_type<Type> prior, Type log_sigma){
    Type penalty;

    if(prior.logNormal == 1){
        penalty = dnorm(log_sigma, prior.par1, prior.par2, true);
    }

    if(prior.logNormal == 0){
        penalty = dgamma(exp(log_sigma), prior.par1, prior.par2, true);
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
  DATA_IVECTOR(geo_re_i);    // GEO identifiers
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

  // Options
  DATA_STRUCT(options, option_list);  // boolean vector of options to be used to select different models/modelling options: see above
  DATA_STRUCT(prior_log_nugget_sigma, prior_type);
  DATA_STRUCT(prior_log_cre_sigma, prior_type);
  DATA_STRUCT(prior_log_nidre_sigma, prior_type);
  DATA_STRUCT(prior_log_geore_sigma, prior_type);
  
  // Parameters
  PARAMETER_VECTOR(alpha_j);   // fixed effect coefs, including intercept as first index
  PARAMETER(logtau);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(trho);             // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho);             // Z autocorrelation parameter for AR1, natural scale
  PARAMETER(log_nugget_sigma); // log of the standard deviation of the normal error nugget term
  PARAMETER(log_cre_sigma);    // log of the standard deviation of the country random effect (later do as vec if using Random SLOPE TODO)
  PARAMETER(log_nidre_sigma);  // log of the standard deviation of the nid random effect
  PARAMETER(log_geore_sigma);  // log of the standard deviation of the geo random effect
  PARAMETER(log_gauss_sigma);  // log of sigma for any gaussian observations

  // Random effects
  PARAMETER_ARRAY(Epsilon_stz);  // Random effects for each STZ mesh location. Should be 3D array of dimensions num_s by num_t by num_z
  PARAMETER_VECTOR(nug_i);       // Random effects of the nugget
  PARAMETER_VECTOR(cntry_re);    // Random effects values for country intercept
  PARAMETER_VECTOR(nid_re);      // Random effects values for nid intercept
  PARAMETER_VECTOR(geo_re);      // Random effects values for nid intercept

  printf("Epsilon_stz size: %ld \n", Epsilon_stz.size());

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

  // Make spatial, t and z precision matrices
  SparseMatrix<Type> Q_ss   = spde_Q(logkappa, logtau, M0, M1, M2);
//  SparseMatrix<Type> Q_t    = ar_Q(num_t, trho, Type(1.));
//  SparseMatrix<Type> Q_z    = ar_Q(num_z, zrho, Type(1.));
  printf("Q_ss size: %ld \n", Q_ss.size());

  // Make transformations of some of our parameters
  Type range         = sqrt(8.0) / exp(logkappa);
  Type sigma         = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
  Type trho_trans    = (exp(trho) - 1) / (exp(trho) + 1); // TRANSOFRM from -inf, inf to -1, 1.. //log((1.1 + trho) / (1.1 - trho));
  Type zrho_trans    = (exp(zrho) - 1) / (exp(zrho) + 1); //TRANSOFRM from -inf, inf to -1, 1.. // log((1.1 + zrho) / (1.1 - zrho));
  Type nugget_sigma  = exp(log_nugget_sigma);
  Type cre_sigma     = exp(log_cre_sigma);
  Type nidre_sigma   = exp(log_nidre_sigma);
  Type geore_sigma   = exp(log_geore_sigma);
  Type gauss_sigma   = exp(log_gauss_sigma);


  // Define objects for derived values
  vector<Type> fe_i(num_i);                         // main effect X_ij %*% t(alpha_j)
  vector<Type> epsilon_stz(num_s * num_t * num_z);  // Epsilon_stz unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);                // value of gmrf at data points
  vector<Type> prob_i(num_i);                       // Logit estimated prob for each point i
  vector<Type> nid_re_star(nid_re.size() + 1);      // adjusted nid re effects
  vector<Type> geo_re_star(geo_re.size() + 1);      // adjusted nid re effects

  nid_re_star[0] = Type(0.0);
  for (int i = 1; i < int(nid_re_star.size()); i++){
    nid_re_star[i] = nid_re[i-1];
  }
  
  geo_re_star[0] = Type(0.0);
  for (int i = 1; i < int(geo_re_star.size()); i++){
    geo_re_star[i] = geo_re[i-1];
  }

  // define values for sums of random effects to be able to ensure sum to one
  Type nug_i_sum    = nug_i.sum();
  Type cntry_re_sum = cntry_re.sum();
  Type nid_re_sum   = nid_re.sum();
  Type geo_re_sum   = geo_re.sum();

  // Prior contribution to likelihood. Values are defaulted (for now). Only run if options[0]==1
  if(options.use_priors == 1) {
  // TODO: CHECK THESE HYPERPRIORS ALIGN WITH R-INLA IMPLEMENTATION
   PARALLEL_REGION jnll -= dnorm(logtau,    Type(0.0), Type(1.0),   true);  // N(0,1) prior for logtau
   PARALLEL_REGION jnll -= dnorm(logkappa,  Type(0.0), Type(1.0),   true);  // N(0,1) prior for logkappa
   if(num_t > 1) {
     PARALLEL_REGION jnll -= dnorm(trho, Type(1.0), Type(1), true);  // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   }
   if(num_z > 1) {
     PARALLEL_REGION jnll -= dnorm(zrho, Type(1.0), Type(1), true);  // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   }
   for( int j = 0; j < alpha_j.size(); j++){
     PARALLEL_REGION jnll -= dnorm(alpha_j(j), Type(0.0), Type(3), true); // N(0, sqrt(1/.001)) prior for fixed effects.
   }
   // if using nugget (option in 3rd index)
   if(options.nugget == 1 ){
     // TODO ADD NUGGET PRIOR, MAKE SURE SAME AS R-INLA IMPLEMENTATION
     PARALLEL_REGION jnll -= eval_prior(prior_log_nugget_sigma, log_nugget_sigma);
   }
   // if using nugget (option in 4th index)
   if(options.country_random == 1 ){
     PARALLEL_REGION jnll -= eval_prior(prior_log_cre_sigma, log_cre_sigma);
   }
   // if using nugget (option in 5th index)
   if(options.NID_random == 1 ){
     PARALLEL_REGION jnll -= eval_prior(prior_log_nidre_sigma, log_nidre_sigma);
   }
   if(options.GEO_random == 1 ){
     PARALLEL_REGION jnll -= eval_prior(prior_log_geore_sigma, log_geore_sigma);
   }
  }


  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: S, ST, SZ, and STZ
  if (num_t == 1 & num_z == 1)  {
    printf("GP FOR SPACE  ONLY \n");
    PARALLEL_REGION jnll += GMRF(Q_ss,false)(epsilon_stz);
  } else if(num_t > 1 & num_z == 1) {
    printf("GP FOR SPACE-TIME \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(trho_trans),GMRF(Q_ss,false))(Epsilon_stz);
  //  PARALLEL_REGION jnll += SEPARABLE(GMRF(Q_t,false),GMRF(Q_ss,false))(Epsilon_stz);
  } else if (num_t == 1 & num_z > 1) {
    printf("GP FOR SPACE-Z \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_trans),GMRF(Q_ss,false))(Epsilon_stz);
    //PARALLEL_REGION jnll += SEPARABLE(GMRF(Q_z,false),GMRF(Q_ss,false))(Epsilon_stz);
  } else if (num_t > 1 & num_z > 1) {
    printf("GP FOR SPACE-TIME-Z \n");
    PARALLEL_REGION jnll += SEPARABLE(AR1(zrho_trans),SEPARABLE(AR1(trho_trans),GMRF(Q_ss,false)))(Epsilon_stz);
    //PARALLEL_REGION jnll += SEPARABLE(GMRF(Q_z,false),SEPARABLE(GMRF(Q_t,false),GMRF(Q_ss,false)))(Epsilon_stz);
  }
  

  // nugget contribution to the likelihood
  if(options.nugget == 1 ){
    printf("Nugget \n");
    for (int i = 0; i < num_i; i++){
      // binomial models with sd_i, the additional variance gets put into the nugget.
      // for gaussian models nuggets seem unidentifiable so, the sd_i is in the data likelihood
      PARALLEL_REGION jnll -= dnorm(nug_i(i), Type(0.0), sqrt( pow(sd_i(i),2) + pow(nugget_sigma,2) ), true);
    }
    // soft fix to sum to 1 (should help with identifiability issues)
  //  PARALLEL_REGION jnll -= dnorm(nug_i_sum, Type(0.0), Type(0.0001), true);
  }


  // country random intercept
  if(options.country_random == 1 ){
    printf("Country RE \n");
    for(int i=0; i<cntry_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(cntry_re(i), Type(0.0), cre_sigma, true);
    }
    // soft fix to sum to 1
    // PARALLEL_REGION jnll -= dnorm(cntry_re_sum, Type(0.0), Type(0.0001), true);
  }


  // nid random intercept
  if(options.NID_random == 1 ){
    printf("NID RE \n");
    for(int i=0; i<nid_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(nid_re(i), Type(0.0), nidre_sigma, true);
    }
    // soft fix to sum to 1
   // PARALLEL_REGION jnll -= dnorm(nid_re_sum, Type(0.0), Type(0.0001), true);
  }
  
  // nid random intercept
  if(options.GEO_random == 1 ){
    printf("GEO RE \n");
    for(int i=0; i<geo_re.size(); i++){
      PARALLEL_REGION jnll -= dnorm(geo_re(i), Type(0.0), geore_sigma, true);
    }
    // soft fix to sum to 1
    // PARALLEL_REGION jnll -= dnorm(nid_re_sum, Type(0.0), Type(0.0001), true);
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

  
  // soft-fix sum of epsilons to be 1
//  Type projepsilon_i_sum   = projepsilon_i.sum();
//  PARALLEL_REGION jnll    -= dnorm(projepsilon_i_sum, Type(0.0), Type(0.0001), true);
  
  // evaluate fixed effects for alpha_j values
  vector<Type> calpha_j = constrain_pars(alpha_j, fconstraints);
  fe_i = X_ij * calpha_j.matrix();


  // Return un-normalized density on request
  if (flag == 0) return jnll;

  // Likelihood contribution from each datapoint i
  printf("Data likelihood \n");
  for (int i = 0; i < num_i; i++){

    // mean model
    if(options.useGP==1){
      prob_i(i) = fe_i(i) + projepsilon_i(i) + nug_i(i) + cntry_re(c_re_i(i)) + nid_re_star(nid_re_i(i)) + geo_re_star(geo_re_i(i));
    } else {
      prob_i(i) = fe_i(i) +  nug_i(i) + cntry_re(c_re_i(i)) + nid_re_star(nid_re_i(i)) + geo_re_star(geo_re_i(i));
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
 //   ADREPORT(Epsilon_stz);
  }

  return jnll;
}
