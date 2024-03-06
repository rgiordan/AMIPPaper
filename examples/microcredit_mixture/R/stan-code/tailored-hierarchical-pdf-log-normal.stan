



data {
int M; // number of categories. (negative, zero, positive) = (1, 2, 3)
int N; // number of observations
int P; // dimentionality of beta parameter
int K; // number of sites
int Q; // number of quantiles
int cat[N]; // category indicator
int N_neg; // number of obvs in negative tail
int N_pos; // number of obvs in positive tail
vector[P] x[N]; // covariates, used only to predict category.
int site[N]; // site indicator

// The *_neg and *_pos are the category 1 resp. 3 components of the data,
// in the original order.
int treatment_neg[N_neg]; // treatment in neg tail (category 1)
int treatment_pos[N_pos]; // treatment in pos tail (category 3)
int site_neg[N_neg]; // site indicator in neg tail (category 1)
int site_pos[N_pos]; // site indicator in pos tail (category 3)

// Note that y_neg and y_pos are the absolute values of the negative resp.
// positive responses.
real y_neg[N_neg]; // the negative tail (category 1)
real y_pos[N_pos]; // the positive tail (category 3)

vector[Q] quantile_vec ; // vector of quantiles of interest.  Unused.
}

parameters {

real mu[2];
real tau[2];
real<lower=0> sd_mu[2];
real<lower=0> sd_tau[2];
real sigma_control[2];
real sigma_TE[2];
real<lower=0> sd_sigma_control[2];
real<lower=0> sd_sigma_TE[2];
matrix[K,2] mu_k;
matrix[K,2] tau_k;
matrix[K,2] sigma_control_k;
matrix[K,2] sigma_TE_k;
matrix[M-1,P] beta; // the parent parameters minus the Mth category

// beta_k = sigma * beta_k_raw + beta.  So sigma is a scale, not a variance.
matrix<lower=0>[M,P] sigma; // the set of M*P parent variances (not a covariance matrix)
matrix[M,P] beta_k_raw[K]; // the hierarchical increments

}
transformed parameters{
matrix[M,P] beta_full;
matrix[M,P] beta_k[K];

// The last category is taken to have 0 log odds.
beta_full = append_row(beta, rep_row_vector(0, P));
for (m in 1:M){
 for (k in 1:K){
  for (p in 1:P){
   beta_k[k,m,p] = beta_full[m,p] + sigma[m,p] * beta_k_raw[k,m,p];
}}}
}
model {

// Priors
to_vector(beta) ~ normal(0,5);
to_vector(sigma) ~ cauchy(0,2);
for (m in 1:M){
  for (k in 1:K){
    beta_k_raw[k,m] ~ normal(0,1);
  }}

mu ~ normal(0,100);
tau ~ normal(0,100);
sd_mu ~ cauchy(0,2);
sd_tau ~ cauchy(0,2);
sigma_control ~ normal(0,100);
sd_sigma_control ~ cauchy(0,2);
sigma_TE ~ normal(0,100);
sd_sigma_TE ~ cauchy(0,2);

// Hierachical parameters

// *[sign] = Location for the variable *_k[:, sign]
// sd_*[sign] = Standard deviation for the variable *_k[:, sign]

for (k in 1:K){ // Site
  for (i in 1:2){ // Sign
    mu_k[k,i] ~ normal(mu[i], sd_mu[i]);
    tau_k[k,i] ~ normal(tau[i], sd_tau[i]);
    sigma_control_k[k,i] ~ normal(sigma_control[i], sd_sigma_control[i]);
    sigma_TE_k[k,i] ~ normal(sigma_TE[i], sd_sigma_TE[i]);
  }
}

// beta_k[site] = regressor for category (negative, zero, positive)
// mu_k[site, sign] = Offset for log(y)
// tau_k[site, sign] = Treatment effect for log(y)
// sigma_control_k[site, sign] = log variance offset for log(y)
// sigma_TE_k[site, sign] = log variance treatment effect for log(y)

// Data generation

// This is the only place `cat` and `x` are used.  Note that x[n] is a P-vector
// and beta_k[site[n]] is an M x P matrix, where M=3 is the number of
// categories, and P=2 is the dimension of x.
for (n in 1:N)
  cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);

for (n in 1:N_neg){
    y_neg[n] ~ lognormal(
      mu_k[site_neg[n],1] + tau_k[site_neg[n],1] * treatment_neg[n],
      exp(sigma_control_k[site_neg[n],1] +
          sigma_TE_k[site_neg[n],1] * treatment_neg[n]));
}
for(n in 1:N_pos){
    y_pos[n] ~ lognormal(
      mu_k[site_pos[n],2] + tau_k[site_pos[n],2] * treatment_pos[n],
      exp(sigma_control_k[site_pos[n],2] +
          sigma_TE_k[site_pos[n],2] * treatment_pos[n]));
}
} // Model block
