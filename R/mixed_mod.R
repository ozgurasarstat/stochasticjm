mixed_mod = "

data{
int<lower = 1> ntot;        // total number observations
vector[ntot] y;             // response matrix
int<lower = 1> p;           // number of covariates in the fixed effects design matrix
int<lower = 1> q;           // number of covariates in the random effects design matrix
int<lower = 1> ngroup;      // number of subjects/clusters/groups
matrix[ntot, p] x;          // fixed effects design matrix
matrix[ntot, q] d; // random effects design matrix, block diagonal
vector[4] priors; // prior hyperparameters, order: theta, Omega, sigma_B, sigma_Z
int ind[ngroup, 2];
}

transformed data{
vector[q] zero_B = rep_vector(0, q);
}

parameters{
vector[p] alpha;              // fixed effects coefficients - qr
matrix[ngroup, q] B;          // random effects coefficients
corr_matrix[q] Omega;             // correlation matrix for random effects
vector<lower = 0>[q] sigma_B; // scale parameters for random effects
real<lower = 0> sigma_Z;      // scale parameter of measurement error
}

transformed parameters{
cov_matrix[q] Sigma;
vector[ntot] linpred;
vector[ntot] d_B;

for(i in 1:ngroup){
d_B[ind[i,1]:ind[i,2]] = to_vector(d[ind[i,1]:ind[i,2],] * to_matrix(B[i,],q,1));
}

linpred = x * alpha + d_B;
Sigma = quad_form_diag(Omega, sigma_B);
}

model{
alpha ~ cauchy(0, priors[1]);
for(i in 1:ngroup){
B[i] ~ multi_normal(zero_B, Sigma);
}
Omega ~ lkj_corr(priors[2]);
sigma_B ~ cauchy(0, priors[3]);
sigma_Z ~ cauchy(0, priors[4]);
y ~ normal(linpred, sigma_Z);
}

generated quantities{
real sigmasq;
sigmasq = sigma_Z^2;
}
"
