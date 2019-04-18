bm_mod = "

data{
int<lower = 1> ntot;
vector[ntot] y;
int<lower = 1> p;
int<lower = 1> q;
int<lower = 1> ngroup;
matrix[ntot, p] x;
matrix[ntot, q] d;
int ind[ngroup, 2];
vector[4] priors;
}

transformed data{
vector[q] zero_B = rep_vector(0.0, q);
}

parameters{
vector[p] alpha;
matrix[ngroup, q] B;
corr_matrix[q] Omega;
vector<lower = 0>[q] sigma_B;
real<lower = 0> sigma_Z;
}

transformed parameters{
cov_matrix[q] Sigma;
vector[ntot] linpred;
vector[ntot] d_B;

for(i in 1:ngroup){
d_B[ind[i, 1]:ind[i, 2]] = to_vector(d[ind[i, 1]:ind[i, 2], ] * to_matrix(B[i, ], q, 1));
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

"
