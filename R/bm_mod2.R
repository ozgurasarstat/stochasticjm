bm_mod2 = "

//functions{

//matrix create_covmat_bm(vector t, int mi, real sigmasq_W){

//matrix[mi, mi] out;

//for (i in 1:(mi-1)){
//out[i, i] = sigmasq_W * t[i];
//for (j in (i+1):mi){
//out[i, j] = sigmasq_W * fmin(t[i], t[j]);
//out[j, i] = out[i, j];
//}
//}
//out[mi, mi] = sigmasq_W * t[mi];

//return out;

//}

//}

data{
int<lower = 1> ntot;
vector[ntot] y;
int<lower = 1> p;
int<lower = 1> q;
int<lower = 1> ngroup;
matrix[ntot, p] x;
matrix[ntot, q] d;
vector[ntot] locs;
int nrepeat[ngroup];
int ind[ngroup, 2];
vector[5] priors;
}

transformed data{
vector[q] zero_B = rep_vector(0.0, q);
}

parameters{
vector[p] alpha;
matrix[ngroup, q] B;
//vector[ntot] W;
corr_matrix[q] Omega;
vector<lower = 0>[q] sigma_B;
//real<lower = 0> sigma_W;
real<lower = 0> sigma_Z;
}

transformed parameters{
cov_matrix[q] Sigma;
vector[ntot] linpred;
vector[ntot] d_B;
//real sigmasq_W = square(sigma_W);

for(i in 1:ngroup){
d_B[ind[i, 1]:ind[i, 2]] =
      to_vector(d[ind[i, 1]:ind[i, 2], ] * to_matrix(B[i, ], q, 1));
}

//linpred = x * alpha + d_B + W;
linpred = x * alpha + d_B;

Sigma = quad_form_diag(Omega, sigma_B);

}

model{

alpha ~ cauchy(0, priors[1]);
Omega ~ lkj_corr(priors[2]);
sigma_B ~ cauchy(0, priors[3]);
//sigma_W ~ cauchy(0, priors[4]);
sigma_Z ~ cauchy(0, priors[5]);

for(i in 1:ngroup){
B[i] ~ multi_normal(zero_B, Sigma);
//W[ind[i, 1]:ind[i, 2]] ~
//  multi_normal(rep_vector(0.0, nrepeat[i]),
//    create_covmat_bm(locs[ind[i, 1]:ind[i, 2]], nrepeat[i], sigma_W));
}

y ~ normal(linpred, sigma_Z);

}

"
