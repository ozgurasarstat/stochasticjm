library(stochasticjm)
data(renal_data)

gfr_data <- renal_data$gfr_surv_short$gfr_short
surv_data <- renal_data$gfr_surv_short$surv_short

## fit the normal model
fit <- wrap(fixed = gfr ~ age + years_after1year,
            random = ~ years_after1year,
            data = gfr_data,
            id = "id",
            model = "exp",
            chains = 4,
            timeVar = "years_after1year",
            #priors = list(sigma_W = 1),
            cores = 4,
            iter = 2000,
            warmup = 1000,
            control = list(adapt_delta = 0.8),
            pars = c("alpha", "Sigma", "sigma_W", "phi", "sigma_Z"))

print(fit, pars = c("alpha", "Sigma", "sigma_W", "phi", "sigma_Z"), digits_summary = 2)
traceplot(fit, pars = c("alpha", "Sigma", "sigma_W", "phi", "sigma_Z"))

exp_jm_fit <- wrap_jm(fixed_long = egfr ~ age + years_after1year,
                      random_long = ~ 1,
                      fixed_surv = cbind(fuyears_after1year, failure) ~ age,
                      data_long = gfr_data,
                      data_surv = surv_data,
                      id_long = "id",
                      id_surv = "id",
                      #model,
                      timeVar = "years_after1year",
                      #deriv = NULL,
                      bh = "weibull",
                      #bh_nknots = 2, #number of knots for baseline hazard
                      #spline_tv = list("time", 2), #spline for tv dof - same in fit_ld
                      Q = 15)
