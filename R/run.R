library(stochasticjm)
library(lmenssp)

data("data.sim.ibm")

data <- data.sim.ibm[unique(data.sim.ibm$id)[1:50], ]

## fit the normal model
fit <- wrap(fixed = log.egfr ~ bage + fu,
            random = ~ fu,
            data = data,
            id = "id",
            model = "exp",
            chains = 1,
            timeVar = "fu",
            priors = list(sigma_W = 1),
            cores = 1,
            iter = 2000,
            warmup = 1000,
            control = list(adapt_delta = 0.99))
print(fit, pars = c("alpha", "Sigma", "sigma_W", "range", "sigma_Z"), digits_summary = 2)
