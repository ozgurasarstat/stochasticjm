# library(nlme)
# library(rstan)
# library(tidyverse)
# library(lmenssp)
#
#
# data("data.sim.ibm")
# data <- data.sim.ibm
#
# ## fit the normal model
# fit <- wrap(fixed = log.egfr ~ bage + fu,
#             random = ~ fu,
#             data = data,
#             id = "id",
#             model = "exp",
#             chains = 1,
#             timeVar = "fu",
#             cores = 1,
#             iter = 2000,
#             warmup = 1000,
#             control = list(adapt_delta = 0.8))
# print(fit, pars = c("alpha", "Sigma", "sigma_W", "range", "sigma_Z"))
