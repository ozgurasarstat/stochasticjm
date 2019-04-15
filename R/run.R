# head(longitudinal_data)
# head(survival_data)
#
# library(JM)
#
# kfit <- survfit(Surv(stime, death) ~ 1, data = survival_data)
# plot(kfit)
#
# lme1 <- lme(log.eGFR ~ bage + fu, random = ~ fu | ID, data = longitudinal_data, method = "ML")
# summary(lme1)
#
# cox1 <- coxph(Surv(stime, death) ~ bage, data = survival_data, x = TRUE)
# summary(cox1)
#
# jm1 <- jointModel(lme1, cox1, timeVar = "fu", method = "weibull-PH-aGH", verbose = TRUE)
# summary(jm1)
#
# dform <- list(fixed = ~ 1, indFixed = 3, random = ~ 1, indRandom = 2)
# jm2 <- update(jm1, parameterization = "both", derivForm = dform)
#
# summary(jm2)
#
# alpha1_true
# sqrt(varpar_true)
# alpha2_true
#
# joint.fit <- joint.ibm(formula.l = log.eGFR ~ bage + fu,
#                        formula.s = cbind(stime, death) ~ bage,
#                        data.l = longitudinal_data, data.s = survival_data,
#                        id.l = longitudinal_data$ID, timeVar = longitudinal_data$fu,
#                        init.em.l = c(c(alpha1_true), varpar_true),
#                        init.em.s = c(k_true, c(alpha2_true), c(gamma_true)),
#                        tol.em = 0.001, maxiter.em = 1000, discrete = c(0.1, 1),
#                        P = 500, tol.nr = 0.001, maxiter.nr = 1000)
