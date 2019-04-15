simulate_data_joint_ibm <- function(){

  #######################################################################################
  ##                                                                                   ##
  ## R codes for the following joint model                                             ##
  ##                                                                                   ##
  ## Y_{ij} = \boldsymbol X_{ij} \boldysmbol \alpha + U_i + W_i(t_{ij}) + Z_{ij}       ##
  ## \lambda_i(t) = \lambda_0(t) exp(D_i(t) \beta + \gamma_1 U_i + \gamma_2 W_i(t))    ##
  ##                                                                                   ##
  ## with Weibull baseline hazard and W is integrated random Walk                      ##
  ##									             ##
  ##                                                                                   ##
  #######################################################################################

  ###########################
  ##                       ##
  ## Simulating a data set ##
  ##                       ##
  ###########################

  ## packages needed
  library(mvnfast)
  library(mvtnorm)
  library(MASS)

  ## follow-up period
  fu_period <- seq(0, 5, by = 0.1)

  ## number of subjects
  n_subj <- 250

  ## random intercept, integrated Random-Walk and measurement error variances
  varpar_true <- c(0.11, 0.001, 0.05)

  ## longitudinal model fixed effects
  alpha1_true <- matrix(c(4.6, -0.05, -0.1), ncol = 1)

  ## shape parameter of Weibull
  k_true <- 1.2

  ## survival model fixed effects, first element is the logarithm of the scale parameter for Weibull
  alpha2_true <- matrix(c(log(0.01), 0.001), ncol = 1)

  ## association parameters
  gamma_true  <- matrix(c(-0.3, -0.5))

  ## empty longitudinal and survival datasets
  longitudinal_data <- survival_data <- c()

  for(i in 1 : n_subj){

    bage_i <- sample(20 : 80, 1)

    X2i_sim <- matrix(c(1, bage_i), nrow = length(fu_period), ncol = 2, byrow = T)

    Ui_sim <- rnorm(1, 0, sqrt(varpar_true[1]))

    cov_W  <- varpar_true[2] *
      outer(fu_period, fu_period, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))

    Wi_sim <- matrix(rmvnorm(1, rep(0, ncol(cov_W)), cov_W), ncol = 1)

    surv_prob_save <- c()

    for(ii in 2 : length(fu_period)){

      inside_left  <- k_true * fu_period[2 : ii]^(k_true - 1)
      inside_right <- exp(as.numeric(X2i_sim[2 : ii, ] %*% alpha2_true) +
                            Ui_sim * gamma_true[1] +
                            Wi_sim[2 : ii, ] * gamma_true[2])

      inside <- inside_left * inside_right

      surv_prob_save <- c(surv_prob_save, exp(-sum(inside)))

    }#for(ii in 1 : length(fui))

    surv_prob_save <- c(1, surv_prob_save)

    failure_prob <- 1 - surv_prob_save

    unif_rv <- runif(1)

    if(unif_rv <= max(failure_prob)) ind <- which((failure_prob >= unif_rv) == TRUE)[1]
    if(unif_rv >  max(failure_prob)) ind <- length(failure_prob)

    Ti_sim <- fu_period[ind]

    if(ind < length(fu_period)) deltai_sim <- 1

    if(ind == length(fu_period)){
      if(unif_rv <= max(failure_prob)) deltai_sim <- 1
      if(unif_rv >  max(failure_prob)) deltai_sim <- 0
    }

    fu_i        <- fu_period[fu_period <= min(4, Ti_sim)]
    length_fu_i <- length(fu_i)

    W_i <- Wi_sim[1 : length_fu_i]

    X1i_sim <- cbind(1, rep(bage_i, length_fu_i), fu_i)

    Zi_sim <- rnorm(length_fu_i, 0, sqrt(varpar_true[3]))

    Yi_sim <- X1i_sim %*% alpha1_true + rep(Ui_sim, length_fu_i) + W_i + Zi_sim

    longitudinal_data_i <- cbind(rep(i, length_fu_i), X1i_sim[, -1, drop = F], Yi_sim)
    longitudinal_data   <- rbind(longitudinal_data, longitudinal_data_i)

    survival_data_i <- cbind(i, X2i_sim[1, -1, drop = F], Ti_sim, deltai_sim)
    survival_data   <- rbind(survival_data, survival_data_i)

  }#for(i in 1 : n)

  longitudinal_data <- as.data.frame(longitudinal_data)
  survival_data     <- as.data.frame(survival_data)

  colnames(longitudinal_data) <- c("ID", "bage", "fu", "log.eGFR")
  colnames(survival_data) <- c("ID",  "bage", "stime", "death")

}
