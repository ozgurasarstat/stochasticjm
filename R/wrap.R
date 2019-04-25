wrap <- function(fixed,
                   random,
                   data,
                   id,
                   model,
                 timeVar,
                   priors = list(),
                   ...){

   if(model == "mixed"){
     priors_full <- list(alpha = 5,
                         Omega = 2,
                         sigma_B = 5,
                         sigma_Z = 5)
     for(i in 1:4){
       if(!(names(priors_full)[i] %in% names(priors))){
         priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
       }
     }
   }else if(model %in% c("exp", "bm")){
     priors_full <- list(alpha = 5,
                         Omega = 2,
                         sigma_B = 5,
                         sigma_W = 5,
                         phi = 5,
                         sigma_Z = 5)
     for(i in 1:5){
       if(!(names(priors_full)[i] %in% names(priors))){
         priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
       }
     }
   }

    ## number of subjects and repeated measures
    ngroup  <- data[, id] %>% unique %>% length
    nrepeat <- data[, id] %>% table %>% as.numeric

    ## be sure that id is: 1, 2, 3, ...
    data[, id] <- rep(1:ngroup, nrepeat)

    ## time points
    if(model %in% c("bm", "exp")){
      locs <- data[, timeVar]
    }

    ## x and y matrices
    x <- model.matrix(fixed, data)
    y <- model.frame(fixed, data)[, 1]

    ## create random effects design matrix
    d <- model.matrix(random, data)

    ## total number of observations in the longitudinal data
    ## number of covariates in the x and d matrices
    ntot <- nrow(x)
    p <- ncol(x)
    q <- ncol(d)

    # prepare a matrix of indices to select rows of d in for loop in stan
    cumsum_nrepeat <- cumsum(nrepeat)
    indices <- cbind(c(1, (cumsum_nrepeat[-ngroup] + 1)), cumsum_nrepeat)

    ## Fit the mixed model
    if(model == "mixed"){
      dat_mixed_mod <- list(ntot = ntot,
                            y = y,
                            x = x,
                            d = d,
                            p = p,
                            q = q,
                            ngroup = ngroup,
                            priors = unlist(priors),
                            ind = indices)
      res <- stan(model_code = mixed_mod, data = dat_mixed_mod, ...)
    }else if(model == "bm"){
      dat_bm <- list(ntot = ntot,
                     y = y,
                     x = x,
                     d = d,
                     p = p,
                     q = q,
                     nrepeat = nrepeat,
                     locs = locs,
                     ngroup = ngroup,
                     priors = unlist(priors),
                     ind = indices)
      res <- stan(model_code = bm_mod, data = dat_bm, ...)
    }else if(model == "exp"){
      dat_exp <- list(ntot = ntot,
                     y = y,
                     x = x,
                     d = d,
                     p = p,
                     q = q,
                     nrepeat = nrepeat,
                     locs = locs,
                     ngroup = ngroup,
                     priors = unlist(priors),
                     ind = indices)
      res <- stan(model_code = exp_cor_mod, data = dat_exp, ...)
    }

  return(res)

}
