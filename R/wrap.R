wrap <- function(fixed, 
                   random, 
                   data, 
                   id, 
                   model, 
                   priors = list(), 
                   ...){
  

    priors_full <- list(alpha = 5, 
                        Omega = 2, 
                        sigma_B = 5, 
                        sigma_Z = 5)
    for(i in 1:4){
      if(!(names(priors_full)[i] %in% names(priors))){
        priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
      }
    } 
  
  
    ## number of subjects and repeated measures
    ngroup  <- data[, id] %>% unique %>% length
    nrepeat <- data[, id] %>% table %>% as.numeric
    
    ## be sure that id is: 1, 2, 3, ...
    data[, id] <- rep(1:ngroup, nrepeat)

    ## time points
    #locs <- data[, timeVar]
    
    ## x and y matrices
    x <- model.matrix(fixed, data)
    y <- model.frame(fixed, data)[, 1]
    
    ## create blok-diagonal random effects design matrix
    d <- model.matrix(random, data)
    
    ## total number of observations in the longitudinal data
    ## number of covariates in the x and d matrices
    ntot <- nrow(x)
    p <- ncol(x)
    q <- ncol(d)
    
    # prepare a matrix of indices to select rows of d in for loop in stan
    cumsum_nrepeat <- cumsum(nrepeat)
    indices <- cbind(c(1, (cumsum_nrepeat[-ngroup] + 1)), cumsum_nrepeat)
  
  ## Fit the normal - normal model
    dat_mixed_mod <- list(ntot = ntot,                        
                        y = y,
                        x = x,
                        d = d,
                        p = p,
                        q = q,
                        ngroup = ngroup,
                        priors = unlist(priors),
                        ind = indices
                        )
    res <- stan(model_code = mixed_mod, data = dat_mixed_mod, ...)
    
  return(res)
  
}