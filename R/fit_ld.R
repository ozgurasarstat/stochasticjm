fit_ld <- function(fixed,
                   random,
                   data,
                   id,
                   timeVar,
                   process,
                   priors,
                   ...){

  ## number of subjects and repeated measures
  ngroup  <- data[, id] %>% unique %>% length
  nrepeat <- data[, id] %>% table %>% as.numeric

  ## be sure that id is: 1, 2, 3, ...
  data[, id] <- rep(1:ngroup, nrepeat)
  ids <- data[, id]

  ## time points
  locs <- data[, timeVar]

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

  ## prepare the data
  data_stan <- list(ntot = ntot,
                    id = ids,
                    y = y,
                    p = p,
                    q = q,
                    ngroup = ngroup,
                    nrepeat = nrepeat,
                    x = x,
                    d = d,
                    locs = locs,
                    ind = indices,
                    priors = priors)

  res <- stan(model_code = bm_mod, data = data_stan, ...)

}
