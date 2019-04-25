fit_jm <- function(fixed_long,
                   random_long,
                   fixed_surv,
                   data_long,
                   data_surv,
                   id_long,
                   id_surv,
                   #model,
                   timeVar,
                   #deriv = NULL,
                   bh = "weibull",
                   #bh_nknots = 2, #number of knots for baseline hazard
                   #spline_tv = list("time", 2), #spline for tv dof - same in fit_ld
                   Q = 15,
                   priors = list(),
                   ...){

  if(length(priors) != 10){

      priors_full <- list(alpha = 5,
                          Omega = 2,
                          sigma_B = 5,
                          sigma_W = 5,
                          phi = 5,
                          sigma_Z = 5,
                          log_lambda = 5,
                          log_nu = 5,
                          omega = 5,
                          eta = 5)
      for(i in 1:10){
        if(!(names(priors_full)[i] %in% names(priors))){
          priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
        }
      }
    }

  ## number of subjects and repeated measures
  ngroup  <- data_long[, id_long] %>% unique %>% length
  nrepeat <- data_long[, id_long] %>% table %>% as.numeric

  ## time variable
  locs_list <- split(data[, timeVar], data[, id_long])

  ## be sure that id is: 1, 2, 3, ...
  data_long[, id_long] <- rep(1:ngroup, nrepeat)
  #l_id <- data_long[, id_long]

  ## creat survival data
  data_surv[, id_surv] <- rep(1:nrow(data_surv))
  #s_id <- data_surv[, id_surv]

  ## x and y matrices
  x <- model.matrix(fixed_long, data_long)
  y <- model.frame(fixed_long, data_long)[, 1]

  ## create blok-diagonal random effects design matrix
  # dmat <- model.matrix(random_long, data_long)
  # id_dmat <- data.frame(l_id, dmat)
  # id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  # d <- do.call(magic::adiag, id_dmat_list)
  d <- model.matrix(random_long, data_long)

  ## total number of observations in the longitudinal data
  ## number of covariates in the x and d matrices
  ntot <- nrow(x)
  p <- ncol(x)
  q <- ncol(d)#ncol(id_dmat) - 1

  ## Gauss - Legendre weights and abscissas
  gl_quad <- statmod::gauss.quad(Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes

  ## total number of observations for quadrature approx.
  ntot_quad <- ngroup * Q

  ## extract survival times and event indicator
  mf_surv <- model.frame(fixed_surv, data_surv)
  S <- mf_surv[, 1][, 1]
  E <- mf_surv[, 1][, 2]

  ## calculate times for hazard function for quadrature approx
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))

  ## fixed effects for survival sub-model
  c <- model.matrix(fixed_surv, data_surv)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))

  ## x and d matrices at survival times and quadrature times
  data_long_base <- data_long[!duplicated(l_id), ]
  data_long_quad <- data_long_base[rep(1:ngroup, times = Q), ]

  data_long_base[, timeVar] <- S
  data_long_quad[, timeVar] <- t_quad

  x_T    <- model.matrix(fixed_long, data_long_base)
  x_quad <- model.matrix(fixed_long, data_long_quad)

  d_T <- model.matrix(random_long, data_long_base)

  d_quad <- model.matrix(random_long, data_long_quad)

  ## extend the weights for quadrature approx.
  wt_quad <- rep(wt, ngroup)

  ## prior hyperparameters
  priors_long <- unlist(priors)[1:6]
  priors_long <- unlist(priors)[7:10]

  # prepare a matrix of indices to select rows of d in for loop in stan
  cumsum_nrepeat <- cumsum(nrepeat)
  d_ind <- cbind(c(1, (cumsum_nrepeat[-ngroup] + 1)), cumsum_nrepeat)

  # prepare a matrix of indices for matrices for quadratures
  Q_ind <- cbind((0:(ngroup-1))*Q+1, (1:ngroup)*Q)

  #combine and t_quad and locs
  t_quad_locs_comb <- c()
  for(i in 1:ngroup){
    t_quad_locs_comb <- c(t_quad_locs_comb, c(t_quad[((i-1)*Q+1):(i*Q)]), locs_list[[i]])
  }

  nW <- length(t_quad_locs_comb)

  t_quad_locs_ind <- d_ind + Q
  t_quad_locs_ind[1, 1] <- 1

  ## prepare data as a list to be passed to stan

  data_stan <- list(ntot = ntot,
                    id = l_id,
                    y = y,
                    p = p,
                    q = q,
                    ngroup = ngroup,
                    nrepeat = nrepeat,
                    x = x,
                    d = d,
                    d_ind = d_ind,
                    Q_ind = Q_ind,
                    priors_long = priors_long,
                    priors_surv = priors_surv,
                    Q = Q,
                    ntot_quad = ntot_quad,
                    S = S,
                    E = E,
                    ncol_c = ncol_c,
                    c = c,
                    c_quad = c_quad,
                    x_T = x_T,
                    x_quad = x_quad,
                    d_T = d_T,
                    d_quad = d_quad,
                    wt_quad = wt_quad,
                    t_quad = t_quad,
                    t_quad_locs_ind = t_quad_locs_ind,
                    nW = nW,
                    t_quad_locs_comb = t_quad_locs_comb
                    )

  res <- stan(model_code = exp_cor_jm, data = data_stan, ...)

  return(res)

  }
