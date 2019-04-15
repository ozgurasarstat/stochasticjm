joint_ibm <- function(formula.l, formula.s, data.l, data.s, id.l, timeVar,
                        init.em.l, init.em.s, tol.em = 0.001, maxiter.em = 1000,
                        discrete, P = 500, tol.nr = 0.001, maxiter.nr = 1000){

    #####
    ##### options of joint.ibm
    #####
    ## formula.l is a one-sided formula for the longitudinal model
    ## formula.s is a one-sided formula for the survival model
    ## data.l is a data frame that contains the longitudinal data
    ## data.s is a data frame that contains the survival data
    ## id.l is a vector containing the id from the longitudinal data
    ## timeVar is a vector containing the time variable from the longitudinal data
    ## discrete is a two-element vector for the precision of discretisation and rounding, respectively
    ## P is the number of Monte-Carlo samples
    ## tol.nr is the tolerance value for the Newton-Raphson procedure for survival parameters
    ## maxiter.nr is the maximum number of iterations for Newton-Raphson

    # id for survival data
    id.s <- unique(id.l)

    # longitudinal covariates, responses, and time variable
    mf.l <- model.frame(formula = formula.l, data = data.l)
    y.l  <- model.extract(mf.l, "response")
    x.l  <- as.data.frame(cbind(id.l, as.matrix(model.matrix(attr(mf.l, "terms"), data = mf.l))))
    colnames(x.l)[1] <- gsub("[[:punct:]]", "", colnames(x.l)[1])

    # converting longitudinal data into lists
    X1M <- split(x.l[, -1], x.l[,1])
    YM  <- tapply(y.l, id.l, function(x) x)
    time.var <- tapply(timeVar, id.l, function(x) x)

    # survival covariates, time and event indicator
    mf.s <- model.frame(formula = formula.s, data = data.s)
    y.s  <- model.extract(mf.s, "response")
    x.s  <- as.data.frame(cbind(id.s, as.matrix(model.matrix(attr(mf.s, "terms"), data = mf.s))))
    colnames(x.s)[1] <- gsub("[[:punct:]]", "", colnames(x.s)[1])

    # converting survival data into lists
    X2M <- split(x.s[, -1], x.s[, 1])
    TM  <- tapply(y.s[, 1], id.s, function(x) x)
    DM  <- tapply(y.s[, 2], id.s, function(x) x)

    # number of longitudinal and survival model covariates, respectively
    ncov1 <- ncol(x.l) - 1
    ncov2 <- ncol(x.s) - 1

    # number of subjects
    n  <- length(unique(id.l))

    # number of repeated measures per patient
    ni.long <- as.numeric(table(id.l))

    # total number of observations
    N  <- sum(ni.long)

    # number of observations after discretising time
    #mi.long <- ni.long

    # total number of observations after discretising time
    #M  <- sum(mi.long)

    # to start the E-M
    iter.em      <- 0

    # initials
    #theta.new <- c(long.res, pp.res)
    theta.new <- c(init.em.l, init.em.s)
    theta.old <- theta.new * 20

    while(sqrt((theta.old - theta.new) %*% (theta.old - theta.new)) > tol.em & maxiter.em >= iter.em){

      theta.old <- theta.new

      alpha1  <- as.matrix(theta.old[1 : ncov1], ncol = 1)
      omegasq <- theta.old[ncov1 + 1]
      sigmasq <- theta.old[ncov1 + 2]
      tausq   <- theta.old[ncov1 + 3]
      k       <- theta.old[ncov1 + 4]
      alpha2  <- as.matrix(theta.old[(ncov1 + 5) : (ncov1 + 5 + ncov2 - 1)], ncol = 1)
      gamma1  <- theta.old[length(theta.old) - 1]
      gamma2  <- theta.old[length(theta.old)]

      #####################
      ### expectations

      ## total number of time points after discretisation
      M <- 0

      sum.expec.alpha1.hat <- sum.expec.tausq.hat <- sum.expec.omegasq.hat <- sum.expec.sigmasq.hat <- 0

      for (i in 1 : n){

        Yi  <- matrix(YM[[i]], ncol = 1)
        X1i <- as.matrix(X1M[[i]])

        Ti <- round(TM[[i]], discrete[2])
        Di <- DM[[i]]
        X2i <- as.matrix(X2M[[i]], nrow = 1)

        t.i      <- round(time.var[[i]], discrete[2])
        t.i.long <- round(seq(0, Ti, discrete[1]), discrete[2])

        M <- M + length(t.i.long)

        ni <- ni.long[i]
        #mi <- mi.long[i]

        cov.W.Y <- sigmasq * outer(t.i.long, t.i, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))
        cov.Y   <- omegasq * matrix(1, ncol = ni, nrow = ni) +
          sigmasq * outer(t.i, t.i, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) +
          tausq * diag(ni)
        cov.Y.inv <- solve(cov.Y)
        cov.W     <- sigmasq * outer(t.i.long, t.i.long, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))

        cov.U.Y <- omegasq * matrix(1, ncol = ni, nrow = 1)
        cov.U   <- omegasq

        ri <- Yi - X1i %*% alpha1

        mean.W.given.Y <- cov.W.Y %*% cov.Y.inv %*% ri
        cov.W.given.Y  <- cov.W - cov.W.Y %*% cov.Y.inv %*% t(cov.W.Y)

        mean.U.given.Y <- cov.U.Y %*% cov.Y.inv %*% ri
        cov.U.given.Y  <- cov.U - cov.U.Y %*% cov.Y.inv %*% t(cov.U.Y)

        expr.alpha1.hat <- expr.tausq.hat <- expr.omegasq.hat <- expr.sigmasq.hat <- expr.bottom <- 0

        for(p in 1 : P){

          W.mc      <- cbind(t.i.long, as.numeric(rmvnorm(1, mean.W.given.Y, cov.W.given.Y)))
          W0.mc     <- matrix(W.mc[W.mc[, 1] %in% t.i, ], ncol = 2)
          U.mc      <- rnorm(1, mean.U.given.Y, sqrt(cov.U.given.Y))
          U.mc.long <- matrix(U.mc, ncol = 1, nrow = nrow(W0.mc))

          expr.surv <- (k * Ti^(k-1) * exp(X2i %*% alpha2 + U.mc * gamma1 + W.mc[nrow(W.mc), -1] * gamma2))^Di *
            exp(-discrete[1] * sum(k * t.i.long[-1]^(k-1) * exp(as.numeric(X2i %*% alpha2) + U.mc * gamma1 + W.mc[-1, -1] * gamma2)))
          expr.surv <- as.numeric(expr.surv)

          expr.alpha1.hat <- expr.alpha1.hat + t(X1i) %*% (Yi - U.mc.long - W0.mc[, -1]) * expr.surv
          expr.tausq.hat  <- expr.tausq.hat +
            t(Yi - X1i %*% alpha1 - U.mc.long - W0.mc[, -1]) %*% (Yi - X1i %*% alpha1 - U.mc.long - W0.mc[, -1]) * expr.surv
          expr.omegasq.hat <- expr.omegasq.hat + U.mc^2 * expr.surv
          expr.sigmasq.hat <- expr.sigmasq.hat + t(W.mc[, -1]) %*% ginv(cov.W/sigmasq) %*% W.mc[, -1] * expr.surv

          expr.bottom <- expr.bottom + expr.surv

        }#P

        sum.expec.alpha1.hat  <- sum.expec.alpha1.hat  + expr.alpha1.hat/expr.bottom
        sum.expec.tausq.hat   <- sum.expec.tausq.hat   + expr.tausq.hat/expr.bottom
        sum.expec.omegasq.hat <- sum.expec.omegasq.hat + expr.omegasq.hat/expr.bottom
        sum.expec.sigmasq.hat <- sum.expec.sigmasq.hat + expr.sigmasq.hat/expr.bottom

      }#nsubj

      x.l.m <- as.matrix(x.l[, -1])

      alpha1.hat  <- solve(t(x.l.m) %*% x.l.m) %*% sum.expec.alpha1.hat
      tausq.hat   <- sum.expec.tausq.hat/N
      omegasq.hat <- sum.expec.omegasq.hat/n
      sigmasq.hat <- sum.expec.sigmasq.hat/M


      ##################################
      ### estimating survival parameters

      surv.par.new <- c(k, alpha2, gamma1, gamma2)
      surv.par.old <- surv.par.new * 20

      iter.nr <- 0

      while(sqrt((surv.par.old - surv.par.new) %*% (surv.par.old - surv.par.new)) > tol.nr & iter.nr <= maxiter.nr){

        surv.par.old <- surv.par.new

        k.hat      <- surv.par.old[1]
        alpha2.hat <- matrix(surv.par.old[2 : (2 + ncov2 - 1)], ncol = 1)
        gamma1.hat <- surv.par.old[2 + ncov2]
        gamma2.hat <- surv.par.old[2 + ncov2 + 1]

        sum.expec.derv.k  <- sum.expec.derv.alpha2 <- sum.expec.derv.gamma1.1 <- sum.expec.derv.gamma1.2 <-
          sum.expec.derv.gamma2.1 <- sum.expec.derv.gamma2.2 <- sum.expec.derv.k.k <- sum.expec.derv.k.alpha2 <-
          sum.expec.derv.k.gamma1 <- sum.expec.derv.k.gamma2 <- sum.expec.derv.alpha2.alpha2 <- sum.expec.derv.alpha2.gamma1 <-
          sum.expec.derv.alpha2.gamma2 <- sum.expec.derv.gamma1.gamma1 <- sum.expec.derv.gamma1.gamma2 <-
          sum.expec.derv.gamma2.gamma2 <- sum.left.derv.k <- sum.left.derv.alpha2 <- sum.left.derv.k.k <- 0

        for (i in 1 : n){

          Yi  <- matrix(YM[[i]], ncol = 1)
          X1i <- as.matrix(X1M[[i]])

          Ti <- TM[[i]]
          Di <- DM[[i]]
          X2i <- as.matrix(X2M[[i]], nrow = 1)

          t.i      <- round(time.var[[i]], discrete[2])
          t.i.long <- round(seq(0, Ti, discrete[1]), discrete[2])

          ni <- ni.long[i]
          #mi <- mi.long[i]

          cov.W.Y <- sigmasq * outer(t.i.long, t.i, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))
          cov.Y   <- omegasq * matrix(1, ncol = ni, nrow = ni) +
            sigmasq * outer(t.i, t.i, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) +
            tausq * diag(ni)
          cov.Y.inv <- solve(cov.Y)
          cov.W     <- sigmasq * outer(t.i.long, t.i.long, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))

          cov.U.Y <- omegasq * matrix(1, ncol = ni, nrow = 1)
          cov.U   <- omegasq

          ri <- Yi - X1i %*% alpha1

          mean.W.given.Y <- cov.W.Y %*% cov.Y.inv %*% ri
          cov.W.given.Y  <- cov.W - cov.W.Y %*% cov.Y.inv %*% t(cov.W.Y)

          mean.U.given.Y <- cov.U.Y %*% cov.Y.inv %*% ri
          cov.U.given.Y  <- cov.U - cov.U.Y %*% cov.Y.inv %*% t(cov.U.Y)

          derv.k.upper <- derv.alpha2.upper <- derv.gamma1.upper1 <- derv.gamma1.upper2 <- derv.gamma2.upper1 <- derv.gamma2.upper2 <-
            derv.k.k.upper <- derv.k.alpha2.upper <- derv.k.gamma1.upper <- derv.k.gamma2.upper <-
            derv.alpha2.alpha2.upper <- derv.alpha2.gamma1.upper <- derv.alpha2.gamma2.upper <-
            derv.gamma1.gamma1.upper <- derv.gamma1.gamma2.upper <- derv.gamma2.gamma2.upper <- 0

          for(p in 1 : P){

            W.mc <- cbind(t.i.long, as.numeric(rmvnorm(1, mean.W.given.Y, cov.W.given.Y)))
            U.mc <- rnorm(1, mean.U.given.Y, sqrt(cov.U.given.Y))

            expr.surv <- (k * Ti^(k-1) * exp(X2i %*% alpha2 + U.mc * gamma1 + W.mc[nrow(W.mc), -1] * gamma2))^Di *
              exp(-discrete[1] * sum(k * t.i.long[-1]^(k-1) * exp(as.numeric(X2i %*% alpha2) + U.mc * gamma1 + W.mc[-1, -1] * gamma2)))
            expr.surv <- as.numeric(expr.surv)

            ## derv_k_upper
            derv.k.upper0 <- sum(discrete[1] * (t.i.long[-1]^(k.hat-1) + k.hat*(k.hat-1)*t.i.long[-1]^(k.hat-2)) *
                                   exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.k.upper  <- derv.k.upper + derv.k.upper0 * expr.surv

            ## derv_alpha2_upper
            derv.alpha2.upper0 <- 0
            for(ii in 2 : length(t.i.long)){
              derv.alpha2.upper0 <- derv.alpha2.upper0 +
                k.hat*(t.i.long[ii])^(k.hat-1)*t(X2i) * as.numeric(exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[ii, -1]*gamma2.hat))
            }
            derv.alpha2.upper <- derv.alpha2.upper + discrete[1] * derv.alpha2.upper0 * expr.surv

            ## derv_gamma1_upper1
            derv.gamma1.upper1 <- derv.gamma1.upper1 + U.mc * expr.surv

            ## derv_gamma1_upper2
            derv.gamma1.upper20 <- sum(discrete[1] * k.hat*t.i.long[-1]^(k.hat-1)*U.mc * exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.gamma1.upper2  <- derv.gamma1.upper2 + derv.gamma1.upper20

            ## derv_gamma2_upper1
            derv.gamma2.upper1 <- derv.gamma2.upper1 + W.mc[nrow(W.mc), -1] * expr.surv

            ## derv_gamma2_upper2
            derv.gamma2.upper20 <- sum(discrete[1] * k.hat*t.i.long[-1]^(k.hat-1)* W.mc[-1, -1] * exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.gamma2.upper2  <- derv.gamma2.upper2 + derv.gamma2.upper20 * expr.surv

            ## derv_k_k_upper
            derv.k.k.upper0 <- sum(discrete[1] * ((k.hat-1)*t.i.long[-1]^(k.hat-2)+(2*k.hat-1)*t.i.long[-1]^(k.hat-2)+k.hat*(k.hat-1)*(k.hat-2)*t.i.long[-1]^(k.hat-3)) *
                                     exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.k.k.upper  <- derv.k.k.upper + derv.k.k.upper0 * expr.surv

            ## derv_k_alpha2_upper
            derv.k.alpha2.upper0 <- 0
            for (ii in 2 : length(t.i.long)){
              derv.k.alpha2.upper0 <- derv.k.alpha2.upper0 +
                (t.i.long[ii]^(k.hat-1)+k.hat*(k.hat-1)*t.i.long[ii]^(k.hat-2)) * t(X2i) * as.numeric(exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[ii, -1]*gamma2.hat))
            }
            derv.k.alpha2.upper <- derv.k.alpha2.upper + discrete[1] * derv.k.alpha2.upper0 * expr.surv

            ## derv_k_gamma1_upper
            derv.k.gamma1.upper0 <- sum(discrete[1] * (t.i.long[-1]^(k.hat-1) + k.hat*(k.hat-1)*t.i.long[-1]^(k.hat-2)) * U.mc *
                                          exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.k.gamma1.upper  <- derv.k.gamma1.upper + derv.k.gamma1.upper0 * expr.surv

            ## derv_k_gamma2_upper
            derv.k.gamma2.upper0 <- sum(discrete[1] * (t.i.long[-1]^(k.hat-1) + k.hat*(k.hat-1)*t.i.long[-1]^(k.hat-2)) * W.mc[-1, -1] *
                                          exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.k.gamma2.upper  <- derv.k.gamma2.upper + derv.k.gamma2.upper0 * expr.surv

            ## derv_alpha2_alpha2_upper
            derv.alpha2.alpha2.upper0 <- 0
            for(ii in 2 : length(t.i.long)){
              derv.alpha2.alpha2.upper0 <- derv.alpha2.alpha2.upper0 +
                k.hat*t.i.long[ii]^(k.hat-1) * t(X2i) %*% X2i * as.numeric(exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[ii, -1]*gamma2.hat))
            }
            derv.alpha2.alpha2.upper <- derv.alpha2.alpha2.upper + discrete[1] * derv.alpha2.alpha2.upper0 * expr.surv

            ## derv_alpha2_gamma1_upper
            derv.alpha2.gamma1.upper0 <- 0
            for(ii in 2 : length(t.i.long)){
              derv.alpha2.gamma1.upper0 <- derv.alpha2.gamma1.upper0 +
                k.hat*t.i.long[ii]^(k.hat-1) * t(X2i) * U.mc * as.numeric(exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[ii, -1]*gamma2.hat))
            }
            derv.alpha2.gamma1.upper <- derv.alpha2.gamma1.upper + discrete[1] * derv.alpha2.gamma1.upper0 * expr.surv

            ## derv_alpha2_gamma2_upper
            derv.alpha2.gamma2.upper0 <- 0
            for(ii in 2 : length(t.i.long)){
              derv.alpha2.gamma2.upper0 <- derv.alpha2.gamma2.upper0 +
                k.hat*t.i.long[ii]^(k.hat-1) * t(X2i) * W.mc[ii, -1] * as.numeric(exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[ii, -1]*gamma2.hat))
            }
            derv.alpha2.gamma2.upper <- derv.alpha2.gamma2.upper + discrete[1] * derv.alpha2.gamma2.upper0 * expr.surv

            ## derv_gamma1_gamma1_upper
            derv.gamma1.gamma1.upper0 <- sum(discrete[1] * k.hat*t.i.long[-1]^(k.hat-1)*U.mc^2 * exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.gamma1.gamma1.upper  <- derv.gamma1.gamma1.upper + derv.gamma1.gamma1.upper0 * expr.surv

            ## derv_gamma1_gamma2_upper
            derv.gamma1.gamma2.upper0 <- sum(discrete[1] * k.hat*t.i.long[-1]^(k.hat-1)* U.mc * W.mc[-1, -1] * exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.gamma1.gamma2.upper  <- derv.gamma1.gamma2.upper + derv.gamma1.gamma2.upper0 * expr.surv

            ## derv_gamma2_gamma2_upper
            derv.gamma2.gamma2.upper0 <- sum(discrete[1] * k.hat*t.i.long[-1]^(k.hat-1)* W.mc[-1, -1]^2 * exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat))
            derv.gamma2.gamma2.upper  <- derv.gamma2.gamma2.upper + derv.gamma2.gamma2.upper0 * expr.surv

            expr.bottom <- expr.bottom + expr.surv

          }#P

          sum.expec.derv.k             <- sum.expec.derv.k             + derv.k.upper/expr.bottom
          sum.expec.derv.alpha2        <- sum.expec.derv.alpha2        + derv.alpha2.upper/expr.bottom
          sum.expec.derv.gamma1.1      <- sum.expec.derv.gamma1.1      + Di * derv.gamma1.upper1/expr.bottom
          sum.expec.derv.gamma1.2      <- sum.expec.derv.gamma1.2      + derv.gamma1.upper2/expr.bottom
          sum.expec.derv.gamma2.1      <- sum.expec.derv.gamma2.1      + Di * derv.gamma2.upper1/expr.bottom
          sum.expec.derv.gamma2.2      <- sum.expec.derv.gamma2.2      + derv.gamma2.upper2/expr.bottom
          sum.expec.derv.k.k           <- sum.expec.derv.k.k           + derv.k.k.upper/expr.bottom
          sum.expec.derv.k.alpha2      <- sum.expec.derv.k.alpha2      + derv.k.alpha2.upper/expr.bottom
          sum.expec.derv.k.gamma1      <- sum.expec.derv.k.gamma1      + derv.k.gamma1.upper/expr.bottom
          sum.expec.derv.k.gamma2      <- sum.expec.derv.k.gamma2      + derv.k.gamma2.upper/expr.bottom
          sum.expec.derv.alpha2.alpha2 <- sum.expec.derv.alpha2.alpha2 + derv.alpha2.alpha2.upper/expr.bottom
          sum.expec.derv.alpha2.gamma1 <- sum.expec.derv.alpha2.gamma1 + derv.alpha2.gamma1.upper/expr.bottom
          sum.expec.derv.alpha2.gamma2 <- sum.expec.derv.alpha2.gamma2 + derv.alpha2.gamma2.upper/expr.bottom
          sum.expec.derv.gamma1.gamma1 <- sum.expec.derv.gamma1.gamma1 + derv.gamma1.gamma1.upper/expr.bottom
          sum.expec.derv.gamma1.gamma2 <- sum.expec.derv.gamma1.gamma2 + derv.gamma1.gamma2.upper/expr.bottom
          sum.expec.derv.gamma2.gamma2 <- sum.expec.derv.gamma2.gamma2 + derv.gamma2.gamma2.upper/expr.bottom

          sum.left.derv.k      <- sum.left.derv.k + Di * (1/k.hat + log(Ti))
          sum.left.derv.alpha2 <- sum.left.derv.alpha2 + Di * t(X2i)
          sum.left.derv.k.k    <- sum.left.derv.k.k + Di/k.hat^2

        }#nsubj

        score.k      <- sum.left.derv.k - sum.expec.derv.k
        score.alpha2 <- sum.left.derv.alpha2 - sum.expec.derv.alpha2
        score.gamma1 <- sum.expec.derv.gamma1.1 - sum.expec.derv.gamma1.2
        score.gamma2 <- sum.expec.derv.gamma2.1 - sum.expec.derv.gamma2.2

        score <- c(score.k, score.alpha2, score.gamma1, score.gamma2)

        inf.k.k           <- -sum.left.derv.k.k - sum.expec.derv.k.k
        inf.k.alpha2      <- -sum.expec.derv.k.alpha2
        inf.k.gamma1      <- -sum.expec.derv.k.gamma1
        inf.k.gamma2      <- -sum.expec.derv.k.gamma2
        inf.alpha2.alpha2 <- -sum.expec.derv.alpha2.alpha2
        inf.alpha2.gamma1 <- -sum.expec.derv.alpha2.gamma1
        inf.alpha2.gamma2 <- -sum.expec.derv.alpha2.gamma2
        inf.gamma1.gamma1 <- -sum.expec.derv.gamma1.gamma1
        inf.gamma1.gamma2 <- -sum.expec.derv.gamma1.gamma2
        inf.gamma2.gamma2 <- -sum.expec.derv.gamma2.gamma2

        inf.mat.row1 <- c(inf.k.k, inf.k.alpha2, inf.k.gamma1, inf.k.gamma2)
        inf.mat.row2 <- cbind(inf.k.alpha2, inf.alpha2.alpha2, inf.alpha2.gamma1, inf.alpha2.gamma2)
        inf.mat.row3 <- c(inf.k.gamma1, inf.alpha2.gamma1, inf.gamma1.gamma1, inf.gamma1.gamma2)
        inf.mat.row4 <- c(inf.k.gamma2, inf.alpha2.gamma2, inf.gamma1.gamma2, inf.gamma2.gamma2)

        inf.mat <- rbind(inf.mat.row1, inf.mat.row2, inf.mat.row3, inf.mat.row4)

        surv.par.new <- as.numeric(surv.par.old - solve(inf.mat) %*% score)

        print(surv.par.new)

      }#while nr

      ###################################################
      ### estimating survival parameters - by nelder-mead

      P2 <- P/5

      surv.par.est <- function(surv.par){

        print(surv.par)

        k.hat      <- surv.par[1]
        alpha2.hat <- matrix(surv.par[2 : (2 + ncov2 - 1)], ncol = 1)
        gamma1.hat <- surv.par[length(surv.par) - 1]
        gamma2.hat <- surv.par[length(surv.par)]


        sum.expr1 <- sum.expr2 <- sum.expr3 <- sum.expr4 <- 0

        for (i in 1 : n){

          Yi  <- matrix(YM[[i]], ncol = 1)
          X1i <- as.matrix(X1M[[i]])

          Ti <- TM[[i]]
          Di <- DM[[i]]
          X2i <- as.matrix(X2M[[i]], nrow = 1)

          t.i      <- round(time.var[[i]], discrete[1])
          t.i.long <- round(seq(0, Ti, discrete[1]), discrete[2])

          ni <- ni.long[i]
          mi <- mi.long[i]

          cov.W.Y <- sigmasq * outer(t.i.long, t.i, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))
          cov.Y   <- omegasq * matrix(1, ncol = ni, nrow = ni) +
            sigmasq * outer(t.i, t.i, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) +
            tausq * diag(ni)
          cov.Y.inv <- solve(cov.Y)
          cov.W     <- sigmasq * outer(t.i.long, t.i.long, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))

          cov.U.Y <- omegasq * matrix(1, ncol = ni, nrow = 1)
          cov.U   <- omegasq

          ri <- Yi - X1i %*% alpha1

          mean.W.given.Y <- cov.W.Y %*% cov.Y.inv %*% ri
          cov.W.given.Y  <- cov.W - cov.W.Y %*% cov.Y.inv %*% t(cov.W.Y)

          mean.U.given.Y <- cov.U.Y %*% cov.Y.inv %*% ri
          cov.U.given.Y  <- cov.U - cov.U.Y %*% cov.Y.inv %*% t(cov.U.Y)

          U.upper <- W.upper <- integral.upper <- bottom <- 0

          for(p in 1 : P2){

            W.mc  <- cbind(t.i.long, as.numeric(rmvnorm(1, mean.W.given.Y, cov.W.given.Y)))
            #W0.mc <- matrix(W.mc[W.mc[, 1] %in% t.i, ], ncol = 2)
            U.mc  <- rnorm(1, mean.U.given.Y, sqrt(cov.U.given.Y))
            #U.mc.long <- matrix(U.mc, ncol = 1, nrow = nrow(W0.mc))

            #Ki  <- cbind(U.mc, W.mc[-1, -1])

            expr.surv <- (k * Ti^(k-1) * exp(X2i %*% alpha2 + U.mc * gamma1 + W.mc[nrow(W.mc), -1] * gamma2))^Di *
              exp(-discrete[1] * sum(k * t.i.long[-1]^(k-1) * exp(as.numeric(X2i %*% alpha2) + U.mc * gamma1 + W.mc[-1, -1] * gamma2)))
            expr.surv <- as.numeric(expr.surv)

            U.upper <- U.upper + U.mc * expr.surv
            W.upper <- W.upper + W.mc[nrow(W.mc), -1]* expr.surv

            integral.upper0 <- discrete[1] * sum(k.hat*t.i.long[-1]^(k.hat-1) *exp(X2i %*% alpha2.hat + U.mc * gamma1.hat + W.mc[-1, -1]*gamma2.hat)) * expr.surv
            integral.upper  <- integral.upper + integral.upper0

            bottom <- bottom + expr.surv

          }#P

          sum.expr2 <- sum.expr2 + Di * U.upper/bottom
          sum.expr3 <- sum.expr3 + Di * W.upper/bottom
          sum.expr4 <- sum.expr4 + integral.upper/bottom

          sum.expr1 <- sum.expr1 + Di * (log(k.hat) + (k.hat - 1) * log(Ti) + X2i %*% alpha2.hat)

        }#nsubj

        -(sum.expr1 + gamma1.hat * sum.expr2 + gamma2.hat * sum.expr3 - sum.expr4)

      }#surv.par.est

      #fit <- optim(par = c(k, alpha2, gamma1, gamma2), fn = surv.par.est, method = "Nelder-Mead", control = list(maxit = 30))
      #surv.par.hat <- fit$par

      theta.new <- c(alpha1.hat, omegasq.hat, sigmasq.hat, tausq.hat, surv.par.hat)

      iter.em <- iter.em + 1

      print(theta.new)

    }#while.em

    theta.new

  }#joint.ibm
