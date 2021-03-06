#### template_gp_warp.R ####

#' Template Gaussian Process Warping
#'
#' Performs hierarchical curve warping by modeling the warping functions with GP.
#'
#' @param y_list list of curves with equal length
#' @param feat_list list of feature vectors with equal length
#' @param template_feats a feature vector with length mathing those in feat_list
#' @param niter,nburn MCMC iterations and burns. Total run time is the sum.
#' @param int_p number of internal knots specified for the data curves
#' @param asig,bsig,at,bt,al,bl hyperparameters for inverse gamma priors on model variances
#' @param aa,ba hyperparameters for the uniform prior on inverse range parameter
#' @param tune tuning parameter for the proposal normal distribution of alpha.
#' @param progress logical; indicates whether or not to print percent progress.
#' @param debug logical; stops function right before return and enters debug mode.
#' @return a list of warped curves, the estimated mean curve, and a vector of MCMC
#'   acceptance rates. Curve list elements are vectors with length matching
#'   the elements of y_list.
#' @export

template_gp_warp <- function(y_list, feat_list, template_feats,
                             niter = 1000, nburn = 1000, int_p = 20,
                             asig = .1, bsig = .1, at = .1, bt = .1, al = .1, bl = .1,
                             aa = 1, ba = 30, be = 2,tune = 2, progress = TRUE, debug = FALSE){

  #----- Fixed Values -----#

  n <- length(y_list)
  m <- unique(sapply(y_list,length))
  r <- 1 + 3
  if(length(m) != 1) stop("All y_list elements must have equal length")

  time <- seq(0, 1, len = m)
  y_vec <- unlist(y_list)
  template_inds <- c(1, template_feats, m)
  template_feats <- c(time[1], time[template_feats], time[m])
  feat_inds <- lapply(1:n, function(i) c(1, feat_list[[i]], m))
  feat_list <- lapply(1:n, function(i) c(time[1], time[feat_list[[i]]], time[m]))

  p <- int_p + r
  knot_loc_p <- seq(time[1], time[m], len = int_p+2)[-c(1,int_p+2)]
  Hp <- cbs(time, int_p)

  P <- K1(p)
  P[1,1] <- 2

  mb <- rep(0,p)
  diag_nm <- diag(n*m)

  accepts <- 0 #numeric(n)

  nrun <- nburn + niter

  beta_save <- matrix(NA, nrow = nrun, ncol = p)
  beta_save[1,] <- mb

  sig2_save <- numeric(nrun)
  sig2_save[1] <- 1
  tau2_save <- numeric(nrun)
  tau2_save[1] <- 1
  #lam2_save <- lapply(1:n, function(i) numeric(nrun))
  #for(i in 1:n) lam2_save[[i]][1] <- 1
  lam2_save <- numeric(nrun)
  lam2_save[1] <- 1
  #alpha_save <- lapply(1:n, function(i) numeric(nrun))
  #for(i in 1:n) alpha_save[[i]][1] <- (aa + ba)/2
  alpha_save <- numeric(nrun)
  alpha_save[1] <- (aa + ba)/2

  eta2_save <- matrix(NA, nrow = nrun, ncol = n)
  eta2_save[1,] <- be/2

  wtime <- otime <- lapply(1:n, function(i) time)
  wtime_save <- otime_save <- lapply(1:n, function(i){
    out <- matrix(NA, nrow = nrun, ncol = m)
    out[1,] <- time
    return(out)
  })
  H_list <- lapply(1:n, function(i) Hp)
  H_stack <- stack_Matrix(H_list)

  create_M <- function(x_pts, y_pts = x_pts, alpha){
    fields::Matern(fields::rdist(x_pts, y_pts), alpha = alpha, smoothness = 5.1)}

  M_list <- lapply(1:n, function(i) create_M(x_pts = feat_list[[i]], alpha = (aa + ba)/2))
  Minv_list <- lapply(M_list, function(m) chol2inv(chol(m)))


  #----- MCMC Loop -----#

  if(progress == TRUE) bar <- txtProgressBar(min = 2, max = nrun, style = 3)
  for(it in 2:nrun){
    if (progress == TRUE) {
      setTxtProgressBar(bar, it)
    }

    #-- Update alpha, M, otime, wtime --#
    #for(i in 1:n){
    #  current_llik <- try(dmnorm(y = template_feats, mu = feat_list[[i]],
    #                         prec = 1 / lam2_save[it - 1] * Minv_list[[i]], log = TRUE, unnorm = FALSE), silent = TRUE)
    current_llik <- try(dmnorm(y = rep(template_feats,n), mu = stack(feat_list)$stack,
                               prec = 1 / lam2_save[it - 1] * Matrix::bdiag(Minv_list), log = TRUE, unnorm = FALSE), silent = TRUE)
    #  if(class(current_llik) == "try-error"){
    #    current_llik <- dmnorm(y = template_feats, mu = feat_list[[i]],
    #                               cov = lam2_save[it - 1] * M_list[[i]], log = TRUE, unnorm = FALSE)
    #  }
    if(class(current_llik) == "try-error"){
      current_llik <- dmnorm(y = rep(template_feats,n), mu = stack(feat_list)$stack,
                             cov = lam2_save[it - 1] * Matrix::bdiag(M_list), log = TRUE, unnorm = FALSE)
    }
    current_lprior <- dunif(alpha_save[it - 1], min = aa, max = ba, log = TRUE)

    #  #cand_alpha <- runif(1, aa, ba)
    cand_alpha <- rnorm(1, alpha_save[it - 1], tune)
    cand_lprior <- dunif(cand_alpha, min = aa, max = ba, log = TRUE)
    if(cand_lprior == -Inf){
      alpha_save[it] <- alpha_save[it - 1]
      #otime_save[[i]][it,] <- otime[[i]]
      #next
    } else{
      #  cand_M <- create_M(x_pts = feat_list[[i]], alpha = cand_alpha)
      cand_M_list <- lapply(1:n, function(i) create_M(x_pts = feat_list[[i]], alpha = cand_alpha))
      #  cand_Minv <- try(chol2inv(chol(cand_M)), silent = TRUE)
      cand_Minv_list <- try(lapply(1:n, function(i) chol2inv(chol(cand_M_list[[i]]))), silent = TRUE)
      if(class(cand_Minv_list) == "try-error"){
        alpha_save[it] <- alpha_save[it - 1]
        #otime_save[[i]][it,] <- otime[[i]]
        #next
      } else{
        #  cand_llik <- try(dmnorm(y = template_feats, mu = feat_list[[i]],
        #                      prec = 1 / lam2_save[it - 1] * cand_Minv, log = TRUE, unnorm = TRUE), silent = TRUE)
        cand_llik <- try(dmnorm(y = rep(template_feats,n), mu = stack(feat_list)$stack,
                                prec = 1 / lam2_save[it - 1] * Matrix::bdiag(cand_Minv_list), log = TRUE, unnorm = FALSE), silent = TRUE)
        #
        #  if(class(cand_llik) == "try-error"){
        #    cand_llik <- dmnorm(y = template_feats, mu = feat_list[[i]],
        #                        cov = lam2_save[it - 1] * cand_M, log = TRUE, unnorm = FALSE)
        #  }
        if(class(cand_llik) == "try-error"){
          cand_llik <- dmnorm(y = rep(template_feats,n), mu = stack(feat_list)$stack,
                              cov = lam2_save[it - 1] * Matrix::bdiag(cand_M_list), log = TRUE, unnorm = FALSE)
        }

        lratio <- cand_llik + cand_lprior - current_llik - current_lprior
        #otime_tmp <- as.numeric(time + create_M(x_pts = time, y_pts = feat_list[[i]], alpha = cand_alpha) %*%
        #                          cand_Minv%*%(template_feats - feat_list[[i]]))
        if(log(runif(1)) < lratio){ #&& is_monotone(otime_tmp, strict = TRUE)){
          accepts <- accepts + 1
          alpha_save[it] <- cand_alpha
          M_list <- cand_M_list
          Minv_list <- cand_Minv_list
          #otime_save[[i]][it,] <- otime[[i]] <- otime_tmp
        }
        else{
          alpha_save[it] <- alpha_save[it - 1]
          #otime_save[[i]][it,] <- otime[[i]]
        }
        #}
      }
    }
    #-- Update Lam2 --#
    #for(i in 1:n){#3
    #  lam2_save[[i]][it] <- update_normal_invgamma(y = template_feats, a = al, b = bl,
    #                                               mu = feat_list[[i]], R_inv = Minv_list[[i]])
    #}
    lam2_save[it] <- update_normal_invgamma(y = rep(template_feats, n), a = al, b = bl,
                                            mu = stack(feat_list)$stack, R_inv = Matrix::bdiag(Minv_list))

    #-- Update otime
    for(i in 1:n){
      otime[[i]] <- otime_save[[i]][it,] <- monotonize(as.numeric(time + create_M(x_pts = time, y_pts = feat_list[[i]], alpha = alpha_save[it]) %*%
                                                                    Minv_list[[i]]%*%(template_feats - feat_list[[i]])), forced = feat_inds[[i]])
      #chol2inv(chol(create_M(feat_list[[i]], alpha = alpha_save[[i]][it])))%*%(template_feats - feat_list[[i]]))
    }

    #-- Update wtime, eta2
    for(i in 1:n){
      curll <- sum(dnorm(wtime[[i]], otime[[i]], eta2_save[it-1,i],log = TRUE))
      curlp <- dunif(eta2_save[it-1,i],0,be,log = TRUE)
      cand_eta2 <- abs(runif(1,-be/10,be/10) + eta2_save[it-1,i])
      candlp <- dunif(cand_eta2,0,be, log = TRUE)
      candll <- sum(dnorm(wtime[[i]], otime[[i]], cand_eta2,log = TRUE))
      lrat <- candll + candlp - curll - curlp
      if(lrat > log(runif(1))){
        eta2_save[it,i] <- cand_eta2
      } else {
        eta2_save[it,i] <- eta2_save[it-1,i]
      }
      #eta2_save[it,i] <- 0.5#####also change default
      wtime[[i]] <- wtime_save[[i]][it,] <- monotonize(otime[[i]] + c(0,rnorm(m-2, 0, eta2_save[it,i]),0), forced = c(1,m))
    }


    #-- Update H_stack --#
    for(i in 1:n){
      H_list[[i]] <- bs(wtime[[i]], knots = knot_loc_p, intercept = T)
    }
    H_stack <- stack_Matrix(H_list)

    #-- Update Beta --#
    beta_save[it,] <- update_normal_normal(y = y_vec, X = H_stack, mu = mb,
                                           Sig_inv = 1 / sig2_save[it - 1] * diag_nm,
                                           V_inv = 1 / tau2_save[it - 1] * P)

    #-- Update Sig2 --#
    sig2_save[it] <- update_normal_invgamma(y = y_vec, a = asig, b = bsig,
                                            mu = H_stack %*% beta_save[it,],
                                            R_inv = diag_nm)

    #-- Update Tau2 --#
    tau2_save[it] <- update_normal_invgamma(y = beta_save[it,], a = at, b = bt,
                                            mu = mb, R_inv = P)

  }
  close(bar)
  accepts <- accepts / nrun
  wtime_post <- lapply(wtime_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  Hlist_post <- lapply(wtime_post, function(w) cbs(w, int_p))
  beta_post <- apply(beta_save[-c(1:nburn),], 2, mean)
  eta2_post <- apply(eta2_save[-c(1:nburn),], 2, mean)
  sig2_post <- mean(sig2_save[-c(1:nburn)])
  tau2_post <- mean(tau2_save[-c(1:nburn)])
  lam2_post <- lapply(lam2_save, function(lvec) mean(lvec[-c(1:nburn)]))
  alpha_post <- lapply(alpha_save, function(avec) mean(avec[-c(1:nburn)]))

  y_post <- lapply(1:n, function(i) as.numeric(Hlist_post[[i]] %*% beta_post))
  y_reg <- lapply(1:n, function(i) interp_spline(wtime_post[[i]], y_list[[i]]))
  mean_post <- as.numeric(Hp %*% beta_post)

  if(debug == TRUE) browser()
  return(list(y_post = y_post, y_reg = y_reg, mean_post = mean_post, accepts = accepts, wtime_post = wtime_post))
}

