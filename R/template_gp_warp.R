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
                             aa = 1, ba = 30, tune = 2, progress = TRUE, debug = FALSE){

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
  Hp <- splines::bs(time, knots = knot_loc_p, intercept = T)

  P <- K1(p)
  P[1,1] <- 2

  mb <- rep(0,p)

  accepts <- numeric(n)

  nrun <- nburn + niter

  beta_save <- matrix(NA, nrow = nrun, ncol = p)
  beta_save[1,] <- mb

  sig2_save <- numeric(nrun)
  sig2_save[1] <- 1
  tau2_save <- numeric(nrun)
  tau2_save[1] <- 1
  lam2_save <- lapply(1:n, function(i) numeric(nrun))
  for(i in 1:n) lam2_save[[i]][1] <- 1
  alpha_save <- lapply(1:n, function(i) numeric(nrun))
  for(i in 1:n) alpha_save[[i]][1] <- (aa + ba)/2


  wtime <- lapply(1:n, function(i) time)
  wtime_save <- lapply(1:n, function(i){
    out <- matrix(NA, nrow = nrun, ncol = m)
    out[1,] <- time
    return(out)
  })
  H_list <- lapply(1:n, function(i) Hp)
  H_stack <- stack(H_list)$stack

  create_M <- function(x_pts, y_pts = x_pts, alpha){
    fields::Matern(fields::rdist(x_pts, y_pts), alpha = alpha, smoothness = 5.1)
  }

  M_list <- lapply(1:n, function(i) create_M(x_pts = feat_list[[i]], alpha = (aa + ba)/2))
  Minv_list <- lapply(M_list, solve)


  #----- MCMC Loop -----#

  cat("Progress:  0 %")
  for(it in 2:nrun){
    if (progress == TRUE) {
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun < 0.1)
        cat("\b\b\b\b", round(it/nrun * 100), "%")
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun >= 0.1)
        cat("\b\b\b\b\b", round(it/nrun * 100), "%")
    }

    #-- Update alpha, M, wtime --#
    for(i in 1:n){
      current_llik <- try(dmnorm(y = template_feats, mu = feat_list[[i]],
                             prec = 1 / lam2_save[[i]][it - 1] * Minv_list[[i]], log = TRUE, unnorm = FALSE), silent = TRUE)
      if(class(current_llik) == "try-error"){
        current_llik <- dmnorm(y = template_feats, mu = feat_list[[i]],
                                   cov = lam2_save[[i]][it - 1] * M_list[[i]], log = TRUE, unnorm = FALSE)
      }
      current_lprior <- dunif(alpha_save[[i]][it - 1], min = aa, max = ba, log = TRUE)

      #cand_alpha <- runif(1, aa, ba)
      cand_alpha <- rnorm(1, alpha_save[[i]][it - 1], tune)
      cand_lprior <- dunif(cand_alpha, min = aa, max = ba, log = TRUE)
      if(cand_lprior == -Inf){
        alpha_save[[i]][it] <- alpha_save[[i]][it - 1]
        #wtime_save[[i]][it,] <- wtime[[i]]
        next
      }
      cand_M <- create_M(x_pts = feat_list[[i]], alpha = cand_alpha)
      cand_Minv <- try(chol2inv(chol(cand_M)), silent = TRUE)
      if(class(cand_Minv) == "try-error"){
        alpha_save[[i]][it] <- alpha_save[[i]][it - 1]
        #wtime_save[[i]][it,] <- wtime[[i]]
        next
      }
      cand_llik <- try(dmnorm(y = template_feats, mu = feat_list[[i]],
                          prec = 1 / lam2_save[[i]][it - 1] * cand_Minv, log = TRUE, unnorm = TRUE), silent = TRUE)

      if(class(cand_llik) == "try-error"){
        cand_llik <- dmnorm(y = template_feats, mu = feat_list[[i]],
                            cov = lam2_save[[i]][it - 1] * cand_M, log = TRUE, unnorm = FALSE)
      }

      lratio <- cand_llik + cand_lprior - current_llik - current_lprior
      #wtime_tmp <- as.numeric(time + create_M(x_pts = time, y_pts = feat_list[[i]], alpha = cand_alpha) %*%
      #                          cand_Minv%*%(template_feats - feat_list[[i]]))
      if(log(runif(1)) < lratio){ #&& is_monotone(wtime_tmp, strict = TRUE)){
        accepts[i] <- accepts[i] + 1
        alpha_save[[i]][it] <- cand_alpha
        M_list[[i]] <- cand_M
        Minv_list[[i]] <- cand_Minv
        #wtime_save[[i]][it,] <- wtime[[i]] <- wtime_tmp
      }
      else{
        alpha_save[[i]][it] <- alpha_save[[i]][it - 1]
        #wtime_save[[i]][it,] <- wtime[[i]]
      }
    }

    #-- Update Lam2 --#
    for(i in 1:n){#3
      lam2_save[[i]][it] <- update_normal_invgamma(y = template_feats, a = al, b = bl,
                                                   mu = feat_list[[i]], R_inv = Minv_list[[i]])
    }

    #-- Update wtime
    for(i in 1:n){
      wtime[[i]] <- wtime_save[[i]][it,] <- monotonize(as.numeric(time + create_M(x_pts = time, y_pts = feat_list[[i]], alpha = alpha_save[[i]][it]) %*%
                                                         Minv_list[[i]]%*%(template_feats - feat_list[[i]])), forced = feat_inds[[i]])
                                                         #chol2inv(chol(create_M(feat_list[[i]], alpha = alpha_save[[i]][it])))%*%(template_feats - feat_list[[i]]))
    }

    #-- Update H_stack --#
    for(i in 1:n){
      H_list[[i]] <- splines::bs(wtime[[i]], knots = knot_loc_p, intercept = T)
    }
    H_stack <- stack(H_list)$stack

    #-- Update Beta --#
    beta_save[it,] <- update_normal_normal(y = y_vec, X = H_stack, mu = mb,
                                           Sig_inv = 1 / sig2_save[it - 1] * diag(n*m),
                                           V_inv = 1 / tau2_save[it - 1] * P)

    #-- Update Sig2 --#
    sig2_save[it] <- update_normal_invgamma(y = y_vec, a = asig, b = bsig,
                                            mu = H_stack %*% beta_save[it,],
                                            R_inv = diag(n*m))

    #-- Update Tau2 --#
    tau2_save[it] <- update_normal_invgamma(y = beta_save[it,], a = at, b = bt,
                                            mu = mb, R_inv = P)

  }
  accepts <- accepts / nrun
  wtime_post <- lapply(wtime_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  Hlist_post <- lapply(wtime_post, function(w) bs(w, knots = knot_loc_p, intercept = TRUE))
  beta_post <- apply(beta_save[-c(1:nburn),], 2, mean)
  sig2_post <- mean(sig2_save[-c(1:nburn)])
  tau2_post <- mean(tau2_save[-c(1:nburn)])
  lam2_post <- lapply(lam2_save, function(lvec) mean(lvec[-c(1:nburn)]))
  alpha_post <- lapply(alpha_save, function(avec) mean(avec[-c(1:nburn)]))

  y_post <- lapply(1:n, function(i) as.numeric(Hlist_post[[i]] %*% beta_post))
  y_reg <- lapply(1:n, function(i) interp_spline(wtime_post[[i]], y_list[[i]]))
  mean_post <- as.numeric(Hp %*% beta_post)

  if(debug == TRUE) browser()
  return(list(y_post = y_post, mean_post = mean_post, accepts = accepts))
}

