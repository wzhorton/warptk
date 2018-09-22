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
#' @param asig,bsig,at,bt,al,bl,aa,ba hyperparameters for inverse gamma priors on model variances
#' @param progress logical; indicates whether or not to print percent progress.
#' @return a list of warped curves, the estimated mean curve, and a vector of MCMC
#'   acceptance rates. Curve list elements are vectors with length matching
#'   the elements of y_list.
#' @export

template_gp_warp <- function(y_list, feat_list, template_feats,
                             niter = 5000, nburn = 10000, int_p = 20,
                             asig = .1, bsig = .1, at = .1, bt = .1, al = .1, bl = .1,
                             aa = .01, ba = .01, progress = TRUE){

  #----- Fixed Values -----#

  n <- length(y_list)
  m <- unique(sapply(y_list,length))
  r <- 1 + 3
  if(length(m) != 1) stop("All y_list elements must have equal length")

  time <- seq(0, 1, len = m)
  y_vec <- unlist(y_list)

  p <- int_p + r
  knot_loc_p <- seq(time[1], time[m], len = int_p+2)[-c(1,int_p+2)]
  Hp <- bs(time, knots = knot_loc_p, intercept = T)

  P <- K1(p)
  P[1,1] <- 2

  mb <- rep(0,p)

  tune <- 5
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
  for(i in 1:n) alpha_save[[i]][1] <- 50


  wtime <- lapply(1:n, function(i) time)
  wtime_save <- lapply(1:n, function(i){
    out <- matrix(NA, nrow = nrun, ncol = m)
    out[1,] <- time
    return(out)
  })
  H_list <- lapply(1:n, function(i) Hp)
  H_stack <- stack(H_list)$stack

  R0 <- fields::rdist(time)
  M_list <- lapply(1:n, function(i) fields::Matern(R0, alpha = 50, smoothness = 5.1))

  mn_gp <- function(x) x

  #----- MCMC Loop -----#

  cat("Progress:  0 %")
  for(it in 2:nrun){
    if (progress == TRUE) {
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun < 0.1)
        cat("\b\b\b\b", round(it/nrun * 100), "%")
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun >= 0.1)
        cat("\b\b\b\b\b", round(it/nrun * 100), "%")
    }

    #-- Update wtime, M, and alpha --#
    for(i in 1:n){
      current_llik <- dmnorm(y = y_list[[i]], mu =  H_list[[i]] %*% beta_save[it - 1,],
                             prec = 1/ sig2_save[it - 1] * diag(m), log = TRUE)
      current_lprior <- dinvgamma(x = alpha_save[[i]][it - 1], shape = aa, scale = ba, log = TRUE)

      cand_alpha <- abs(rnorm(1, mean = alpha_save[[i]][it - 1], sd = tune))
      cv_gp <- function(d) lam2_save[[i]][it - 1] * fields::Matern(d, alpha = cand_alpha, smoothness = 5.1)

      cand_up <- update_gaussian_process(x = time[feat_list[[i]]], y = time[template_feats],
                                         time = time, mnfun = mn_gp, covfun = cv_gp)
      cand_wtime <- as.numeric(cand_up$up_mean)
      cand_M <- cand_up$up_var
      cand_H <- bs(cand_wtime, knots = knot_loc_p, intercept = TRUE)

      cand_lprior <- dinvgamma(x = cand_alpha, shape = aa, scale = ba, log = TRUE)
      cand_llik <- dmnorm(y = y_list[[i]], mu = cand_H %*% beta_save[it - 1,],
                          prec = 1 / sig2_save[it - 1] * diag(m), log = TRUE)

      if(log(runif(1)) < cand_llik + cand_lprior - current_llik - current_lprior){
        accepts[i] <- accepts[i] + 1/nrun
        alpha_save[[i]][it] <- cand_alpha
        M_list[[i]] <- cand_M
        H_list[[i]] <- cand_H
        wtime[[i]] <- cand_wtime

      }
      else{
        alpha_save[[i]][it] <- alpha_save[[i]][it - 1]
      }
      wtime_save[[i]][it,] <- wtime[[i]]
    }

    #-- Update H_stack --#
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

    #-- Update Lam2 --#
    for(i in 1:n){
      lam2_save[[i]][it] <- update_normal_invgamma(y = wtime[[i]], a = al, b = bl,
                                                   mu = time, R = M_list[[i]])
    }

  }
  #accepts <- accepts / nrun
  wtime_post <- lapply(wtime_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  beta_post <- apply(beta_save[-c(1:nburn),], 2, mean)
  sig2_post <- mean(sig2_save[-c(1:nburn)])
  tau2_post <- mean(tau2_save[-c(1:nburn)])
  #lam2_post <- lapply(lam2_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  #alpha_post <- mean(alpha_save[-c(1:nburn)])

  y_post <- lapply(1:n, function(i) interp_spline(x = wtime_post[[i]], y = y_list[[i]]))
  mean_post <- as.numeric(Hp %*% beta_post)

  return(list(y_post = y_post, mean_post = mean_post, accepts = accepts))
}
