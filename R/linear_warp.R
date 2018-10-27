#### linear_warp.R ####

#' Linear Warping
#'
#' Warps a list of curves by fitting a linear warping function with provided template landmarks.
#'
#' @param y_list list of curves with equal length
#' @param feat_list list of feature vectors with equal length
#' @param template_feats a feature vector with length mathing those in feat_list
#' @param niter,nburn MCMC iterations and burns. Total run time is the sum.
#' @param int_p number of internal knots specified for the data curves
#' @param asig,bsig,at,bt hyperparameters for inverse gamma priors on model variances
#' @param progress logical; indicates whether or not to print percent progress.
#' @param debug logical; stops function right before return and enters debug mode.
#' @return a list of warped curves and an estimated mean curve. List elements are
#'   vectors with length matching the elements of y_list.
#' @export

linear_warp <- function(y_list, feat_list, template_feats, niter = 1000, nburn = 1000,
                        int_p = 20, asig = .1, bsig = .1, at = .1, bt = .1,
                        progress = TRUE, debug = FALSE) {

  #----- Fixed Values -----#

  n <- length(y_list)
  m <- unique(sapply(y_list,length))
  r <- 3 + 1
  if(length(m) != 1) stop("All y_list elements must have equal length")

  time <- seq(0, 1, len = m)
  y_vec <- unlist(y_list)

  p <- int_p + r
  knot_loc_p <- seq(time[1], time[m], len = int_p+2)[-c(1,int_p+2)]
  Hp <- bs(time, knots = knot_loc_p, intercept = T)

  P <- K1(p)
  P[1,1] <- 2
  mb <- rep(0,p)


  #----- Save Structures -----#

  nrun <- nburn + niter

  beta_save <- matrix(NA, nrow = nrun, ncol = p)
  beta_save[1,] <- mb

  sig2_save <- numeric(nrun)
  sig2_save[1] <- 1
  tau2_save <- numeric(nrun)
  tau2_save[1] <- 1


  #----- Construct wtime and H -----#

  wtime <- lapply(1:n, function(i) lsfit(time[feat_list[[i]]],
                                         time[template_feats],
                                         intercept = FALSE)$coefficients * time)

  H_list <- lapply(wtime, function(w) bs(w, knots = knot_loc_p, intercept = TRUE))
  Hstack <- stack(H_list)$stack

  #----- MCMC Loop -----#

  cat("Progress:  0 %")
  for(it in 2:nrun){
    if (progress == TRUE) {
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun < 0.1)
        cat("\b\b\b\b", round(it/nrun * 100), "%")
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun >= 0.1)
        cat("\b\b\b\b\b", round(it/nrun * 100), "%")
    }

    #-- Update Beta --#
    beta_save[it,] <- update_normal_normal(y = y_vec, X = Hstack, mu = mb,
                                           Sig_inv = 1 / sig2_save[it - 1] * diag(n*m),
                                           V_inv = 1 / tau2_save[it - 1] * P)

    #-- Update Sig2 --#
    sig2_save[it] <- update_normal_invgamma(y = y_vec, a = asig, b = bsig,
                                            mu = Hstack %*% beta_save[it,],
                                            R_inv = diag(n*m))

    #-- Update Tau2 --#
    tau2_save[it] <- update_normal_invgamma(y = beta_save[it,], a = at, b = bt,
                                            mu = mb, R_inv = P)
  }
  beta_post <- apply(beta_save[-c(1:nburn),], 2, mean)
  sig2_post <- mean(sig2_save[-c(1:nburn)])
  tau2_post <- mean(tau2_save[-c(1:nburn)])

  y_post <- lapply(1:n, function(i) as.numeric(H_list[[i]] %*% beta_post))
  y_reg <- lapply(1:n, function(i) interp_spline(wtime[[i]], y_list[[i]]))
  mean_post <- as.numeric(Hp %*% beta_post)

  if(debug == TRUE) browser()
  return(list(y_post = y_post, mean_post = mean_post))
}
