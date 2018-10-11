#### telesca_warp.R ####

#' Telesca Hierarchical Warping
#'
#' Implements a method of bayesian hierarchical curve registration described in a paper
#' by Telesca (date).
#'
#' @param y_list list of curves with equal length
#' @param niter,nburn MCMC iterations and burns. Total run time is the sum.
#' @param int_q number of internal knots specified for the warping function
#' @param int_p number of internal knots specified for the data curves
#' @param asig,bsig,at,bt,al,bl hyperparameters for inverse gamma priors on model variances
#' @param progress logical; indicates whether or not to print percent progress.
#' @return a list of warped curves, the estimated mean curve, and a vector of MCMC
#'   acceptance rates. Curve list elements are vectors with length matching
#'   the elements of y_list.
#' @export



telesca_warp <- function(y_list, niter = 5000, nburn = 10000, int_q = 5, int_p = 20,
                         asig = .1, bsig = .1, at = .1, bt = .1, al = .1, bl = .1,
                         progress = TRUE){

  #----- Fixed Values -----#

  n <- length(y_list)
  m <- unique(sapply(y_list,length))
  r <- 3 + 1
  if(length(m) != 1) stop("All y_list elements must have equal length")

  time <- seq(0, 1, len = m)
  y_vec <- unlist(y_list)

  q <-  int_q + r
  knot_loc_q <- seq(time[1], time[m], len = int_q+2)[-c(1,int_q+2)]
  Hq <- bs(time, knots = knot_loc_q, intercept = T)
  p <- int_p + r
  knot_loc_p <- seq(time[1], time[m], len = int_p+2)[-c(1,int_p+2)]
  Hp <- bs(time, knots = knot_loc_p, intercept = T)

  nu <- c(rep(time[1],r),knot_loc_q,rep(time[m],r))
  Upsilon <- (nu[r] - nu[1])/(r-1)
  for(i in 1:(int_q+r-1)){
    Upsilon[i+1] <- (nu[i+r] - nu[i+1])/(r-1) + Upsilon[i]
  }

  P <- K1(p)
  Q <- K1(q)
  P[1,1] <- 2
  Q[1,1] <- 2
  mb <- rep(0,p)

  tune <- .005
  accepts <- numeric(n*(q-2))


  #----- Save Structures -----#

  nrun <- nburn + niter

  beta_save <- matrix(NA, nrow = nrun, ncol = p)
  beta_save[1,] <- mb

  sig2_save <- numeric(nrun)
  sig2_save[1] <- 1
  tau2_save <- numeric(nrun)
  tau2_save[1] <- 1
  lam2_save <- numeric(nrun)
  lam2_save[1] <- 1

  phi <- lapply(1:n, function(i) Upsilon)
  wtime <- lapply(1:n, function(i) time)
  wtime_save <- lapply(1:n, function(i){
    out <- matrix(NA, nrow = nrun, ncol = m)
    out[1,] <- time
    return(out)
  })
  H_list <- lapply(1:n, function(i) Hp)
  H_stack <- stack(H_list)$stack


  #----- MCMC Loop -----#

  cat("Progress:  0 %")
  for(it in 2:nrun){
    if (progress == TRUE) {
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun < 0.1)
        cat("\b\b\b\b", round(it/nrun * 100), "%")
      if (((it%%round(nrun * 0.01)) == 0) && it/nrun >= 0.1)
        cat("\b\b\b\b\b", round(it/nrun * 100), "%")
    }

    #-- Update Phi --#
    for(i in 1:n){
      tmp_phi <- phi[[i]]
      current_llik <- dmnorm(y = y_list[[i]], mu =  H_list[[i]] %*% beta_save[it - 1,],
                             prec = 1 / sig2_save[it - 1] * diag(m), log = TRUE, unnorm = TRUE)
      current_lprior <- dmnorm(y = phi[[i]], mu = Upsilon,
                               prec = 1 / lam2_save[it - 1] * Q, log = TRUE, unnorm = TRUE)
      for(j in 2:(q-1)){
        tmp_phi[j] <- runif(1, min = max(tmp_phi[j] - tune, tmp_phi[j-1]),
                            max = min(tmp_phi[j] + tune, tmp_phi[j+1]))
        tmp_Hp <- bs(Hq %*% tmp_phi, knots = knot_loc_p, intercept = TRUE)
        cand_llik <- dmnorm(y = y_list[[i]], mu =  tmp_Hp %*% beta_save[it - 1,],
                            prec = 1 / sig2_save[it - 1] * diag(m), log = TRUE, unnorm = TRUE)
        cand_lprior <- dmnorm(y = tmp_phi, mu = Upsilon,
                              prec = 1 / lam2_save[it - 1] * Q, log = TRUE, unnorm = TRUE)
        lratio <- cand_llik + cand_lprior - current_llik - current_lprior

        if(log(runif(1)) < lratio){
          current_llik <- cand_llik
          current_lprior <- cand_lprior
          accepts[(i-1)*(q-2)+(j-1)] <- accepts[(i-1)*(q-2)+(j-1)] + 1
        }
        else{
          tmp_phi[j] <- phi[[i]][j]
        }
      }
      phi[[i]] <- tmp_phi
    }

    #-- Update wtime --#
    wtime <- lapply(1:n, function(i) Hq %*% phi[[i]])
    for(i in 1:n){
      wtime_save[[i]][it,] <- wtime[[i]]
    }

    #-- Update H_list --#
    H_list <- lapply(1:n, function(i) bs(wtime[[i]], knots = knot_loc_p, intercept = TRUE))

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
    lam2_save[it] <- update_normal_invgamma(y = unlist(phi), a = al, b = bl,
                                            mu = rep(Upsilon, n),
                                            R_inv = diag(n) %x% Q)

  }

  accepts <- accepts / nrun
  wtime_post <- lapply(wtime_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  Hlist_post <- lapply(wtime_post, function(w) bs(w, knots = knot_loc_p, intercept = TRUE))
  beta_post <- apply(beta_save[-c(1:nburn),], 2, mean)
  sig2_post <- mean(sig2_save[-c(1:nburn)])
  tau2_post <- mean(tau2_save[-c(1:nburn)])
  lam2_post <- mean(lam2_save[-c(1:nburn)])

  y_post <- lapply(1:n, function(i) as.numeric(Hlist_post[[i]] %*% beta_post))
  mean_post <- as.numeric(Hp %*% beta_post)

  return(list(y_post = y_post, mean_post = mean_post, accepts = accepts))
}


