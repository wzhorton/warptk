#### warp.R ####

#' Bayesian Landmark Warping
#' @useDynLib warptk
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#' @export
#'
warp_landmark <- function(y, lmk_inds, ref_lmk_inds, p=10, n_iter=15000, n_burn=10000, n_thin = 1,
                          a_eps=length(y), b_eps=.1, a_a =.1, b_a=.1,
                          a_c=.1, b_c=.1, a_tau=.1, b_tau=.1){
  time <- seq(0,1,len=nrow(y))
  P <- K1(p)
  P[1,1] <- 2
  wtime <- sapply(1:ncol(y), function(i) approx(x = c(0, time[lmk_inds[,i]], 1),
                                          y = c(0, time[ref_lmk_inds], 1), xout = time)$y)
  mcmc_list <- .two_step_warp_C(y, wtime, P, n_iter, n_burn, n_thin, a_eps,
                                b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau)
  names(mcmc_list) <- c("a", "c", "beta", "sig2e", "sig2a", "sig2c", "tau2", "wtime", "mu")
  return(mcmc_list)
}

#' Bayesian Hierarchical Curve Registration
#'
#' q must be at least 6 (the min is 5 in reality, but my fast fromula assumes at least 6)
#'
#' @export

warp_bhcr <- function(y, p = 10, q = 6,  n_iter=15000, n_burn=10000, n_thin = 1,
                      a_eps=length(y), b_eps=.1, a_a =.1, b_a=.1, a_c = .1, b_c = .1,
                      a_tau=.1, b_tau=.1, a_lam = .1, b_lam = .1){
  time <- seq(0,1,len=nrow(y))
  knots_q <- c(0,0,0,seq(0,1,len = q-2),1,1,1)
  Upsilon <- numeric(q)
  Upsilon[1] <- 0
  for(i in 1:(q-1)){
    Upsilon[i+1] <- (knots_q[i+4] - knots_q[i+1])/(4-1) + Upsilon[i]
  }
  P <- K1(p)
  P[1,1] <- 2
  Q <- K1(q)
  Q[1,1] <- 2
  mcmc_list <- .bhcr_warp_C(y, time, P, Q, Upsilon, n_iter, n_burn, n_thin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam)
  names(mcmc_list) <- c("a", "c", "beta", "sig2e", "sig2a", "sig2c", "tau2", "lam2", "phi", "wtime", "mu")
  return(mcmc_list)
}

#' Template Prior Curve Registration
#' @export

warp_template <- function(y, lmk_inds = NULL, ref_lmk_inds = NULL, p = 10,
                      n_iter=15000, n_burn=10000, n_thin = 1,
                      a_eps=length(y)/2, b_eps=.1, a_a =.1, b_a=.1, a_c = .1, b_c = .1,
                      a_tau=.1, b_tau=.1, a_lam = .1, b_lam = .1){
  time <- seq(0,1,len=nrow(y))
  P <- K1(p)
  P[1,1] <- 2

  if(is.null(lmk_inds) && is.null(ref_lmk_inds)){
    ref_lmk <- 0:1
    feats <- matrix(0:1, nrow = 2, ncol = ncol(y))
    wtime_init <- sapply(1:ncol(y), function(i) time)
  } else {
    ref_lmk <- c(0, time[ref_lmk_inds], 1)
    feats <- apply(as.matrix(lmk_inds),2,function(cc) c(0,time[cc],1))
    wtime_init <- sapply(1:ncol(y), function(i) {
      spline(x = c(0, time[lmk_inds[,i]], 1), y = c(0, time[ref_lmk_inds], 1),
             xout = time, method = "hyman")$y})
  }
  mcmc_list <- .template_warp_C(y, time, wtime_init, feats, ref_lmk, P, n_iter, n_burn, n_thin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam)
  names(mcmc_list) <- c("a", "c", "beta", "sig2e", "sig2a", "sig2c", "tau2", "lam2", "alpha", "eta", "wtime", "mu")
  return(mcmc_list)
}

