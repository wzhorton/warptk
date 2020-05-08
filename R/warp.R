#### warp.R ####

#' Filler Testing Function
#' @export

core <- function(){
  Y <- as.matrix(read.csv("~/Downloads/fdat.csv", header = FALSE))[,1:18]
  warp_bhcr(Y)
}


#' Bayesian Landmark Warping
#' @useDynLib warptk
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#' @export
#'
warp_landmark <- function(y, lmk_indx, ref_lmk_indx, p=15, n_iter=10000, n_burn=5000, n_thin = 1,
                          a_eps=.1, b_eps=.1, a_a =.1, b_a=.1,
                          a_c=.1, b_c=.1, a_tau=.1, b_tau=.1){
  time <- seq(0,1,len=nrow(y))
  P <- K1(p)
  P[1,1] <- 2
  wtime <- sapply(1:nrow(y), function(i) approx(x = c(0, time[lmk_indx[,i]], 1),
                                          y = c(0, time[ref_lmk_indx], 1), xout = time)$y)
  .two_step_warp_C(y, wtime, P, n_iter, n_burn, n_thin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau)
}

#' Bayesian Hierarchical Curve Registration
#' @export

warp_bhcr <- function(y, p = 15, q = 5,  n_iter=15000, n_burn=50000, n_thin = 1,
                      a_eps=.1, b_eps=.1, a_a =.1, b_a=.1, a_c = .1, b_c = .1,
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
  .bhcr_warp_C(y, time, P, Q, Upsilon, n_iter, n_burn, n_thin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam)
}
