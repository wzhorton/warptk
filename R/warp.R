#### warp.R ####

#' Filler Testing Function
#' @useDynLib warptk
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#' @export

core <- function(){
  P <- K1(15)
  P[1,1] <- 2
  x <- seq(0,1,len = 100)
  Y <- cbind(sin(x*10)+rnorm(100,0,.1)+1, sin(x*9.8)+rnorm(100,0,.1)+2, sin(x*10.2)+rnorm(100,0,.1)+3,
             sin(x^1.1*10)+rnorm(100,0,.1)+4, sin(x^1.1*9.8)+rnorm(100,0,.1)+5, sin(x^1.1*10.2)+rnorm(100,0,.1)+6)
  #.spline_core_C(Y, 1:100, P, 1, 1, 10000, 1000)
  .two_step_warp_C(Y, cbind(x,x,x,x,x,x), P, 10000, 1000, a_eps=.1, b_eps=1, a_a =.1, b_a=1,
                   a_c=.1, b_c=1, a_tau=.1, b_tau=1)
}

warp_landmark <- function(y, lmk_indx, ref_lmk_indx, p=20, n_iter=10000, n_burn=1000,
                          a_eps=.1, b_eps=1, a_a =.1, b_a=1,
                          a_c=.1, b_c=1, a_tau=.1, b_tau=1){
  time <- seq(0,1,len=nrow(y))
  P <- K1(p)
  wtime <- sapply(1:nrow(y), function(i) approx(x = c(0, time[lmk_indx[,i]], 1),
                                          y = c(0, time[ref_lmk_indx], 1), xout = time)$y)
  .two_step_warp_C(y, wtime, P, niter, nburn,a_eps=.1, b_eps=1, a_a =.1, b_a=1,
                   a_c=.1, b_c=1, a_tau=.1, b_tau=1)
}

