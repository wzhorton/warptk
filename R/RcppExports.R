# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.bs_even_C <- function(time, nk) {
    .Call('_warptk_bs_even', PACKAGE = 'warptk', time, nk)
}

.det_sympd_C <- function(x, Log = FALSE) {
    .Call('_warptk_det_sympd', PACKAGE = 'warptk', x, Log)
}

.dist_C <- function(x, y) {
    .Call('_warptk_dist', PACKAGE = 'warptk', x, y)
}

.monotonize_C <- function(x, y) {
    invisible(.Call('_warptk_curve_monotonize', PACKAGE = 'warptk', x, y))
}

.two_step_warp_C <- function(ymat, wtime, P, niter, nburn, nthin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau) {
    .Call('_warptk_two_step_warp', PACKAGE = 'warptk', ymat, wtime, P, niter, nburn, nthin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau)
}

.bhcr_warp_C <- function(ymat, time, P, Q, U, niter, nburn, nthin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam) {
    .Call('_warptk_bhcr_warp', PACKAGE = 'warptk', ymat, time, P, Q, U, niter, nburn, nthin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam)
}

.template_warp_C <- function(ymat, time, wtime_init, lmk_time, ref_time, P, niter, nburn, nthin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam) {
    .Call('_warptk_template_warp', PACKAGE = 'warptk', ymat, time, wtime_init, lmk_time, ref_time, P, niter, nburn, nthin, a_eps, b_eps, a_a, b_a, a_c, b_c, a_tau, b_tau, a_lam, b_lam)
}

