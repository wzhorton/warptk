#### utility.R ####


#' Spline Curve Interpolation
#'
#' Interpolates out equally spaced values by fitting a cubic spline on given points and
#' matching the endpoints. Often used to redefine a curve with warped time axis on an
#' even time basis. Note that quality of interpolation decreases with inadequate density
#' of defining points. If ties are found, linear interpolation is used.
#'
#' @param x,y coordinate vectors on which to fit a cubic spline.
#' @param nout number of points to interpolate out
#' @importFrom splines bs
#' @export

interp_spline <- function(x, y, nout = length(y)) {
  nout <- force(nout)
  runlen <- rle(x)
  reps <- which(runlen$lengths != 1)
  if(length(reps) > 0){
    #runpos <- cumsum(runlen$lengths) - runlen$lengths + 1
    return(approx(x, y, n = nout)$y)
    #for(i in reps){
    #  run_inds <- runpos[i]:(runpos[i]+runlen$lengths[i] - 1)
    #  y[runpos[i]] <- mean(y[run_inds])
    #  x[run_inds[-1]] <- NA
    #  y[run_inds[-1]] <- NA
    #}
    #x <- na.omit(x)
    #y <- na.omit(y)
  }
  return(spline(x, y, n = nout)$y)
}


#' Second Order Penalty Matrix
#'
#' Creates a second order penalty matrix used in fitting penalized b-splines.
#' Courtesy of CITATION NEEDED. Note that this function produces singular matrices.
#'
#' @param dim dimension of output matrix. Must be at least four.
#' @export

K2 <- function(dim) {
  tmp <- list(c(-2, rep(-4, dim - 3), -2), rep(1, dim))
  K <- Matrix::bandSparse(dim, k = -c(1:2), diag = tmp, symmetric = TRUE)
  diag(K) <- c(1, 5, rep(6, dim - 4), 5, 1)
  K <- matrix(as.numeric(K), nrow = dim, byrow = TRUE)
  K
}


#' First Order Penalty Matrix
#'
#' Creates a first order penalty matrix used in fitting penalized b-splines.
#' Courtesy of CITATION NEEDED. Note that this function produces singular matrices.
#'
#' @param dim dimension of output matrix. Must be at least three.
#' @export

K1 <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(1,rep(2,dim-2),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}
