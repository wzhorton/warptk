#### utility.R ####


#' Evenly Spaced Curve Interpolation using Splines
#'
#' Interpolates curve values using cubic interpolation splines.
#' The resulting evaluations correspond to an evenly spaced grid of input points.
#' Often used to time-normalize registered curves onto a common time domain.
#' Note the quality of interpolation decreases with inadequate density of defining points.
#' If ties are found among input times, linear interpolation is used instead of splines.
#'
#' @param y numeric vector giving curve values
#' @param x numeric vector giving input times. Defaults to evenly spaced values
#' @param nout number of points defining resultant interpolated curve. Defaults to length of y
#' @export

interp_spline <- function(y, x = NULL, nout = NULL) {
  if(is.null(nout))
    nout <- length(y)
  if(is.null(x))
    x <- seq(0,1,len = length(y))

  if(anyDuplicated(x))
    return(approx(x, y, n = nout)$y)

  return(spline(x, y, n = nout)$y)
}


#' First Order Penalty Matrix
#'
#' Generates a first order penalty matrix discussed by Lang and Brezger (2012).
#' In short, when used as a prior precision matrix for B-splines, this matrix imposes
#' a roughness penalty. Note that this matrix is not full rank.
#'
#' @param dim dimension of output matrix. Must be at least three.
#' @export

K1 <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(1,rep(2,dim-2),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}

#' Register Curves Given Warping Functions
#'
#' Computes registered curves from either given time warping vectors or an
#' MCMC array of time warping vectors though the inverse warping function.
#'
#' @param y numeric matrix; columnns should correspond to observed time normalized curves
#' @param wtime numeric matrix or three dimensional array; columns should correspond to warping function realizations and MCMC iterates should correspond to the third array dimension (if applicable).
#' @export

register <- function(y, wtime){
  if(!is.na(dim(wtime)[3])){
    wtime <- as.matrix(apply(wtime,3,mean))
  }
  sapply(1:ncol(y), function(i){
    interp_spline(wtime[,i], y[,i])
  })
}


