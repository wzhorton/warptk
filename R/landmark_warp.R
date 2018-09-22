#### landmark_warp.R ####

#' Landmark Warping
#'
#' Warps a list of curves by matching landmarks to provided template landmarks.
#'
#' @param y_list list of curves with equal length
#' @param feat_list list of feature vectors with equal length
#' @param template_feats a feature vector with length mathing those in feat_list
#' @return a list of warped curves. List elements are vectors with length matching
#'   the elements of y_list.
#' @export

landmark_warp <- function(y_list, feat_list, template_feats){
  m <- length(y_list[[1]])
  time <- seq(0, 1, len = m)
  y_post <- lapply(1:length(y_list), function(i){
    wtime <- approx(x = c(0, time[feat_list[[i]]], 1),
                    y = c(0, time[template_feats], 1), xout = time)$y
    return(interp_spline(x = wtime, y = y_list[[i]]))
  })
  mean_post <- as.numeric(apply(matrix(unlist(y_post), nrow = m), 1, mean))

  return(list(y_post = y_post, mean_post = mean_post))
}
