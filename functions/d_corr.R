#' Convert event matrix to correlation matrix
#'
#' @param D Event matrix.
#'
#' @return Correlation matrix.
#' 

d_corr <- function(D){
  B <- matrix(0, nrow=nrow(D), ncol=nrow(D))
  diag(B) <- 1 / sqrt(diag(D))
  return(B %*% D %*% B)
}