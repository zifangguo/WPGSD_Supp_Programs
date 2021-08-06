#' Utility function for root-finding to compute inflation factor xi with the separate alpha spending approach
#' 
#' @param a sum of cumulative alpha spending from the Bonferroni approach
#' @param alpha_prev alpha boundary at previous interim analyses using the MTP approach
#' @param aprime nominal alpha boundary from the Bonferroni approach
#' @param xi inflation factor
#' @param sig correlation matrix of previous and current analyses test statistics
#' @return Difference. Should be 0 with xi identified.
#' @examples


xi_find <- function(a, alpha_prev=NULL, aprime, xi, sig, ...){

  # Remove column name for proper pmvnorm run
  colnames(sig) <- NULL
  
  if (is.null(alpha_prev)) {
    res <- 1 - a - pmvnorm(lower = -Inf, 
                           upper = qnorm(1 - xi*aprime),
                           sigma = sig,
                           algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
  } else {
    res <-  1 - a - pmvnorm(lower = -Inf, 
                            upper =c(qnorm(1 - alpha_prev), qnorm(1 - xi*aprime)),
                            sigma = sig,
                            algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
  }
  return(res)
}
