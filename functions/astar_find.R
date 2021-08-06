#' Utility function for root-finding to compute crossing probabilities with the overall alpha spending approach
#' 
#' @param a cumulative overall alpha spending up to current analysis
#' @param alpha_prev alpha boundary at previous interim analyses using the WPGSD approach
#' @param astar total nominal alpha level at current analysis from the WPGSD approach 
#' @param w vector of alpha weights at current analysis
#' @param sig correlation matrix of previous and current analyses test statistics
#' @return Difference. Should be 0 with a and astar identified.
#' @examples
  

astar_find <- function(a, alpha_prev=NULL, astar, w, sig, ...){
   # Remove column name for proper pmvnorm run
   colnames(sig) <- NULL
   
   if (is.null(alpha_prev)) {
    res <- 1 - a - mvtnorm::pmvnorm(lower = -Inf, 
                             upper = qnorm(1 - w * astar), 
                             sigma = sig,
                             algorithm = mvtnorm::GenzBretz(maxpts=50000,abseps=0.00001))
    } else {
    res <-  1 - a - mvtnorm::pmvnorm(lower = -Inf, 
                              upper = c(qnorm(1 - alpha_prev), qnorm(1 - w * astar)), 
                              sigma = sig,
                              algorithm = mvtnorm::GenzBretz(maxpts=50000, abseps=0.00001))
  }
  return(res)
}
