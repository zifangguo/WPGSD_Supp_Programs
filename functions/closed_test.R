#' Perform closed testing procedure
#' 
#' @param bounds a tibble object of nominal p-value boundaries from generate_bounds() containing columns Analysis, Hypotheses, H1, H2 etc. 
#' @param p_obs a tibble object of observed p-values containing columns Analysis, H1, H2 etc.
#' @return  An outcome matrix summarizing the testing results.
#' @examples
#' 
#'   # Example p_obs as input to the closed_test() function
#'   p_obs <- bind_rows(tibble(Analysis = 1, H1 = 0.001, H2 = 0.001, H3 = 0.003),
#'                      tibble(Analysis = 2, H1 = 0.001, H2 = 0.001, H3 = 0.005))
#'                      
#'
#'

closed_test <- function(bounds, p_obs) {
   
   n_analyses <- max(p_obs$Analysis)
   n_hypotheses <- ncol(p_obs)-1
   
   result <- NULL
   
   for (i in 1:n_analyses) {
      
      # results comparing p-value with bound at current analysis
      p_tmp <- p_obs %>% filter(Analysis==i) %>% 
                         select(num_range("H", 1:n_hypotheses))
      bounds_tmp <- bounds %>% filter(Analysis==i) %>% 
                               select(num_range("H", 1:n_hypotheses))
      test_raw <- c(unlist(p_tmp)) < t(bounds_tmp)
      
      # number of intersection hypothesis
      n_inter <- ncol(test_raw)
      
      # initial testing result of each intersection hypothesis
      test_inter <- apply(test_raw, 2, any, na.rm=TRUE)
      
      # if a hypothesis was rejected in a previous analysis, then all 
      # intersection hypothesis including that hypothesis is rejected
      if (i!=1) {
         # previous testing results
         prev_res <- apply(result %>% select(num_range("H", 1:n_hypotheses)), 2, any)
         # hypothesis number that was rejected in any previous analyses
         prev_reject <- c(1:n_hypotheses)[prev_res]
         # intersection hypothesis that includes previous rejected hypothesis
         inter_reject <- matrix(!is.na(test_raw[prev_reject,]), ncol=n_inter)
         indx_inter_reject <- c(1:n_inter)[apply(inter_reject,2,sum)>0]
         # convert testing result to TRUE for above intersection hypothesis
         test_inter[indx_inter_reject] <- TRUE
      }
      
      # testing result of each elementary hypothesis
      test_tmp <- rep(NA, n_hypotheses)
      for (j in 1:n_hypotheses){
         indx <- !is.na(test_raw[j,])
         test_elem <- all(test_inter[indx])
         test_tmp[j] <- test_elem
      }
      names(test_tmp) <- paste("H", 1:n_hypotheses, sep="")
      test_tmp <- data.frame(t(test_tmp))
      test_tmp$Analysis <- paste("Analysis", i)
      result <- dplyr::bind_rows(result, test_tmp)
   }

   result[result==TRUE] <- "Success"
   result[result==FALSE] <- "Fail"
   rownames(result) <- NULL
   return(result)
}

