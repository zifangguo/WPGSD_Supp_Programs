#' Compute p value boundaries of the parametric MTP method with overall alpha spending for all hypotheses
#' 
#' @param type Boundary type. 
#'             0 = Bonferroni. Separate alpha spending for each hypotheses.
#'             1 = Fixed alpha spending for all hypotheses. Method 3a in the manuscript.
#'             2 = Overall alpha spending for all hypotheses. Method 3b in the manuscript.
#'             3 = Separate alpha spending for each hypotheses. Method 3c in the manuscript.
#'            
#' @param k Number of analyses up to the current analysis.
#' @param w Initial weights
#' @param m Transition matrix
#' @param corr Correlation matrix of all test statistics up to the current analysis. dim = k*length(w)
#' @param alpha Overall alpha
#' @param cum_alpha Cumulative alpha spent at each analysis. Only required for type=1.  
#' @param sf A list of alpha spending functions to spend alpha for each hypotheses. 
#'           If type = 0 or 3 then length equals to number of hypotheses.  
#'           If type = 1 then sf is not needed.
#'           If type = 2 then only the first component is used. 
#'          
#' @param sfparm A list of parameters to be supplied to sfs.  
#'           If type = 0 or 3 then length equals to number of hypotheses. 
#'           If type = 1 then sfparm is not needed.
#'           If type = 2 then only the first component is used. 
#'           
#' @param t A list of information fraction used for alpha spending, may be different from the actual information fraction. Each component corresponds to a hypothesis.
#'           If type = 0 or 3 then length equals to number of hypotheses.
#'           If type = 1 then t is not needed.
#'           If type = 2 then only the first component is used. 
#'             
#'        
#' @return A \code{tibble} with \code{k*(2^(n_hypotheses-1))} rows of p-value boundaries. 
#'         Inflation factor is also provided if type=3.
#' @examples

generate_bounds <- function(type=1, k=2, w=w, m=m, corr=corr, alpha=0.025, cum_alpha=NULL,
                            sf=sfHSD, sfparm=-4, t=c(0.5,1), ...){
  
  if (type==1 & is.null(cum_alpha)) {
    stop("Boundary type is 1 (fixed alpha spending) but no cummulative alpha was provided.")
  }
  
  if (type==2) {
    if (is.list(sf)) {sf <- sf[[1]]}  
    if (is.list(sfparm)) {sfparm <- sfparm[[1]]}
    if (is.list(t)) {t <- t[[1]]}  
  }
  
  # Number of hypotheses
  n_hypotheses <- length(w)
  
  # Get weights for all intersection hypotheses
  graph <- gMCP::matrix2graph(m)
  graph <- gMCP::setWeights(graph, w)
  
  # Set up hypothetical pvalues (0 or 1) to obtain all combinations
  pvals <- NULL
  for (i in 1:n_hypotheses) {
    if (i==1) {
      pvals <- data.frame(x=c(0,1))
      names(pvals) <- paste("pval_H",i,sep="")
    } else {
      tmp <- data.frame(x=c(0,1))
      names(tmp) <- paste("pval_H",i,sep="")
      pvals <- merge(pvals, tmp)
    }
  }
  
  # Weights for each intersection hypothesis
  inter_weight <- NULL
  for (i in 1:nrow(pvals)){
    pval_tmp <- as.numeric(pvals[i,])
    graph_tmp <- gMCP::gMCP(graph=graph, pvalues=pval_tmp, alpha=alpha)
    weight_tmp <- gMCP::getWeights(graph_tmp)
    inter_weight <- dplyr::bind_rows(inter_weight, weight_tmp)
  }
  
  inter_weight <- inter_weight[-1,]
  inter_weight <- replace(inter_weight, inter_weight == 0, NA)
  
  # Get boundaries
  bounds <- NULL
  for (j in 1:nrow(inter_weight)) {
    
    w_tmp0 <- inter_weight[j,]
    # Hypotheses included in the intersection hypothesis
    hypotheses <- col(w_tmp0)[!is.na(w_tmp0)]
    
    # Remove NA from weight
    w_tmp <-  w_tmp0[(!is.na(w_tmp0))]
    w_tmp0 <- as.numeric(w_tmp0)
    
    if (type==0) { # Bonferroni
      
      bounds_tmp <- tibble(Analysis = 1:k, 
                           Hypotheses = paste("H", hypotheses, sep="", collapse=", "))
      
      for (h in 1:n_hypotheses) {
        
        if (!h %in% hypotheses) {p_tmp <- NA} else {
          # Index to select from the correlation matrix
          indx <- expand.grid(h, (1:k))
          indx <- indx[,1]+(indx[,2]-1)*n_hypotheses
          corr_tmp <- corr[indx, indx]
          # Boundary for a single hypothesis across k for the intersection hypothesis
          p_tmp <- 1 - pnorm(gsDesign(k=k, test.type=1, 
                                      usTime=t[[h]], 
                                      n.I=corr_tmp[,ncol(corr_tmp)]^2,
                                      alpha=alpha*w_tmp0[h], 
                                      sfu=sf[[h]], sfupar=sfparm[[h]])$upper$bound)
        }
        # Record results
        h_var <- paste("H",h, sep="")
        bounds_tmp <-  bounds_tmp %>%
          mutate(!!h_var := p_tmp)
      }
      bounds <- bind_rows(bounds, bounds_tmp)
    } else { # WPGSD Methods
      for (i in 1:k){
        
        if (type %in% c(1,2)) {
          if (is.null(cum_alpha)) {
            alpha_tmp <- sf(alpha = alpha, t = t, param=sfparm)$spend[i]
          } else {
            alpha_tmp <- cum_alpha[i]
          }
          
          if (i==1) {alpha_prev=NULL}
          
          # index to select from the correlation matrix
          indx <- expand.grid(hypotheses, (1:i))
          indx <- indx[,1]+(indx[,2]-1)*n_hypotheses
          corr_tmp <- corr[indx, indx]
          
          p_tmp <-  w_tmp*uniroot(astar_find, a = alpha_tmp,
                                  alpha_prev = alpha_prev,
                                  w = w_tmp, sig = corr_tmp,
                                  lower = 0, upper = alpha_tmp*5,  
                                  tol = 1e-10)$root
        }
        
        if (type==3) {
          
          if (i==1) {alpha_prev=NULL}
          
          # First find Bonferroni spending
          cum_alpha_B <- NULL
          bounds_B <- NULL
          for (h in hypotheses) {
            
            indx_B <- expand.grid(h, (1:k))
            indx_B <- indx_B[,1]+(indx_B[,2]-1)*n_hypotheses
            corr_B_tmp <- corr[indx_B, indx_B]
            # Cummulative Bonferroni spending for a single hypothesis at anlaysis k
            cum_alpha_B_tmp <- sf[[h]](alpha = alpha*w_tmp0[h], t=t[[h]], param=sfparm[[h]])$spend[i]
            cum_alpha_B <- c(cum_alpha_B, cum_alpha_B_tmp)
            
            # Bonferroni nominal boundary for a single hypothesis at analysis i
            p_B_tmp <- 1 - pnorm(gsDesign(k=k, test.type=1, 
                                          usTime=t[[h]], 
                                          n.I=corr_B_tmp[,ncol(corr_B_tmp)]^2,
                                          alpha=alpha*w_tmp0[h], 
                                          sfu=sf[[h]], sfupar=sfparm[[h]])$upper$bound)[i]
            bounds_B <- c(bounds_B, p_B_tmp)
          }
          
          # Find inflation factor xi
          
          if (length(hypotheses)==1) {
            xi <- 1
          } else {
            # index to select from the correlation matrix
            indx <- expand.grid(hypotheses, (1:i))
            indx <- indx[,1]+(indx[,2]-1)*n_hypotheses
            corr_tmp <- corr[indx, indx]
            
            xi <- uniroot(xi_find, lower = 0.5, upper = 10,
                          a = sum(cum_alpha_B), 
                          alpha_prev = alpha_prev,
                          aprime = bounds_B,
                          sig = corr_tmp, 
                          tol=1e-10)$root           
          } 

          p_tmp <- xi*bounds_B
        }
        
        # record results
        pval_tmp <- rep(NA, n_hypotheses)
        pval_tmp[hypotheses] <- p_tmp 
        names(pval_tmp) <- paste("H", 1:n_hypotheses, sep="")
        
        if (type==3) {
          bounds_tmp <- tibble(Analysis = i, 
                               Hypotheses = paste("H", hypotheses, sep="", collapse=", "),
                               as.data.frame(t(pval_tmp)),
                               xi = xi)          
        } else {
          bounds_tmp <- tibble(Analysis = i, 
                               Hypotheses = paste("H", hypotheses, sep="", collapse=", "),
                               as.data.frame(t(pval_tmp)))          
        }

        bounds <- bind_rows(bounds, bounds_tmp)
        
        # Update alpha_prev
        alpha_prev <- c(alpha_prev, p_tmp)
      }
    }
  }
  bounds <- bounds %>% 
                  arrange(Analysis, Hypotheses, .by_group = FALSE)
  return(bounds)
}

