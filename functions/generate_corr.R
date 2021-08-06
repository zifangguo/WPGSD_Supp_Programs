#' Generate correlation matrix based on event counts 
#' 
#' @param event event count of each hypothesis at each analysis, including event count of the intersection of hypotheses. 
#' It contains 4 columns: H1, H2, Analysis, Event. H1 needs to be listed as 1, 2, 3 etc as numbers
#' 
#' 
#' 

generate_corr <- function(event){
  
  elem <- event %>% filter(H1==H2)
  inter <- event %>% filter(H1!=H2)
  n_hypotheses <- max(as.numeric(elem$H1))
  n_analyses <- max(elem$Analysis)
  
  # Diagonal
  D <- diag(elem$Event)
  
  # Within hypothesis across analyses
  for (i in 1:n_hypotheses) {
    for (j in 2:n_analyses){
      count <- as.numeric(event %>% filter(H1==i & H2==i & Analysis==j-1) %>% select(Event))
      D[i, n_hypotheses*(j-1)+i] <- count
      D[n_hypotheses*(j-1)+i, i] <- count
    }
  }
  
  # Between hypotheses
  for (i in 1:n_hypotheses) {
    for (j in c(1:n_hypotheses)[-i]) {
      for (k in 1:n_analyses){
        count1 <- as.numeric(event %>% filter(((H1==i & H2==j)|(H1==j & H2==i)) & Analysis==k) 
                                  %>% select(Event))
        D[n_hypotheses*(k-1)+i, n_hypotheses*(k-1)+j] <- count1
          for (l in c(1:n_analyses)[-k]){
            count2 <- as.numeric(event %>% filter(((H1==i & H2==j)|(H1==j & H2==i)) & Analysis==min(k,l)) 
                                 %>% select(Event))
            D[n_hypotheses*(k-1)+i, n_hypotheses*(l-1)+j] <- count2
        }
      }
    }
  }
  
  corr_mat <- d_corr(D)
  
  col_names <- NULL
  for (k in 1:n_analyses) {
    for (i in 1:n_hypotheses){
      name_tmp <- paste("H",i,"_A",k, sep="")
      col_names <- c(col_names, name_tmp)
    }
  }
  
  colnames(corr_mat) <- col_names
  
  return(corr_mat)
}

