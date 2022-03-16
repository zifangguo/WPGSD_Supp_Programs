
## Correlation matrix of the simulation in the Appendix Table A8-A10

rm(list=ls())
library(kableExtra)
library(tidyverse)

## Please change working directory before running
setwd("/work/bards/guozi/GSD/WPGSD_Supp_Programs")

# Source user-defined functions
my_path <- "./functions/"              
source_files <- list.files(my_path, "*.R$")   
sapply(paste0(my_path, source_files), source)


n32 <- 1000
n31 <- 500

p_mat <- matrix(c(0.2, 0.2, 0.5, 0.1,
                  0.2, 0.2, 0.4, 0.2,
                  0.3, 0.3, 0.1, 0.3),
                  byrow=TRUE, ncol=4)


for (i in 1:3) {
 
  p_vec <- p_mat[i,]
  
  n11 <- n31*(p_vec[1]+p_vec[3])
  n21 <- n31*(p_vec[2]+p_vec[3])
  n12_1 <- n31*p_vec[3]
  
  n12 <- n32*(p_vec[1]+p_vec[3])
  n22 <- n32*(p_vec[2]+p_vec[3])
  n12_2 <- n32*p_vec[3]
  
  event <- tribble(
    ~H1, ~ H2, ~Analysis, ~Event,
    1, 1, 1, n11,
    2, 2, 1, n21,
    3, 3, 1, n31,
    1, 2, 1, n12_1,
    1, 3, 1, n11,
    2, 3, 1, n21,
    1, 1, 2, n12,
    2, 2, 2, n22,
    3, 3, 2, n32,
    1, 2, 2, n12_2,
    1, 3, 2, n12,
    2, 3, 2, n22
  ) 
  cor_mat <- generate_corr(event)

  print(cor_mat %>% kable("latex", booktabs=TRUE, digits = 3, align = "c"))
}
