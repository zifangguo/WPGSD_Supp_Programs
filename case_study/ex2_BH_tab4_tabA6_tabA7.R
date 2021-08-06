# Example 2 BH weighting results in Table 4, Table A6 and A7
rm(list=ls())
library(gsDesign)
library(gMCP)
library(tibble)
library(mvtnorm)
library(dplyr)

setwd("/work/bards/guozi/GSD/WPGSD_Supp_Programs")

# Source user-defined functions
my_path <- "./functions/"              
source_files <- list.files(my_path, "*.R$")   
sapply(paste0(my_path, source_files), source)

set.seed(1234)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex2 BH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Transition matrix in Figure A2
m <- matrix(c(0, 0.5, 0.5,
              0.5, 0, 0.5,
              0.5,0.5,0),nrow=3,byrow=TRUE)
# Initial weights
w <- c(1/3, 1/3, 1/3)

# Event count of intersection of paired hypotheses - Table 2
event <- tribble(
  ~H1, ~ H2, ~Analysis, ~Event,
  1, 1, 1, 155,
  2, 2, 1, 160,
  3, 3, 1, 165,
  1, 2, 1, 85,
  1, 3, 1, 85,
  2, 3, 1, 85,
  1, 1, 2, 305,
  2, 2, 2, 320,
  3, 3, 2, 335,
  1, 2, 2, 170,
  1, 3, 2, 170,
  2, 3, 2, 170
) 
event  

# Generate correlation from events
corr <- generate_corr(event)
corr # correlation matrix in Table 4

# WPGSD bounds, spending method 3c
bound_WPGSD <- generate_bounds(type=3, k=2, w=w, m=m, corr=corr, alpha=0.025, 
                               sf=list(sfLDOF, sfLDOF, sfLDOF), 
                               sfparm=list(0,0,0), 
                               t=list(c(155/305,1), c(160/320,1), c(165/335,1)))

# Bonferroni bounds
bound_Bonf <- generate_bounds(type=0, k=2, w=w, m=m, corr=corr, alpha=0.025, 
                              sf=list(sfLDOF, sfLDOF, sfLDOF), 
                              sfparm=list(0,0,0), 
                              t=list(c(155/305,1), c(160/320,1), c(165/335,1)))

bounds <- left_join(bound_Bonf, bound_WPGSD, 
                    by=c("Hypotheses","Analysis"), 
                    suffix=c(".B", ".W")) 

# Reorder for output
bounds$order <- rep(c(5,2,1,3,6,4,7), 2)
bounds <- bounds %>% arrange(Analysis,order) 
# Table A6
bounds 

# Z-statistics boundary, Table A7
zbounds <- bounds %>% mutate(zH1.B = -qnorm(H1.B),
                             zH2.B = -qnorm(H2.B),
                             zH3.B = -qnorm(H3.B),
                             zH1.W = -qnorm(H1.W),
                             zH2.W = -qnorm(H2.W),
                             zH3.W = -qnorm(H3.W)) %>% 
                      select(Analysis, Hypotheses, zH1.B, zH2.B, zH3.B, zH1.W, zH2.W, zH3.W) 
zbounds

write.csv(corr, "./case_study/ex2_BH_tab4.csv", row.names = FALSE)
write.csv(bounds, "./case_study/ex2_BH_tabA6.csv", row.names = FALSE)
write.csv(zbounds, "./case_study/ex2_BH_tabA7.csv", row.names = FALSE)


