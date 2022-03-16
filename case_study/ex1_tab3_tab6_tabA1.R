# Example 1 results in Table 3, Table 6, and Table A1
rm(list=ls())
library(gsDesign)
library(gMCP)
library(tibble)
library(mvtnorm)
library(dplyr)

## Please change working directory before running
setwd("/work/bards/guozi/GSD/WPGSD_Supp_Programs")

# Source user-defined functions
my_path <- "./functions/"              
source_files <- list.files(my_path, "*.R$")   
sapply(paste0(my_path, source_files), source)

set.seed(1234)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Transition matrix in Figure 2
m <- matrix(c(0, 0, 1,
              0, 0, 1,
              1/2,1/2,0),nrow=3,byrow=TRUE)
# Initial weights
w <- c(0.3, 0.3, 0.4)

# Event count of intersection of paired hypotheses - Table 1
event <- tribble(
  ~H1, ~ H2, ~Analysis, ~Event,
  1, 1, 1, 100,
  2, 2, 1, 110,
  3, 3, 1, 225,
  1, 2, 1, 80,
  1, 3, 1, 100,
  2, 3, 1, 110,
  1, 1, 2, 200,
  2, 2, 2, 220,
  3, 3, 2, 450,
  1, 2, 2, 160,
  1, 3, 2, 200,
  2, 3, 2, 220
) 
event 

# Generate correlation from events
corr <- generate_corr(event)
corr # Correlation matrix in Table 3

# WPGSD bounds, spending method 3b
bound_WPGSD <- generate_bounds(type=2, k=2, w=w, m=m, corr=corr, alpha=0.025, 
                               sf=sfHSD, 
                               sfparm=-4, 
                               t=c(min(100/200,110/220,225/450),1))

# Bonferroni bounds
bound_Bonf <- generate_bounds(type=0, k=2, w=w, m=m, corr=corr, alpha=0.025, 
                              sf=list(sfHSD, sfHSD, sfHSD), 
                              sfparm=list(-4,-4,-4), 
                              t=list(c(0.5,1), c(0.5,1), c(0.5,1)))

# Combine and back-calculate xi
bounds <- left_join(bound_Bonf, bound_WPGSD, 
                    by=c("Hypotheses","Analysis"), 
                    suffix=c(".B", ".W")) 
bounds <- bounds %>% rowwise() %>% 
                     mutate(xi = sum(H1.W, H2.W, H3.W, na.rm=TRUE) /
                                 sum(H1.B, H2.B, H3.B, na.rm=TRUE))
# Reorder for output
bounds$order <- rep(c(5,2,1,3,6,4,7), 2)
bounds <- bounds %>% arrange(Analysis,order) 
# Table 6
bounds 

# Z-statistics boundary, Table A1
zbounds <- bounds %>% mutate(zH1.B = -qnorm(H1.B),
                             zH2.B = -qnorm(H2.B),
                             zH3.B = -qnorm(H3.B),
                             zH1.W = -qnorm(H1.W),
                             zH2.W = -qnorm(H2.W),
                             zH3.W = -qnorm(H3.W)) %>% 
                      select(Analysis, Hypotheses, zH1.B, zH2.B, zH3.B, zH1.W, zH2.W, zH3.W) 
zbounds

write.csv(corr, "./case_study/ex1_tab3.csv", row.names = FALSE)
write.csv(bounds, "./case_study/ex1_tab6.csv", row.names = FALSE)
write.csv(zbounds, "./case_study/ex1_tabA1.csv", row.names = FALSE)




