# WPGSD_Supp_Programs README

Source code for manuscript entitled "A unified framework for weighted parametric group sequential design (WPGSD)" by Keaven M. Anderson, Zifang Guo, Jing Zhao and Linda Z. Sun

Please contact Zifang Guo (zifang.guo@merck.com) for questions, comments and remarks on the code (and to report bugs). 

The code was written under the following environment:

R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Packages versions:
gsdmvn_0.1.1        gsDesign2_0.1       simtrial_0.1.7.9003 mvtnorm_1.1-1      
testthat_2.3.2      tidyr_1.1.0         gsDesign_3.1.1      ggplot2_3.3.2      
tibble_3.0.3        gMCP_0.8-15         survival_3.1-12     dplyr_1.0.0           


Below is a brief summary of the programs in each subfolder:

1. functions: contains functions used in the computations. 
        - generate_corr(): generates correlation matrix from event counts. Calls an utility function d_corr(). 
        - generate_bounds(): computes nominal p-value boundaries using the weighted Bonferroni or the WPGSD approach. All 3 alpha-spending methods for WPGSD described in the manuscript are implemented. Calls utility functions astar_find() and xi_find().
        - closed_test(): performs the closed testing procedure with nominal p-value boundaries obtained from generate_bounds() and the observed p-values.
        
2. case_study: contains source codes to reproduce all examples in the manuscript.

3. simulations: contains simulation codes to reproduce simulation results in the manuscript. The simulation takes a very long time to run if running on a single PC (>1000 hours). As such, High Performing Computing cluster was used in the simulations and codes provided include the main simulation code (simulations.R), the wrapper file to call the simulation code (wrapper.R), and the code to compile results (compile.R). One can also modify the main simulation code to use loops (currently commented out) to run the simulations on a single PC with reduced replications.


