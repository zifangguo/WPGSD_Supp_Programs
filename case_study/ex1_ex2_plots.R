
### Generate multiplicity plots in the manuscript ###

rm(list=ls())
 
library(gsDesign)
library(ggplot2)

setwd("/work/bards/guozi/GSD/WPGSD_Supp_Programs")

####################### Example 1 ####################

########### Multiplicity ##########
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

nameHypotheses <- c("H1: Population 1",
                    "H2: Population 2", 
                    "H3: Overall Population")
m <- matrix(c(0,0,1,
              0,0,1,
              0.5,0.5,0),nrow=3,byrow=TRUE)
alphaHypotheses <- c(0.3, 0.3, 0.4)

hplot <- hGraph(3,alphaHypotheses=alphaHypotheses,m=m,
                nameHypotheses=nameHypotheses, trhw=.2, trhh=.1, 
                digits=5, trdigits=3, size=5, halfWid=1, halfHgt=0.5,      
                offset=0.2 , trprop= 0.4,  
                fill=as.factor(c(2,3,1)),
                palette=cbPalette[1:3],
                wchar = "w") 
hplot

# jpeg("ex1_multiplicity.jpeg", units="in", width=9, height=5, res=300)
# hplot
# dev.off()

ggsave("./case_study/ex1_multiplicity.eps", units="in", width=9, height=5, dpi=300)

####################### Example 1 Holm ####################

########### Multiplicity ##########
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

nameHypotheses <- c("H1: Population 1",
                    "H2: Population 2", 
                    "H3: Overall Population")
m <- matrix(c(0,3/7, 4/7,
              3/7,0,4/7,
              1/2,1/2,0),nrow=3,byrow=TRUE)
alphaHypotheses <- c(0.3, 0.3, 0.4)

hplot <- hGraph(3,alphaHypotheses=alphaHypotheses,m=m,
                nameHypotheses=nameHypotheses, trhw=.2, trhh=.1, 
                digits=5, trdigits=3, size=5, halfWid=1, halfHgt=0.5,      
                offset=0.2 , trprop= 0.4,  
                fill=as.factor(c(2,3,1)),
                palette=cbPalette[1:3],
                wchar = "w") 
hplot

# jpeg("ex1_Holm_multiplicity.jpeg", units="in", width=9, height=5, res=300)
# hplot
# dev.off()

ggsave("./case_study/ex1_Holm_multiplicity.eps", units="in", width=9, height=5, dpi=300)


####################### Example 2 ####################

########### Multiplicity ##########
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
nameHypotheses <- c("H1: Experimental 1 vs Control",
                    "H2: Experimental 2 vs Control",
                    "H3: Experimental 3 vs Control")
m <- matrix(c(0,0.5,0.5,
              0.5,0,0.5,
              0.5,0.5,0),nrow=3,byrow=TRUE)
alphaHypotheses <- c(1/3, 1/3, 1/3)

hplot <- hGraph(3,alphaHypotheses=alphaHypotheses,m=m,
                nameHypotheses=nameHypotheses, trhw=.2, trhh=.1,
                digits=3, trdigits=4, size=5, halfWid=1.2, halfHgt=0.5,
                offset=0.2 , trprop= 0.35,
                fill=as.factor(c(2,3,1)),
                palette=cbPalette[1:3],
                wchar = "w")
hplot

# jpeg("ex2_multiplicity.jpeg", units="in", width=9, height=5, res=300)
# hplot
# dev.off()

ggsave("./case_study/ex2_multiplicity.eps", units="in", width=9, height=5, dpi=300)
