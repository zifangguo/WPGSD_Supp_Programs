## Collect and output simulation results ##
## Only needed if simulation is run on the cluster ##

rm(list=ls())
library(kableExtra)
library(tidyverse)

setwd("/work/bards/guozi/GSD/WPGSD_Supp_Programs")

result <- NULL
for (i in 1:4){
  for (j in 1:3){
    for (task in 1:20) {
      filenam <- paste("./simulations/outtable/sim_result_HR", i, "_prop", j ,"_task", task, ".csv", sep="")
      tmp <- read_csv(filenam)
      result <-  bind_rows(result, tmp)
    }
  }
}

dim(result)

final <- result %>% 
              group_by(hra,hrb,hrab, hrneg, pa, pb, pab, pneg, method, N, event_ia, event_fa) %>%
              summarise(n_rep = sum(n_rep), 
                        power_H1 = mean(power_H1), power_H2 = mean(power_H2), 
                        power_H3 = mean(power_H3), power_Any = mean(power_Any), .groups='drop')%>%
              arrange(hrab, desc(hrneg), pneg)

final


## output to csv
write.csv(final, file = "./simulations/results/sim_all.csv", row.names = FALSE)

## output to latex format
output <- final %>% 
  select(hra,hrb,hrab, hrneg, pa, pb, pab, pneg, method, power_H1, power_H2, power_H3, power_Any) 

output %>% kable("latex", booktabs=TRUE, digits = 3, align = "c")



