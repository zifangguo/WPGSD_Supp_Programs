
######## Call function of the WPGSD simulations under cluster #############

# call this function in the cluster with:
# cd /work/bards/guozi/GSD/WPGSD_Supp_Programs/simulations
# module load R/4.0.2
# R --vanilla --slave < wrapper.R


for (i in 1:4){
  for (j in 1:3){
    nrep <- 5000
    qsub.call <- paste('qsub -cwd -V -l virtual_free=16G -l mem_free=16G -l mem_reserve=16G -l h_vmem=16G -t 1-20 -b y')
    r.command <- paste('R --vanilla --slave --args', i, j, nrep, "'<'", 'simulations.R')
    qsub.call <- paste(qsub.call,r.command)
    system(qsub.call) 
  }
}

