# Thomsom's problem minimizer



## Description
This code implements a Monte Carlo minimisation routine with annealing to extract an homogeneous distribution of particles on a spherical shell. This provides solutions to the the so-called Thomson's problem. 


## Compilation
Simply call `make` on the main folder. For the best performance, ensure OpenMP is available in your system.


## Execution
To execute the program call the script `exec_mcsphere.sh N M` where `N` is the number of particles on the spherical shell and `M` is the number of Monte Carlo minimisation steps. The number of steps depends on the number of particles. For instance, a good starting point is `exec_mcsphere.sh 500 300000`. 
