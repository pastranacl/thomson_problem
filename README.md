# Thomsom's problem minimizer



## Description
This code implements a Monte Carlo minimisation routine with annealing to extract an homogeneous distribution of particles on a spherical shell. This provides solutions to the the so-called Thomson's problem. 


## Compilation
Simply call `make` on the main folder. For the best performance, ensure OpenMP is available in your system.


## Execution
To execute the program call the script `exec_mcsphere.sh N M` where `N` is the number of particles on the spherical shell and `M` is the number of Monte Carlo minimisation steps. The number of steps depends on the number of particles. For instance, a good starting point is `exec_mcsphere.sh 500 300000`. Other interesting test is using a small number of particles to reproduce the first platonic solids. Ensure you have execution permissions for the script. If the user does not have execution rigths on the script, write the command: `chmod +x exec_mcsphere.sh`.

The bash script excecutes the main program, which produces the file `rsphere_eq.dat` with the coordinates of the particles on the `bin` folder. After minimisation, the Python script `sphere_triangulation.py` produces an optimal triangulation of the surface using the Convex hull. The triangulation is saved on the file `mesh.dat`. In addition, the Python script analyzes for the quality of the data and suggests if more Monte Carlo steps are necessary for an optimal minimization. The larger the number of particles, the larger the number of defects, but the topological charge should be *q* = 12 for every sphere. 
