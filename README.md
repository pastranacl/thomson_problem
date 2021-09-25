# Thomsom's problem minimizers 



## Description
These programs use two different strategies to extract an homogeneous distribution of particles on a spherical shell, providing solutions to the the so-called Thomson's problem. Two different minimisation algorithms:
1. Monte Carlo minimisation routine with simulated annealing
2. Gradient descendent using Barzilai-Borwain's update rule 


## Compilation
Simply call `make` on the main folder. For the best performance, ensure OpenMP is available in your system.


## Execution
Each minimisation routines has its own description of the execution procedure. As a first test, call the minimization of small number of particles. These allows to reproduce the first Platonic solids. Ensure you have execution permissions for exeuction scripts. If the user does not have execution rigths on the script, write the command: `chmod +x script_name.sh`.

The bash scripts excecute the main program and produce the file `rsphere_eq.dat` on the `bin` folder. This file has coordinates of the particles on the spherical shell. After minimisation, the Python script `sphere_triangulation.py` produces an optimal triangulation of the surface from the the Convex hull. The triangulation is saved on the file `mesh.dat`. In addition, the Python script analyzes for the quality of the distribution of veretexes. It indicates if More Monte Carlo steps are necessary for an optimal minimization or if the parameters of the gradient descendent have to be optimized. The larger the number of particles, the larger the number of defects, but the topological charge should be *q* = 12 for every sphere. 
