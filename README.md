# Thomsom's problem minimizeres



## Description
This repository shows minimisation via Monte Carlo or gradient descendent minimisation to extract an homogeneous distribution of particles on a spherical shell. This provides solutions to the the so-called Thomson's problem. 


## Compilation
Simply call `make` on the main folder of each project.


## Execution
An interesting test is using a small number of particles to reproduce the first platonic solids (see examples folder). Ensure you have execution permissions for the script. If the user does not have execution rigths on the script, write the command: `chmod +x exec_mcsphere.sh` or `chmod +x exec_mcthomso.sh`.

The aforementioned bash script excecutes the main program, which produces the file `rsphere_eq.dat` with the coordinates of the particles on the `bin` folder. After minimisation, the Python script `sphere_triangulation.py` produces an optimal triangulation of the surface using the Convex hull. The triangulation is saved on the file `mesh.dat`. In addition, the Python script analyzes for the quality of the data and suggests if more Monte Carlo steps are necessary for an optimal minimization. The larger the number of particles, the larger the number of defects, but the topological charge should be *q* = 12 for every sphere. 
