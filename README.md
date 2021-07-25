# Thomsom problem



# Description
This code implements a Monte Carlo minimisation routine to extract an homogeneous distribution of particles on a spherical shell. This is the so-called Thomson's problem. The code uses OpenMP for parallel execution.


# Compilation
Simply call `make` on the main folder. 


# Execution
To execute the program call the script `exec_mcsphere.sh N M` where `N` is the number of particles on the spherical shell and `M` is the number of Monte Carlo minimisation steps. The number of steps depends on the number of particles. For instance, a good starting point is `exec_mcsphere.sh 500 10000`. Ensure you have execution permissions for the script. If the user is not allow write: `chmod +x exec_mcsphere.sh`.

The bash script excecutes the main program, which produces the file `rsphere_eq.dat` with the coordinates of the particles on the `bin` folder. After minimisation, the Python script `sphere_triangulation.py` produces an optimal triangulation of the surface using the Convex hull. The triangulation is saved on the file `mesh.dat`. In addition, the Python script analyzes for the quality of the data and suggests if more Monte Carlo steps are necessary for an optimal minimization.
