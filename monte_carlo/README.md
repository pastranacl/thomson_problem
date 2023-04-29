# Thomsom's problem minimizer



## Description
These codes implements a Monte Carlo minimisation routine with annealing to extract an homogeneous distribution of particles on a spherical shell. This provides solutions to the the so-called Thomson's problem. 

Two different implementations are provided, a C++ implementation and a Python implementation.


## C++
This implementation should work out-of-the-box in every Linux-based computer. It is recommended for small and medium number of vertices (12 to 2000).

### Compilation
Simply call `make` on the main folder. For the best performance, ensure OpenMP is available in your system.


### Execution
To execute the program call the script `exec_mcsphere.sh N M` where `N` is the number of particles on the spherical shell and `M` is the number of Monte Carlo minimisation steps. The number of steps depends on the number of particles. For instance, a good starting point is `exec_mcsphere.sh 500 300000`. 


## Python
The Python implementation employs JAX to carry-on the particle repulsion energies on the GPUs. This makes the execution much faster if the number of vertices is large. However, it is necessary to properly set-up the Python environment, having a Nvidia graphic card and the required drivers installed.

### Seting the environment
The code has been written using Python 3.9. A `requirements.txt` file is provided to install all the packages via `conda`. It is still convenient to resort for the manual installation of JAX to ensure that the GPU-compatible version of JAX has been installed. Read the JAX documentation [here](https://github.com/google/jax#installation).

### Execution
Set the number of particles and the number of files in the `mcsphere.py`. Then execute it with `python3 mcsphere`, ensuring that the environment is active.
