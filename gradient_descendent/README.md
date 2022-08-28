# Thomsom's problem minimizer


## Description
This code implements a gradient-descendent routine to extract an homogeneous distribution of particles on a spherical shell. This provides solutions to the the so-called Thomson's problem. 

## Compilation
Simply call `make` on the main folder. This program requires the Eigen (version 3.4) library to be installed in your system. For the best performance, ensure OpenMP is available.

## Execution
To execute the program call the script `exec_thomson.sh N R` where `N` is the number of particles on the spherical shell and `R` is spherical radius. 

## Motivation
This was selected as a didactic short project to learn different topics at once:
1. gradient-based methods on surfaces.
2. Actively use fully object-oriented programming.
3. Get in touch and work with the Eigen library.
4. Learn the VTK file formating for ParaView visualisation.

