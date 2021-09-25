#!/bin/bash
# Cesar 2021


if [ "$#" -ne 2 ]; then
    echo "Not sufficient arguments. Call as exec_mcsphere N_VERTEXES N_MC_STEPS"
    exit 2
fi

# Determine number of CPUs to optimize parallel resources
nl=$(lscpu -p | wc -l)
NCPU=$(expr $nl - 4)
export  OMP_NUM_THREADS=$NCPU

# Execution
cd ./bin
echo "(Execution with" $NCPU "CPUs)"
./mcsphere $1 $2
echo

# Call riagulation
if [ -f "./rsphere_eq.dat" ]; then
    python ./sphere_triangulation.py
fi

cd ..


