#!/bin/bash
# Cesar L. Pastrana, 2021

export NUMEXPR_NUM_THREADS=8

if [ "$#" -ne 2 ]; then
    echo "Not sufficient arguments. Call as exec_mcsphere N_VERTEXES N_MC_STEPS"
    exit 2
fi

cd ./bin
./mcsphere $1 $2

if [ -f "./rsphere_eq.dat" ]; then
    python ./sphere_triangulation.py
fi

cd ..


