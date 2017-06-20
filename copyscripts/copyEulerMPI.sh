#!/bin/bash 
PATH_PROGRAM=$(find ../cmake-build-release -name 'app_modules_MPI_Euler')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .