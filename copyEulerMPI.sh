#!/bin/bash 
PATH_PROGRAM=$(find ../build -print0 | grep -FzZ 'Release/app_modules_MPI_Euler.exe')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .