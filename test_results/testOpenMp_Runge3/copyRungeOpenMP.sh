#!/bin/bash 
PATH_PROGRAM=$(find ../cmake-build-release -name 'app_modules_OpenMP_Runge')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .