#!/bin/bash 
PATH_PROGRAM=$(find ../build -name 'app_modules_OpenMP_Runge')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .