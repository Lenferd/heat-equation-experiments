#!/bin/bash 
PATH_PROGRAM=$(find ../build -print0 | grep -FzZ 'app/app_modules_OpenMP_Euler')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .