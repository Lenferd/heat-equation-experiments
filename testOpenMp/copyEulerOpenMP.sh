#!/bin/bash 
PATH_PROGRAM=$(find ../build -print0 | grep -FzZ 'Release/app_modules_OpenMP_Euler.exe')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .