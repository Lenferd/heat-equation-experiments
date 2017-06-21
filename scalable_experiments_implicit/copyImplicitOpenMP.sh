#!/bin/bash 
PATH_PROGRAM=$(find ../build -name 'app_modules_OpenMP_Implicit')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .