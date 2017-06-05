#!/bin/bash 
PATH_PROGRAM=$(find ../../build -name 'app_Sergey_OpenMP_Euler.exe')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .