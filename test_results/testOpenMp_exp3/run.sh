#!/bin/bash

for ((h = 1; h<=17; h+=1))
do 
    echo "*******"
    echo $h
    echo "*******"
    srun ./app_modules_OpenMP_Runge setting2.ini function2-64.txt out.txt $h
done