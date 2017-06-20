for ((dt = 3; dt <= 7; dt++ ))
do
	echo "-------"
	echo $dt
	echo "-------"
	sed -i 's/dt=.*/dt=1e-'"$dt"'/g' setting2.ini
	./app_modules_OpenMP_Euler setting2.ini function2-65.txt out_Euler_65/out_65_$dt.txt  4 >> log_euler65.txt
	python3 fault.py out_Euler_65/function2-129-m65.txt out_Euler_65/out_65_$dt.txt >> faults_euler.txt
done