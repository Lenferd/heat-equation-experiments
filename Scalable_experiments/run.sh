filename="faults_euler.txt"
for h in 129
do
	echo "*********"
	echo $h
	echo $h >> $filename
	echo "*********"
	name="out_Euler_"
	mkdir $name$h
	sed -i 's/NX=.*/NX='"$h"'/g' setting2.ini
	sed -i 's/NY=.*/NY='"$h"'/g' setting2.ini
	sed -i 's/NZ=.*/NZ='"$h"'/g' setting2.ini

	for (( dt = 1; dt <= 6; dt++ ))
	do
		echo "-------"
		echo $dt
		echo "-------"
		sed -i 's/dt=.*/dt=1e-'"$dt"'/g' setting2.ini
		./app_modules_OpenMP_Euler setting2.ini function2-$h.txt $name$h/out_$h\_$dt.txt  4 >> $name$h/log_euler$h.txt
		python3 fault.py function2-129m$h.txt $name$h/out_$h\_$dt.txt >> $filename
	done
done
