filename="faults_runge.txt"
for h in 9 17 33 65 129
do
	echo "*********"
	echo $h
	echo $h >> $filename
	echo "*********"
	name="out_Runge_"
	mkdir $name$h
	sed -i 's/NX=.*/NX='"$h"'/g' setting2.ini
	sed -i 's/NY=.*/NY='"$h"'/g' setting2.ini
	sed -i 's/NZ=.*/NZ='"$h"'/g' setting2.ini

	for (( dt = 2; dt <=7; dt++ ))
	do
		echo "-------"
		echo $dt
		echo "-------"
		sed -i 's/dt=.*/dt=1e-'"$dt"'/g' setting2.ini
		./app_modules_OpenMP_Runge setting2.ini function2-$h.txt $name$h/out_$h\_$dt.txt  16 >> $name$h/log_runge$h.txt
	done
done
