for h in 129
do
	filename="faults_euler_$h.txt"
	echo "*********"
	echo $h
	echo "*********"
	name="out_Euler_"

	for (( dt = 1; dt <= 7; dt++ ))
	do
		echo "-------"
		echo $dt
		echo "-------"
		sed -i 's/dt=.*/dt=1e-'"$dt"'/g' setting2.ini
		python3 fault.py $name$h/out_$h\_7.txt $name$h/out_$h\_$dt.txt >> faults/$filename
	done
done
