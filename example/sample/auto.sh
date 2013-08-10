#!/bin/sh

TIME=450					# DECA EXECUTION TIME FOR ONE GALAXY, CHANGE IT IF THIS PERIOD IS NOT ENOUGH
decapath=/home/mosenkov/diser/DECA_0.1/deca.py	# CHANGE THE PATH TO THE DECA DIRECTORY ON YOUR COMPUTER


T1=$(date +"%s")


rm -f deca_out.txt
rm -f ini_out.txt
rm -f sextr_out.txt

N=$(cat deca_input.dat | wc -l)


for k in $( seq 1 $N ) ; do
	( python $decapath $k 0 ) & pid=$!
	( sleep $TIME && kill -HUP $pid ) 2>/dev/null & watcher=$!
	if wait $pid 2>/dev/null; then
	    pkill -HUP -P $watcher
	    wait $watcher
	else
	    echo "The program was stuck somwhere!"
	    sh cleaner.sh
	fi
	killall -9 galfit
	killall -9 python
	sh cleaner.sh
done

T2=$(date +"%s")
T=$((T2-T1))

Th=$((T/3600))
Tm=$(((T-Th*3600)/60))
Ts=$((T-Th*3600-Tm*60))
echo "TOTAL TIME: $Th h $Tm m $Ts s"
