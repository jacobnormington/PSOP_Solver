#!/bin/bash

# delete the files to output values to
start=0
if [[ $start == 1 ]]; then
	rm -f outfile_tsp.log
	rm -f outfile_tsp_raw.log
	rm -f outfile_joined_tsp.log
fi

# run the tests
for str in $(ls tsplib);
do
	if [[ $start == 1 ]]; then
		if [[ $str != "ESC07_last.sop" && $str != "ESC11_last.sop" ]]; then #these are not "real" instances, ask Taspon about them
			if [[ $str != "ESC78.sop" && $str != "ft53.1.sop" && $str != "ft53.2.sop" && $str != "ft53.3.sop" && $str != "ft70.1.sop" && $str != "ft70.2.sop" && $str != "ft70.3.sop" 
				&& $str != "kro124p.1.sop" && $str != "kro124p.2.sop" && $str != "kro124p.3.sop" && $str != "kro124p.4.sop" && $str != "p43.1.sop" && $str != "p43.2.sop" 
				&& $str != "prob.100.sop" && $str != "rbg323a.sop" && $str != "rbg341a.sop" && $str != "rbg358a.sop" && $str != "rbg378a.sop" 
				&& $str != "ry48p.1.sop" && $str != "ry48p.2.sop" && $str != "ry48p.3.sop" ]]; then #skip the instances that always timeout on 32 threads
				echo $str >> outfile_tsp.log
				output=$(./sop_solver tsplib/$str 3600)
				echo "$str" >> outfile_tsp_raw.log
				echo "$output" >> outfile_tsp_raw.log
				exit_status=$?
				if [[ $exit_status -eq 124 ]]; then # test times out
					echo TIMED OUT >> outfile_tsp.log
				elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # test produces no usable timing data
					echo NO TIMING DATA >> outfile_tsp.log
				else # output timing data to file
					echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile_tsp.log # add " | cut -d, -f2" to find only the time
				fi

				if [[ $( echo $output | egrep -o "Total Progress = 100%" ) == "" ]]; then # if test doesn't end in 100% progress
					echo PROGRESS BAR NOT 100% >> outfile_tsp.log
					echo PROGRESS BAR NOT 100% >> outfile_tsp_raw.log
				fi
			fi
		fi
	elif [[ $str == "ft53.4.sop" ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
#cat outfile_tsp.log | paste -d " " - - >> outfile_joined_tsp.log
