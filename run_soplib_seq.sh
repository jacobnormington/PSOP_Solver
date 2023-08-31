#!/bin/bash

# delete the files to output values to
start=0 # set to 1 in order to start from the beginning
if [[ $start == 1 ]]; then
	rm -f outfile.log
	rm -f outfile_raw.log
	rm -f outfile_joined.log
fi

# run the tests
for str in $(ls soplib);
do
	if [[ $start == 1 ]]; then
		# if [[ $str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" && $str != "R.700.1000.15.sop" ]]; then
			echo $str >> outfile.log
			output=$(./sop_solver soplib/$str 3600)
			echo "$str" >> outfile_raw.log
			echo "$output" >> outfile_raw.log
			exit_status=$?
			if [[ $exit_status -eq 124 ]]; then # test times out
				echo TIMED OUT >> outfile.log
			elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # test produces no usable timing data
				echo NO TIMING DATA >> outfile.log
			else # output timing data to file
				echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile.log # add " | cut -d, -f2" to find only the time
			fi

			if [[ $( echo $output | egrep -o "Total Progress = 100%" ) == "" ]]; then # if test doesn't end in 100% progress
				echo PROGRESS BAR NOT 100% >> outfile.log
				echo PROGRESS BAR NOT 100% >> outfile_raw.log
			fi
		# fi
	elif [[ $str == "R.600.100.60.sop" ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
#cat outfile.log | paste -d " " - - >> outfile_joined.log
