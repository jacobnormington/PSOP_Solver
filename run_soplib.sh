#!/bin/bash

# delete the files to output values to
start=1 					# set to 1 in order to start from the beginning
startInstance=""			# instance before the start
endInstance=""				# instance after the end
dataset=0					# which instances to consider: 1 to exclude very hard instances, 0 for all instances, -1 for only very hard instances
num_threads=32 				# check 8, 16, 32
if [[ $start == 1 ]]; then
	rm -f outfile.log
	rm -f outfile_raw.log
	rm -f outfile_joined.log
fi
echo $num_threads threads >> outfile.log
echo $num_threads threads >> outfile_raw.log

# run the tests
for str in $(ls soplib);
do
	if [[ $start == 1 ]]; then
		if [[ $str == $endInstance ]]; then
			echo "...script ended early"
			exit 1
		fi 

		if [[ 	$dataset == 0 ||
				($dataset -gt 0 && ($str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" 
					&& $str != "R.700.1000.15.sop")) ||
				($dataset -lt 0 && !($str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" 
					&& $str != "R.700.1000.15.sop"))
			]]; then # dont skip the instances that always timeout on 32 threads
			echo $str >> outfile.log
			output=$(./sop_solver soplib/$str 8 soplib_config.txt)
			echo "$str" >> outfile_raw.log
			echo "$output" >> outfile_raw.log
			
			if [[ $( echo "$output" | egrep -o "instance timed out" ) != "" ]]; then # test times out
				echo TIMED OUT >> outfile.log
			elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # time and cost not printed at the end
					echo NO TIMING DATA >> outfile.log
			else # output timing data to file
				echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile.log # add " | cut -d, -f2" to find only the time
			fi
		fi
	elif [[ $str == $startInstance ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
# cat outfile.log | paste -d " " - - >> outfile_joined.log
