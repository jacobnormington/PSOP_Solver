#!/bin/bash

# delete the files to output values to
start=1				# set to 1 in order to start from the beginning
starting_instance="R.300.1000.15.sop"
dataset=0					# which instances to consider: 1 to exclude very hard instances, 0 for all instances, -1 for only very hard instances
num_threads=32 		# check 8, 16, 32
enable=1			# 1 for enabling LKH thread, 0 for disabled
if [[ $enable == 1 ]]; then
	config="soplib_config.txt"
else
	config="soplib_config_disable.txt"
fi

if [[ $start == 1 ]]; then
	rm -f outfile.log
	rm -f outfile_raw.log
	rm -f outfile_joined.log
	if [[ $enable == 1 ]]; then
		echo "lkh + B&B" >> outfile.log
		echo "lkh + B&B" >> outfile_raw.log
	else
		echo "B&B alone" >> outfile.log
		echo "B&B alone" >> outfile_raw.log
	fi
	echo $num_threads threads >> outfile.log
	echo $num_threads threads >> outfile_raw.log
fi

# run the tests
for str in $(ls soplib);
do
	if [[ $start == 1 ]]; then
		if [[ ${str: -6:2} != "BB" ]]; then #ignore the old B&B files that don't have necessary header data
			if [[ 	$dataset == 0 ||
					($dataset -gt 0 && ($str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" 
						&& $str != "R.700.1000.15.sop")) ||
					($dataset -lt 0 && !($str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" 
						&& $str != "R.700.1000.15.sop"))
				]]; then
				echo $str >> outfile.log
				output=$(./sop_solver soplib/$str $num_threads $config)
				echo "$str" >> outfile_raw.log
				echo "$output" >> outfile_raw.log

				if [[ $( echo "$output" | egrep -o "instance timed out" ) != "" ]]; then # test times out
					echo TIMED OUT >> outfile.log
				elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # time and cost not printed at the end
					echo NO TIMING DATA >> outfile.log
				else # output timing data to file
					echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile.log # add " | cut -d, -f2" to find only the time
				fi

				if [[ $( echo $output | egrep -o "Total Progress = 100%" ) == "" ]]; then # if test doesn't end in 100% progress
					echo PROGRESS BAR NOT 100% >> outfile.log
					echo PROGRESS BAR NOT 100% >> outfile_raw.log
				fi
			fi
		fi
	elif [[ $str == $starting_instance ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
# cat outfile.log | paste -d " " - - >> outfile_joined.log
