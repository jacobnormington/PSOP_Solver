#!/bin/bash

# delete the files to output values to
start=1 					# set to 1 in order to start from the beginning
startInstance=""			# instance before the start
endInstance=""				# instance after the end
dataset=0					# which instances to consider: 1 to exclude very hard instances, 0 for all instances, -1 for only very hard instances
num_threads=32 				# check 8, 16, 32
if [[ $start == 1 ]]; then
	rm -f outfile_tsp.log
	rm -f outfile_tsp_raw.log
	rm -f outfile_joined_tsp.log

	echo $num_threads threads >> outfile_tsp.log
	echo $num_threads threads >> outfile_tsp_raw.log
fi

# run the tests
for str in $(ls tsplib);
do
	if [[ $start == 1 ]]; then
		if [[ $str == $endInstance ]]; then
			echo "...script ended early"
			exit 1
		fi 

		if [[ $str != "ESC07_last.sop" && $str != "ESC11_last.sop" ]]; then #these are not "real" instances, ask Taspon about them
			if [[ 	$dataset == 0 ||
					($dataset > 0 && ($str != "ESC78.sop" && $str != "ft53.1.sop" && $str != "ft53.2.sop" && $str != "ft53.3.sop" && $str != "ft70.1.sop" && $str != "ft70.2.sop" 
						&& $str != "ft70.3.sop" && $str != "kro124p.1.sop" && $str != "kro124p.2.sop" && $str != "kro124p.3.sop" && $str != "kro124p.4.sop" && $str != "p43.1.sop" 
						&& $str != "p43.2.sop" && $str != "prob.100.sop" && $str != "rbg323a.sop" && $str != "rbg341a.sop" && $str != "rbg358a.sop" && $str != "rbg378a.sop" 
						&& $str != "ry48p.1.sop" && $str != "ry48p.2.sop" && $str != "ry48p.3.sop")) ||
					($dataset < 0 && !($str != "ESC78.sop" && $str != "ft53.1.sop" && $str != "ft53.2.sop" && $str != "ft53.3.sop" && $str != "ft70.1.sop" && $str != "ft70.2.sop" 
						&& $str != "ft70.3.sop" && $str != "kro124p.1.sop" && $str != "kro124p.2.sop" && $str != "kro124p.3.sop" && $str != "kro124p.4.sop" && $str != "p43.1.sop" 
						&& $str != "p43.2.sop" && $str != "prob.100.sop" && $str != "rbg323a.sop" && $str != "rbg341a.sop" && $str != "rbg358a.sop" && $str != "rbg378a.sop" 
						&& $str != "ry48p.1.sop" && $str != "ry48p.2.sop" && $str != "ry48p.3.sop"))
				]]; then # dont skip the instances that always timeout on 32 threads
				echo $str >> outfile_tsp.log
				output=$(./sop_solver tsplib/$str $num_threads tsplib_config.txt)
				echo "$str" >> outfile_tsp_raw.log
				echo "$output" >> outfile_tsp_raw.log
				
				if [[ $( echo "$output" | egrep -o "instance timed out" ) != "" ]]; then # test times out
					echo TIMED OUT >> outfile_tsp.log
				elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # time and cost not printed at the end
						echo NO TIMING DATA >> outfile_tsp.log
				else # output timing data to file
					echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile_tsp.log # add " | cut -d, -f2" to find only the time
				fi
			fi
		fi
	elif [[ $str == $startInstance ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
# cat outfile_tsp.log | paste -d " " - - >> outfile_joined_tsp.log
