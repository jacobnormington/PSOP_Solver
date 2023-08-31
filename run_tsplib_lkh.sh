#!/bin/bash

# delete the files to output values to
start=0 					# set to 1 in order to start from the beginning
startInstance="ESC78.sop"	# instance before the start
endInstance=""				# instance after the end
dataset=0					# which instances to consider: 1 to exclude very hard instances, 0 for all instances, -1 for only very hard instances
num_threads=32 				# check 8, 16, 32
enable=1					# 1 for enabling LKH thread, 0 for disabled
if [[ $enable == 1 ]]; then
	config="tsplib_config.txt"
else
	config="tsplib_config_disable.txt"
fi

if [[ $start == 1 ]]; then
	rm -f outfile_tsp.log
	rm -f outfile_tsp_raw.log
	rm -f outfile_joined_tsp.log
	if [[ $enable == 1 ]]; then
		echo "lkh + B&B" >> outfile_tsp.log
		echo "lkh + B&B" >> outfile_tsp_raw.log
	else
		echo "B&B alone" >> outfile_tsp.log
		echo "B&B alone" >> outfile_tsp_raw.log
	fi
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

		if [[ ${str: -6:2} != "BB" ]]; then #ignore the old B&B files that don't have necessary header data
			if [[ $str != "ESC07_last.sop" && $str != "ESC11_last.sop" ]]; then #these are not "real" instances, ask Taspon about them
				if [[ ($str != "ESC14.sop" && $str != "ESC98.sop" && $str != "prob.1.sop" && $str != "prob.2.sop" && $str != "prob.3.sop" && $str != "prob.4.sop" && $str != "prob.5.sop" 
					&& $str != "prob.6.sop" && $str != "prob.7.0.sop" && $str != "prob.7.30.sop" && $str != "prob.7.35.sop" && $str != "prob.7.40.sop" && $str != "prob.7.45.sop" 
					&& $str != "prob.7.50.sop" && $str != "prob.7.55.sop" && $str != "prob.7.60.sop" && $str != "prob.7.65.sop" && $str != "prob.7.70.sop" 
					&& $str != "rbg019a.sop" && $str != "rbg019b.sop" && $str != "rbg021a.sop" && $str != "rbg023a.sop" && $str != "rbg029a.sop" && $str != "rbg049a.sop" 
					&& $str != "rbg050a.sop" && $str != "rbg050b.sop" && $str != "rbg068a.sop" && $str != "rbg088a.sop" && $str != "rbg092a.sop" && $str != "rbg094a.sop"
					&& $str != "rbg105a.sop" && $str != "rbg113a.sop" && $str != "rbg117a.sop" && $str != "rbg118a.sop" && $str != "rbg124a.sop" && $str != "rbg126a.sop" 
					&& $str != "rbg143a.sop" && $str != "rbg148a.sop" && $str != "rbg161a.sop" && $str != "rbg190a.sop" && $str != "rbg219a.sop" && $str != "rbg247a.sop"
					&& $str != "rbg285a.sop") ]]; then # skip instances that only have LKH-structured data files
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
						output=$(./sop_solver tsplib/$str $num_threads $config)
						echo "$str" >> outfile_tsp_raw.log
						echo "$output" >> outfile_tsp_raw.log
						
						if [[ $( echo "$output" | egrep -o "instance timed out" ) != "" ]]; then # test times out
							echo TIMED OUT >> outfile_tsp.log
						elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # time and cost not printed at the end
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
			fi
		fi
	elif [[ $str == $startInstance ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
# cat outfile_tsp.log | paste -d " " - - >> outfile_joined_tsp.log
