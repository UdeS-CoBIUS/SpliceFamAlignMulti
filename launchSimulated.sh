#!/bin/bash

for family in small medium large
do
    for i in `seq 1 20`
    do
        echo  examples/input/simulated/$family/_iteration_${i}
	echo ""
	echo  "..." generate pairwise alignments
	python2.7 src/main.py -ce Yes -s 2  -sf examples/input/simulated/$family/_iteration_${i}_initialsource.fasta -tf examples/input/simulated/$family/_iteration_${i}_target.fasta -s2tf examples/input/simulated/$family/_iteration_${i}_initialsource2target.txt -sef examples/input/simulated/$family/_iteration_${i}_initialsourceexonlist.txt -op examples/output/simulated/$family/_iteration_${i}_ -of list
	echo ""
	echo  "..." generate multiple alignment with SFAM_tcoffee_p
	python3 src_tcoffee/main.py -sf examples/input/simulated/$family/_iteration_${i}_initialsource.fasta -tf examples/input/simulated/$family/_iteration_${i}_target.fasta -s2tf examples/input/simulated/$family/_iteration_${i}_initialsource2target.txt -sef examples/input/simulated/$family/_iteration_${i}_initialsourceexonlist.txt -palnf examples/output/simulated/$family/_iteration_${i}_result.txt -psegf examples/output/simulated/$family/_iteration_${i}_segment.txt -op examples/output/simulated/$family/_iteration_${i}_tcoffee_p_

	echo ""
	echo  "..." generate multiple alignment with SFAM_tcoffee_m
	python3 src_tcoffee/main.py -sf examples/input/simulated/$family/_iteration_${i}_initialsource.fasta -tf examples/input/simulated/$family/_iteration_${i}_target.fasta -s2tf examples/input/simulated/$family/_iteration_${i}_initialsource2target.txt -sef examples/input/simulated/$family/_iteration_${i}_initialsourceexonlist.txt -palnf examples/output/simulated/$family/_iteration_${i}_result.txt -op examples/output/simulated/$family/_iteration_${i}_tcoffee_m_

	echo ""
	echo  "..." generate multiple alignment with SFAM_mblock
	python3 src_mblock/main.py -sf examples/input/simulated/$family/_iteration_${i}_initialsource.fasta -tf examples/input/simulated/$family/_iteration_${i}_target.fasta -s2tf examples/input/simulated/$family/_iteration_${i}_initialsource2target.txt -sef examples/input/simulated/$family/_iteration_${i}_initialsourceexonlist.txt -palnf examples/output/simulated/$family/_iteration_${i}_result.txt -op examples/output/simulated/$family/_iteration_${i}_
	echo ""
	echo ""
    done
done

