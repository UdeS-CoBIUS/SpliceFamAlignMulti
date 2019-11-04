#!/bin/bash

# Pairwise FAM86SG
python src/main.py -s 2 -ce Yes -sf examples/input/FAM86SG/FAM86SG_initialsource.fasta -tf examples/input/FAM86SG/FAM86SG_target.fasta -s2tf examples/input/FAM86SG/FAM86SG_initialsource2target.txt -sef examples/input/FAM86SG/FAM86SG_initialsourceexonlist.txt -op examples/output/FAM86SG_ -of list

# Multi FAM86SG
python3 src_multi/main.py -sf examples/input/FAM86SG/FAM86SG_initialsource.fasta -tf examples/input/FAM86SG/FAM86SG_target.fasta -s2tf examples/input/FAM86SG/FAM86SG_initialsource2target.txt -sef examples/input/FAM86SG/FAM86SG_initialsourceexonlist.txt -palnf examples/output/FAM86SG_result.txt -op examples/output/FAM86SG_ -ce Yes

# Pairwise FAM86
python src/main.py -s 2 -ce Yes -sf examples/input/FAM86/FAM86_initialsource.fasta -tf examples/input/FAM86/FAM86_target.fasta -s2tf examples/input/FAM86/FAM86_initialsource2target.txt -sef examples/input/FAM86/FAM86_initialsourceexonlist.txt -op examples/output/FAM86_ -of list

# Multi FAM86
python3 src_multi/main.py -sf examples/input/FAM86/FAM86_initialsource.fasta -tf examples/input/FAM86/FAM86_target.fasta -s2tf examples/input/FAM86/FAM86_initialsource2target.txt -sef examples/input/FAM86/FAM86_initialsourceexonlist.txt -palnf examples/output/FAM86_result.txt -op examples/output/FAM86_ -ce Yes

python src/main.py -s 2 -ce Yes -sf examples/input/MAG/MAG_initialsource.fasta -tf examples/input/MAG/MAG_target.fasta -s2tf examples/input/MAG/MAG_initialsource2target.txt -sef examples/input/MAG/MAG_initialsourceexonlist.txt -op examples/output/MAG_ -of list

python3 src_multi/main.py -sf examples/input/MAG/MAG_initialsource.fasta -tf examples/input/MAG/MAG_target.fasta -s2tf examples/input/MAG/MAG_initialsource2target.txt -sef examples/input/MAG/MAG_initialsourceexonlist.txt -palnf examples/output/MAG_result.txt -op examples/output/MAG_ 

python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource.fasta -tf examples/input/ENSGT00940000158754/ENSGT00940000158754_target.fasta -s2tf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource2target.txt -sef examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsourceexonlist.txt -op examples/output/ENSGT00940000158754_ -of list

python3 src_multi/main.py -sf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource.fasta -tf examples/input/ENSGT00940000158754/ENSGT00940000158754_target.fasta -s2tf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource2target.txt -sef examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000158754_result.txt -op examples/output/ENSGT00940000158754_ -ce Yes

python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource.fasta -tf examples/input/ENSGT00940000162816/ENSGT00940000162816_target.fasta -s2tf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource2target.txt -sef examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsourceexonlist.txt -op examples/output/ENSGT00940000162816_ -of list

python3 src_multi/main.py -sf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource.fasta -tf examples/input/ENSGT00940000162816/ENSGT00940000162816_target.fasta -s2tf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource2target.txt -sef examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000162816_result.txt -op examples/output/ENSGT00940000162816_  -ce Yes

