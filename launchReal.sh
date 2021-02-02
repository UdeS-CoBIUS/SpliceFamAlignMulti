#!/bin/bash

for f in ENSGT00530000063205 ENSGT00390000008371 ENSGT00390000000715 ENSGT00940000157909 ENSGT00530000063187 ENSGT00950000182978 ENSGT00950000182681 ENSGT00950000182875 ENSGT00950000182728 ENSGT00950000182931 ENSGT00950000182705 ENSGT00530000063023 ENSGT00940000153241 ENSGT00950000182956 ENSGT00950000182783 ENSGT00950000183192 ENSGT00950000182727
do
    echo $f
    echo ""
    echo  "..." generate pairwise alignments
    python2.7 src/main.py -ce Yes -s 2  -sf examples/input/real/97/${f}_initialsource.fasta -tf examples/input/real/97/${f}_target.fasta -s2tf examples/input/real/97/${f}_initialsource2target.txt -sef examples/input/real/97/${f}_initialsourceexonlist.txt -op examples/output/real/97/${f}_ -of list

    echo ""
    echo  "..." generate multiple alignment with SFAM_mblock
    python3 src_mblock/main.py -sf examples/input/real/97/${f}_initialsource.fasta -tf examples/input/real/97/${f}_target.fasta -s2tf examples/input/real/97/${f}_initialsource2target.txt -sef examples/input/real/97/${f}_initialsourceexonlist.txt -palnf examples/output/real/97/${f}_result.txt -op examples/output/real/97/${f}_
    echo ""
    echo ""
done

