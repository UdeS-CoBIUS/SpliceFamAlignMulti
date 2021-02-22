# SpliceFamAlignMulti
Program Version 2020

computes multiple spliced alignment of all CDS and genes of a gene family, computes CDS orthology groups, and predicts new CDS by homology
----------------------------------------------------------------

Authors: Safa Jammali, Abigail Djossou et al.
Université de Sherbrooke, Canada
Cobius Lab:  https://cobius.usherbrooke.ca/
for questions email us at Safa.Jammali@USherbrooke.ca


### Requirements:

-Python2.7
-Python3
-argparse
-glob
-numpy
-multiprocessing
-functools
-Bio
-Ete3
-scipy
-skbio
-Muscle standalone
-Blast+ standalone
-Mafft standalone
-TCoffee standalone
-Splign standalone

### Usage to compute the spliced alignments for the examples contained in directory examples/input

####  Command line for simulated gene families contained in examples/input/simulated/
```
./launchSimulated.sh
```

Results are written in directory examples/output/simulated.

#### Command line for real gene families contained in examples/input/real/97/
```
./launchReal.sh
```

Results are written in directory examples/output/real/97.

### Usage to compute pairwise spliced alignments for any gene family
```
usage: src/main.py [-h] [-c CHOICESTRUCTURE] [-pm PAIRWISEMETHOD] [-s STEP]
               [-ce COMPAREEXON] [-sf SOURCEFILE] [-tf TARGETFILE]
               [-s2tf SOURCE2TARGETFILE] [-sef SOURCEEXONFILE]
               [-op OUTPUTPREFIX] [-of OUTPUTFORMAT] [-oa OUTPUTALIGNMENT]

```
  *-h*, --help            show help message and exit

  *-c CHOICESTRUCTURE*, --choiceStructure CHOICESTRUCTURE
                        Method used to infer splicing structure when a
                        splicing structure file is not given: blast or splign

  *-pm PAIRWISEMETHOD*, --pairwiseMethod PAIRWISEMETHOD
                        Method used to compute pairwise spliced alignment: sfa
                        or splign

  *-s STEP*, --step STEP  The method goes until Step 1, 2, or 3: 1, 2 or 3
                        (required)

  *-ce COMPAREEXON*, --compareExon COMPAREEXON
                        The method includes in Step2 a comparison of exons:
                        Yes or No (required)

  *-sf SOURCEFILE*, --sourceFile SOURCEFILE
                        Source (CDS) file name (required)

  *-tf TARGETFILE*, --targetFile TARGETFILE
                        Target (gene) file name (required)

  *-s2tf SOURCE2TARGETFILE*, --source2TargetFile SOURCE2TARGETFILE
                        Association between source and target file name
                        (required)

  *-sef SOURCEEXONFILE*, --sourceExonFile SOURCEEXONFILE
                        Source Exon (splicing structure) file name

  *-op OUTPUTPREFIX*, --outputPrefix OUTPUTPREFIX
                        Output prefix (required)

  *-of OUTPUTFORMAT*, --outputFormat OUTPUTFORMAT
                        Output format : list or aln (required)

  *-oa OUTPUTALIGNMENT*, --outputAlignment OUTPUTALIGNMENT
                        Output alignment : method used to compute pairwise
                        block alignment: nm or fsepsa

### Computing pairwise alignments: 

#### Command line

```
python2.7 src/main.py -ce Yes -s 2  -sf path_to_cds_fasta_file -tf path_to_gene_fasta_file -s2tf path_to_mapping_cds_gene_file  -sef path_to_cds_exon_structure_file  -op output_path_prefix -of list 
```

#### Example of command line
```
python3 src/main.py -ce Yes -s 2  -sf examples/input/simulated/small/_iteration_1_initialsource.fasta -tf examples/input/simulated/small/_iteration_1_target.fasta -s2tf examples/input/simulated/small/_iteration_1_initialsource2target.txt -sef examples/input/simulated/small/_iteration_1_initialsourceexonlist.txt -op examples/output/simulated/small/_iteration_1_ -of list
```
Results are written in directory examples/output/simulated/small with prefix \_iteration_1\_.

### Usage to compute multiple spliced alignments using SFAM_tcoffee for any gene family
```
usage: src_tcoffee/main.py [-h] [-idty IDENTITYTHRESHOLD] [-sf SOURCEFILE]
               [-tf TARGETFILE] [-s2tf SOURCE2TARGETFILE]
               [-sef SOURCEEXONFILE] [-palnf PAIRWISEALNFILE]
               [-psegf PAIRWISESEGMENTFILE] [-ba BLOCKALIGNMENT]
               [-op OUTPUTPREFIX] [-oa OUTPUTALIGNMENT] [-ce COMPAREEXON]
               [-msa MSAMETHOD]

```
  *-h*, --help            show help message and exit

  *-idty IDENTITYTHRESHOLD*, --identityThreshold IDENTITYTHRESHOLD
                        Identity threshold: real between 0.0 and 1.0 (default
                        = 0.3)

  *-sf SOURCEFILE*, --sourceFile SOURCEFILE
                        Source file name (required)

  *-tf TARGETFILE*, --targetFile TARGETFILE
                        Target file name (required)

  *-s2tf SOURCE2TARGETFILE*, --source2TargetFile SOURCE2TARGETFILE
                        Source to target file name (required)

  *-sef SOURCEEXONFILE*, --sourceExonFile SOURCEEXONFILE
                        Source exon file name

  *-palnf PAIRWISEALNFILE*, --pairwiseAlnFile PAIRWISEALNFILE
                        Pairwise alignment file name (required)

  *-psegf PAIRWISESEGMENTFILE*, --pairwiseSegmentFile PAIRWISESEGMENTFILE
                        Pairwise segment file name (required)

  *-ba BLOCKALIGNMENT*, --blockAlignment BLOCKALIGNMENT
                        Method for block alignment

  *-op OUTPUTPREFIX*, --outputPrefix OUTPUTPREFIX
                        Output prefix (required)

  *-oa OUTPUTALIGNMENT*, --outputAlignment OUTPUTALIGNMENT
                        Method for final alignment

  *-ce COMPAREEXON*, --compareExon COMPAREEXON
                        The method includes a final step that compares exons
                        of blocks for further merges: Yes or No

  *-msa MSAMETHOD*, --msaMethod MSAMETHOD
                        Multiple sequence aligner: muscle or mafft

### Computing multiple alignment using SFAM_tcoffee_p: 

#### Command line

```
python3 src_tcoffee/main.py -sf path_to_cds_fasta_file -tf path_to_gene_fasta_file -s2tf path_to_mapping_cds_gene_file  -sef path_to_cds_exon_structure_file -palnf output_path_prefix_result.txt -psegf output_path_prefix_segment.txt -op output_path_prefix -of list 
```

#### Example of command line
```
python3 src_tcoffee/main.py -sf examples/input/simulated/small/_iteration_1_initialsource.fasta -tf examples/input/simulated/small/_iteration_1_target.fasta -s2tf examples/input/simulated/small/_iteration_1_initialsource2target.txt -sef examples/input/simulated/small/_iteration_1_initialsourceexonlist.txt -palnf examples/output/simulated/small/_iteration_1_result.txt -psegf examples/output/simulated/small/_iteration_1_segment.txt -op examples/output/simulated/small/_iteration_1_tcoffee_p_ -of list
```

Results are written in directory examples/output/simulated/small with prefix \_iteration_1_tcoffee_p\_.

### Computing multiple alignment using SFAM_tcoffee_m: 

#### Command line

```
python3 src_tcoffee/main.py -sf path_to_cds_fasta_file -tf path_to_gene_fasta_file -s2tf path_to_mapping_cds_gene_file  -sef path_to_cds_exon_structure_file -palnf output_path_prefix_result.txt -op output_path_prefix -of list 
```

#### Example of command line
```
python3 src_tcoffee/main.py -sf examples/input/simulated/small/_iteration_1_initialsource.fasta -tf examples/input/simulated/small/_iteration_1_target.fasta -s2tf examples/input/simulated/small/_iteration_1_initialsource2target.txt -sef examples/input/simulated/small/_iteration_1_initialsourceexonlist.txt -palnf examples/output/simulated/small/_iteration_1_result.txt -op examples/output/simulated/small/_iteration_1_tcoffee_m_ -of list
```

Results are written in directory examples/output/simulated/small with prefix \_iteration_1_tcoffee_m\_.
 
### Usage to compute multiple spliced alignments using SFAM_mblock  for any gene family
```
usage: src_mblock/main.py [-h] [-idty IDENTITYTHRESHOLD]
               [-sf SOURCEFILE] [-tf TARGETFILE] [-s2tf SOURCE2TARGETFILE]
               [-sef SOURCEEXONFILE] [-palnf PAIRWISEALNFILE]
               [-op OUTPUTPREFIX] [-ce COMPAREEXON] [-msa MSAMETHOD]

```
  *-h*, --help            show help message and exit

  *-idty IDENTITYTHRESHOLD*, --identityThreshold IDENTITYTHRESHOLD
                        Identity threshold: real between 0.0 and 1.0 (default
                        = 0.3)

  *-sf SOURCEFILE*, --sourceFile SOURCEFILE
                        Source file name (required)

  *-tf TARGETFILE*, --targetFile TARGETFILE
                        Target file name (required)

  *-s2tf SOURCE2TARGETFILE*, --source2TargetFile SOURCE2TARGETFILE
                        Source to target file name (required)

  *-sef SOURCEEXONFILE*, --sourceExonFile SOURCEEXONFILE
                        Source exon file name

  *-palnf PAIRWISEALNFILE*, --pairwiseAlnFile PAIRWISEALNFILE
                        Pairwise alignment file name (required)

  *-op OUTPUTPREFIX*, --outputPrefix OUTPUTPREFIX
                        Output prefix (required)

  *-ce COMPAREEXON*, --compareExon COMPAREEXON
                        The method includes a final step that compares exons
                        of blocks for further merges: Yes or No

  *-msa MSAMETHOD*, --msaMethod MSAMETHOD
                        Multiple sequence aligner: muscle or mafft

### Computing multiple alignment using SFAM_mblock: 

#### Command line

```
python3 src_mblock/main.py -sf path_to_cds_fasta_file -tf path_to_gene_fasta_file -s2tf path_to_mapping_cds_gene_file  -sef path_to_cds_exon_structure_file -palnf output_path_prefix_result.txt -op output_path_prefix -of list 
```

#### Example of command line
```
python3 src_mblock/main.py -sf examples/input/simulated/small/_iteration_1_initialsource.fasta -tf examples/input/simulated/small/_iteration_1_target.fasta -s2tf examples/input/simulated/small/_iteration_1_initialsource2target.txt -sef examples/input/simulated/small/_iteration_1_initialsourceexonlist.txt -palnf examples/output/simulated/small/_iteration_1_result.txt -op examples/output/simulated/small/_iteration_1_mblock_ -of list
```

Results are written in directory examples/output/simulated/small with prefix \_iteration_1_mblock\_.

