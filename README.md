# SpliceFamAlignMulti
Program Version 2019

Performs a multiple spliced alignment from the set of all pairs of source CDS and target gene comparison in a progressive way.
----------------------------------------------------------------

Authors: Aida Ouangraoua and Abigail Djossou

Université de Sherbrooke, Canada
CoBIUS Lab:  https://cobius.usherbrooke.ca/

For questions email us at abigail.djossou@usherbrooke.ca

### Requirements:

-PyCogent
-Bio
-Ete toolkit
-Blast+ standalone
-Splign tool
-argparse


### Usage
```
usage: main.py [-h] [-c CHOICE] [-f FORCE] [-m METHOD] [-it INPUTTYPE]
               [-gn GENENUMBER] [-sf SOURCEFILE] [-tf TARGETFILE]
               [-s2tf SOURCE2TARGETFILE] [-sef SOURCEEXONFILE]
               [-s1 FIRSTSPECIES] [-s2 SECONDSPECIES] [-gid1 FIRSTGENEID]
               [-gid2 SECONDGENEID] [-gidlf GENEIDLISTFILE] [-g GENE]
               [-slf SPECIESLISTFILE] [-op OUTPUTPREFIX] [-of OUTPUTFORMAT]

```
 *-h*, --help   show this help message and exit
  
  *-c , --choice \<CHOICE>*         used where there 's no structre information, it can be blast or splign  

  *-f , --force \<FORCE>*           use SFA_G, it can be Yes or No   
  
  *-m , --method \<METHOD>*         SFA  

  *-it , --inputType \<INPUTYPE>*   file id or name   
  
  *-gn , --geneNumber \<GENENUMBER>*         Gene number: pairwise or multiple  
  
  *-sf , --sourceFile \<SOURCEFILE>*         Source file name  
  
  *-tf , --targetFile \<TARGETFILE>*         Target file name  
  
  *-s2tf , --source2TargetFile \<SOURCE2TARGETFILE>*         Source to target file name   
  
  *-sef , --sourceExonFile \<SOURCEEXONFILE>*         Source Exon file name  
  
  *-s1 , --firstSpecies \<FIRSTSPECIES>*         First species common name (required if --inputType = file or name,
                                                  and --geneNumber = pairwise")    
  
  *-s2 , -secondSpecies \<SECONDPECIES>*         Second species common name (required if --inputType = file or name,
                                                  and --geneNumber = pairwise")    
  
  *-gid1 , --firstGeneId \<FIRSTGENEID>*         First gene Ensembl Id (required if --inputType = id, and 
						--geneNumber = pairwise)    

  *-gid2 , --secondGeneId \<SECONGENEID>*         Second gene Ensembl Id (required if --inputType = id, and 
						--geneNumber = pairwise)  

  *-gidlf , --geneIdListFile \<GENEIDLISTFILE>*         Gene Id list file name (required if --inputType = id,
							 and --geneNumber = multiple)  

  *-g , -- gene \<GENE>*        Gene common name (required if --inputType = name) 
   

  *-slf , --speciesListFile \<SPECIESLISTFILE>*         Species list file name (required if --inputType = name,
                                                and --genenumber = multiple)    
  
  *-o , --outfile \<OUTFILE> *      output file   

  *-of , --outformat \<OUTFORMAT> *      output file format (list or aln)   

### Running SpliceFamAlignMulti: examples of command line

#### To compute pairwise:
```
python src/main.py -s 2 -ce Yes -sf _path-to-fasta-file_ -tf _path-to-gene-fasta-file_ -s2tf _path-to-association-cds-vs-gene-file_ -sef _path-to-exon-list-file_  -op _output-path_ -of list
```
#### To compute multiple spliced alignment:
```
python3 src_multi/main.py -ce Yes -sf _path-to-fasta-file_ -tf _path-to-gene-fasta-file_ -s2tf _path-to-association-cds-vs-gene-file_ -sef _path-to-exon-list-file_ -palnf _path-to-result-file_ -op _output-path_
```
