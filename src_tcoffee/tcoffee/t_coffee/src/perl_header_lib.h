char *PerlScriptName[]={"rec_sum.pl","count.pl","p\
rocess_list.pl","make_license.pl","CCsed.script","\
msa2bootstrap.pl","tc_generic_method.pl","dynamic.\
pl","rnapdb2protpdb.pl","generic_method.tc_method"\
,"clustalw_method.tc_method","extract_from_pdb","i\
nstall.pl","clean_cache.pl","nature_protocol.pl","\
mocca","dalilite.pl","wublast.pl","blastpgp.pl","n\
cbiblast_lwp.pl","wublast_lwp.pl","RNAplfold2tclib\
.pl","fasta_seq2RNAplfold_templatefile.pl","fasta_\
seq2hmmtop_fasta.pl","fasta_seq2consan_aln.pl","cl\
ustalw_aln2fasta_aln.pl","seq2name_seq.pl","seq2in\
tersection.pl","msf_aln2fasta_aln.pl","msa.pl","up\
p.pl","clustalo.pl","dca.pl","blast_aln2fasta_aln.\
pl","blast_xml2fasta_aln.pl","fasta_aln2fasta_aln_\
unique_name.pl","newick2name_list.pl","excel2fasta\
.pl","nameseq2fasta.pl","any_file2unix_file.pl","E\
ndList"};char *PerlScriptFile[]={"use File::Copy;\\
nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(\
USER);\n$x_field=0;\n$y_field=1;\n$y_field_set=1;\\
n$nyf=1;\n\n$interval=0;\n$file=\"stdin\";\n\n$pri\
nt_avg=1;\n$print_sd=0;\n$print_sum=0;\n$print_n=0\
;\nforeach $value ( @ARGV)\n    {\n	if ($value ne \
$ARGV[$np]) \n	    {\n	    ;\n	    }\n	elsif($valu\
e eq \"-s\")\n	     {\n	       $step=$ARGV[++$np];\
\n	       $np++;\n	     }\n	elsif($value eq \"-pri\
nt_all\")\n	    {\n	    $print_sd=$print_avg=$prin\
t_n=$print_sum=1;\n	    $np++;\n	    }\n	elsif($va\
lue eq \"-print_sum\")\n	    {\n	    $print_sum=1;\
\n	    $print_avg=0;\n	    $np++;\n	    }\n	elsif(\
$value eq \"-print_n\")\n	    {\n	    $print_n=1;\\
n	    $print_avg=0;\n	    $np++;\n	    }\n	elsif($\
value eq \"-print_avg\")\n	    {\n	    $print_avg=\
1;\n	    $print_avg=0;\n	    $np++;\n	    }\n	elsi\
f($value eq \"-sd\")\n	    {\n	    $print_sd=1;\n	\
    $print_avg=0;\n	    $np++;\n	    }\n	elsif($va\
lue eq \"-h\")\n	    {\n	    $header=1;\n	    $np+\
+;\n	    }\n	elsif ($value eq \"-i\")\n	    {\n	  \
  $interval= $ARGV[++$np];\n	    $np++;\n    	    \
}\n	elsif ($value eq \"-r2\")\n	    {\n	      $r2=\
1;\n	      \n	      $np  = $ARGV[++$np];\n	      $\
nsim= $ARGV[++$np];\n	      $np++;\n    	    }\n	e\
lsif ($value eq \"-r\")\n	    {\n	    $min= $ARGV[\
++$np];\n	    $max= $ARGV[++$np];\n	    $np++;\n  \
  	    }\n	\n	elsif ($value eq \"-x\")\n	    {\n	 \
   $x_field= $ARGV[++$np]-1;\n	    $np++;\n    	  \
  }\n	elsif ($value eq \"-y\")\n	    {\n	    $nyf=\
0;  \n	    while ($ARGV[$np+1] && !($ARGV[$np+1]=~\
/\\-/))\n	      {\n		$y_field[$nyf++]=$ARGV[++$np]\
-1;\n		$y_field_set=1;\n	      }\n\n	    $np++;\n \
   	    }\n	elsif ($value eq \"-file\")\n	    {\n	\
    $file= $ARGV[++$np];\n	    $file_set=1;\n	    \
$np++;\n    	    }       \n	elsif ( $value eq \"h\\
" ||  $value eq \"-h\" || $value eq \"-H\" || $val\
ue eq \"-help\" || $value eq \"help\")\n	  {\n	   \
 print STDOUT \"data_analyse: Analyse and discreti\
zation of data\\n\";\n	    print STDOUT \"       -\
file:    <file containing the data to analyze>,.<d\
ef=STDIN>\\n\";\n	    print STDOUT \"       -x: <f\
ield containing the X>,...............<Def=0>\\n\"\
;\n	    print STDOUT \"       -y: <field containin\
g the Y>,...............<Def=1>\\n\";\n	    print \
STDOUT \"       -i:<Interval size on the X>,......\
.........<Def=0>\\n\";\n	    print STDOUT \"      \
 -i:<0:only one interval>\\n\";\n	    print STDOUT\
 \"       -r:<Range of the X>\\n\";\n	    print ST\
DOUT \"       -s:<Step on the  X, 0 means non slid\
ing bins>\\n\";\n	    print STDOUT \"       -sd: p\
rint standard deviation on the Y\";\n	    print ST\
DOUT \"       -h  : print column header \\n\";\n	 \
   exit (0);\n	  }\n	elsif ($value=~/-/)\n	  {\n	 \
   print \"$value is not a valid FLAG[FATAL]\\n\";\
\n	    exit (0);\n	   } \n	elsif ($list eq \"\") \\
n	    {\n	    $file=$ARGV[$np];\n	    $np++;\n	   \
 }\n	\n	\n      }\n\n\n\n\n\nif ($file eq \"stdin\\
")\n	{\n	$remove_file=1;\n	$file=\"tmp$$\";\n	open\
 (F, \">$file\");\n	while (<STDIN>)\n		{\n		print \
F $_;\n		}\n	close (F);\n	 \n	;}\n\n\n\nif ($inter\
val && $step)\n  {\n    my $nl;\n    open(F,$file)\
;\n    while (<F>)\n      {\n	$line=$_;\n	\n	if (!\
/\\S/){next;}\n	@list=($line=~/(\\S+)/g);\n	$val{$\
nl}{x}=$list[$x_field];\n	$val{$nl}{y}=$list[$y_fi\
eld[0]];\n	$nl++\n      }\n    close (F);\n    \n \
   for (my $a=$min; $a<($max+$interval); $a+=$step\
)\n      {\n	my ($avgx, $avgy, $cn);\n	\n	my $rmin\
=$a-$interval;\n	my $rmax=$a;\n	$cn=0;\n	for (my $\
b=0; $b<$nl; $b++)\n	  {\n	    my $x=$val{$b}{x};\\
n	    my $y=$val{$b}{y};\n	    if ($x<=$rmax && $x\
>=$rmin)\n	      {\n		$avgx+=$x;\n		$avgy+=$y;\n		\
$cn++;\n		$tcn++;\n		$val{$b}{used}=1;\n	      }\n\
	  }\n	if ($cn)\n	  {\n	    $avgx/=$cn;\n	    $avg\
y/=$cn;\n	  }\n	printf \"%.3f %.3f %.3f\\n\", $avg\
x, $avgy, $avgx-$avgy;\n      }\n    for (my $a=0;\
 $a<$nl; $a++)\n      {\n	if ( !$val{$a}{used})\n	\
  {\n	    print \"---$val{$a}{x}; $val{$a}{y}\\n\"\
;\n	  }\n      }\n  }\nelse\n  {\n    if ($interva\
l && $max)\n      {\n	$interval_size=($max-$min)/$\
interval;\n      }\n    elsif ($interval)\n      {\
\n	open(F,$file);  \n	my $set_max=0;\n	my $set_min\
=0;\n	while (<F>)\n	  {\n	    my $v=$_;\n	    chom\
p($v);\n	    print \"--$v--\";\n	    \n	    if ($v\
<$min ||!$set_min){$set_min=1;$min=$v;}\n	    if (\
$v>$max ||!$set_max){$set_max=1;$max=$v;}\n	  }\n	\
close (F);\n	print \"$min $max uuuu\";\n	$interval\
_size=($max-$min)/$interval;\n      }\n    open(F,\
$file);  \n    while (<F>)\n      {\n	$line=$_;\n	\
if (!/\\S/){next;}\n	@list=($line=~/(\\S+)/g);\n	\\
n	if ($interval==0){$bin=0;}\n	else{$bin=int (($li\
st[$x_field]-$min)/($interval_size));}\n	\n	\n	if \
($bin && $bin==$interval){$bin--;}\n	for ( $a=0; $\
a<$nyf; $a++)\n	  {\n	    $sum{$a}{$bin}+=$list[$y\
_field[$a]];\n	    $sum2{$a}{$bin}+=$list[$y_field\
[$a]]*$list[$y_field[$a]];\n	    $n{$a}{$bin}++;\n\
	  }\n      }\n    \n    if (!$interval){$interval\
=1;}\n    for ( $a=0; $a<$interval; $a++)\n      {\
\n	printf ( \"%4d %4d \", $interval_size*$a, $inte\
rval_size*($a+1));\n	for ( $b=0; $b<$nyf; $b++)	\n\
	  {\n	    $i=$interval*$a;\n	    if ( $n{$b}{$a}=\
=0)\n	      {\n		$avg=0;\n		$sd=0;\n	      }\n	   \
 else\n	      {\n		$avg=$sum{$b}{$a}/$n{$b}{$a};\n\
		$sd=sqrt($sum2{$b}{$a}*$n{$b}{$a}-$sum{$b}{$a}*$\
sum{$b}{$a})/($n{$b}{$a}*$n{$b}{$a});\n	      }\n	\
    if ($print_n) {printf \"%15.4f \", $n{$b}{$a};\
}\n	    if ($print_sum){printf \"%15.4f \", $sum{$\
b}{$a};}\n	    if ($print_avg){printf \"%15.4f \",\
 $avg}\n	    if ($print_sd) {printf \"%15.4f \", $\
sd;}\n	  }\n	printf (\"\\n\");\n      }\n  }\n\nif\
 ( $remove_file){unlink $file;}\n","use File::Copy\
;\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env q\
w(USER);\n\nforeach $v (@ARGV){$cl.=$v;}\n\n\nif (\
 $cl=~/-k(\\d+)/){$k=$1;}\nelse {$k=1;}\nif ( $cl=\
~/-w(\\d+)/){$w=$1;}\nelse {$w=-1;}\nif ( $cl=~/-p\
(\\d+)/){$p=$1;}\nelse {$p=-1;}\n\nwhile (<STDIN>)\
\n  {\n    @l=($_=~/(\\S+)/g);\n    $v=$l[$k-1];\n\
    if ( !$h{$v}){@ll=($v, @ll);}\n    \n    if ( \
$w==-1)\n      {$h{$v}++;}\n    else\n      {$h{$v\
}+=$l[$w-1];}\n\n    if ($p!=-1){$print{$v}=$l[$p-\
1];}\n\n  }\nforeach $v (@ll)\n  {\n    print \"$v\
 $print{$v} $h{$v}\\n\";\n  }\n","\nuse Env qw(HOS\
T);\nuse Env qw(HOME);\nuse Env qw(USER);\n$random\
_tag=int (rand 10000)+1;\n$unique_prefix=\"$$.$HOS\
T.$random_tag\";\n$queue=\"distillery.and.mid\";\n\
$monitor=0;\n$stderr_file=\"/dev/null\";\n$stdio_f\
ile=\"/dev/null\";\n$log_file=\"/dev/null\";\n$pau\
se_time=0;\n$max_sub_jobs=60;\n$min_sub_jobs=30;\n\
$output_all=0;\n$var='\\$';\n\nforeach $value ( @A\
RGV)\n    {\n	if ($value ne $ARGV[$np]) \n	    {\n\
	    ;\n	    }\n	elsif ($value eq \"-max_sub_jobs\\
")\n	    {\n	    $max_sub_jobs= $ARGV[++$np];\n	  \
  $np++;\n    	    }	\n	elsif ($value eq \"-min_su\
b_jobs\" )\n	    {\n	    $min_sub_jobs= $ARGV[++$n\
p];\n	    $np++;\n    	    }\n	elsif ($value eq \"\
-para\")\n	    {\n	    $para=1;\n	    $monitor=1;\\
n	    $np++;\n    	    }\n	elsif ($value eq \"-mon\
itor\") \n	    {\n	    $monitor=1;\n	    $np++;\n	\
    }\n	elsif ($value eq \"-no_monitor\") \n	    {\
\n	    $monitor=0;\n	    $np++;\n	    }\n	elsif ($\
value eq \"-queue\")\n	    {\n	    $queue=$ARGV[++\
$np];\n	    $np++;\n	    }	\n	elsif ($value eq \"-\
stderr_file\")\n	    {\n	    $stderr_file=$ARGV[++\
$np];\n	    $np++;\n	    }\n	elsif ($value eq \"-s\
tdio_file\")\n	    {\n	    $stdio_file=$ARGV[++$np\
];\n	    $np++;\n	    }\n	elsif ($value eq \"-outp\
ut_all\")\n	    {\n	    $output_all=1;\n	    $np++\
;\n	    }\n	elsif ($value eq \"-pause\") \n	    {\\
n	    $pause_time=$ARGV[++$np];\n	    $np++;\n	   \
 }\n	elsif ($value eq \"-log\")\n	      {\n	      \
 $log=1;\n	       \n	       if ($ARGV[$np+1]=~/\\-\
\\S+/) \n	          {\n		  $log_file=\"stderr\";\n\
	          }\n	       else \n	          {\n		  $lo\
g_file=$ARGV[++$np]; \n		  ++$np;\n		 \n	         \
 }\n	      }\n	elsif ( $value eq \"-com\")\n	    {\
\n		\n		if (!$ARGV[$np+1]=~/^\\'/) { $com=$ARGV[++\
$np];}\n		else {$com=$ARGV[++$np];}\n\n	     $np++\
;\n	    }\n	elsif ( $value eq \"-check\")\n	  {\n	\
    \n	    if (!$ARGV[$np+1]=~/^\\'/) { $check=$AR\
GV[++$np];}\n	    else {$check=$ARGV[++$np];}\n	  \
  $np++;\n	  }\n	elsif ($com eq \"\") \n	    {\n	 \
   $com_set=1;\n	    $com=$ARGV[$np];\n	    \n	   \
 $np++;\n	    }\n	elsif ($list eq \"\") \n	    {\n\
	    $list_set=1;\n	    $list=$ARGV[$np];\n	    $n\
p++;\n	    }\n	elsif ( $var_set eq \"\")\n	    {\n\
	    $var_set=1;\n	    $var=$ARGV[$np];\n	    $np+\
+;\n	    }\n	}\n\n\n\n\nif ( $com eq \"\"){print \\
"You Need to Provide a Command [FATAL]\\n\";\n	   \
   die;\n	     }\n\n\n\nif ($list_set==0) \n    {\\
n    $x= int (rand 100000)+1;\n    $tmp_file_name=\
\"tmp_file_$x\";\n    open ( TMP, \">$tmp_file_nam\
e\");\n    while (<STDIN>)\n      {\n	print TMP $_\
;\n      }\n    close (TMP);\n    open (F, $tmp_fi\
le_name);\n    }\nelse \n    {\n    open (F, $list\
);\n    }\n\nif ($para==0) \n    {\n\n     @tc_lis\
t= <F>;\n     close (F); \n     \n     foreach $va\
l(@tc_list) \n	    {\n	      \n	      \n	      \n	\
      $loc_com=$com;\n	      if ($check){$loc_chec\
k=$check;}\n	      \n	      @i_val=($val=~/([^\\s]\
+)/g);\n	      \n	      if ( $#i_val==0)\n		{\n		 \
 if ($check){$loc_check=~s/$var/$i_val[0]/g;}\n		 \
 $loc_com=~s/$var/$i_val[0]/g;\n		}\n	      else\n\
		{\n		  for ($n=1; $n<=$#i_val+1;$n++ )\n		    {\\
n		      \n		      $sub=\"$var$n\";\n		      \n		 \
     $loc_com=~s/$sub/$i_val[$n-1]/g;\n		      if \
($check){$loc_check=~s/$var/$i_val[0]/g;}\n		    }\
\n		}\n	      if ( $check && -e $loc_check)\n		{\n\
		  print STDERR \"skipping $loc_com...\\n\";\n		 \
 }\n	      else\n		{\n		  system \"$loc_com\";\n		\
}\n	    }\n    exit;\n    }\n\nelsif ($para==1) \n\
    {\n    print STDERR \"do parallel execution of\
: \\\"$com $list\\\"\\n\";\n    \n    if ($log==1)\
 \n	{\n	if ($log_file eq \"stdout\" || $log_file e\
q \"stderr\" ) \n		{\n		$log_file=\"\";\n	        \
}\n\n        else \n		{\n		system \"echo LOG FILE>\
 $log_file\";\n		\n	        }\n	}\n    else	\n	{\n\
	open ( OUT, \">/dev/null\");\n	}\n	\n    \n    $i\
d=0;\n    $n_sub=0;\n    while ($val=<F>) \n	    {\
	    	    \n	    $job_log[$id]=\"$HOME/tmp/$unique\
_prefix.$id.log_file\";\n	    \n	    $job=$unique_\
prefix.\"_$id\";\n	    open (JOB, \">$job\");\n	  \
  \n	    $loc_com=$com;\n	    chop $val;\n\n	    $\
loc_com=~s/\\$/$val/g;\n	 \n	    print JOB \"#!/bi\
n/csh\\n\";\n	    print JOB \"#\\$ -cwd\\n\";\n	  \
  print JOB \"#\\$ -N $unique_prefix\\n\";\n	    i\
f ($queue && !($queue eq \" \")) {print JOB \"#\\$\
 -l $queue\\n\";}\n	    print JOB \"#\\n\";	    \n\
            print JOB \"$loc_com\\n\";\n	    print\
 JOB \"echo FINISHED  >> $job_log[$id]\\n\";\n	   \
 print JOB \"pwd\\n\";\n	    \n	    close (JOB);\n\
	    if ( $output_all==1)\n		{\n		system \"qsub $j\
ob >  $unique_prefix\";		\n	        }\n	    else\n\
		{system \"qsub $job -e $stderr_file -o $stdio_fi\
le >$unique_prefix\";	        \n	        } \n\n\n\\
n	    print STDERR \"$id: $output_all\\n\";\n	    \
$n_sub++;\n	    if ( $max_sub_jobs && $n_sub==$max\
_sub_jobs) \n		{\n		$n_sub=monitor_process($min_su\
b_jobs,@job_log); 		 \n		\n	        }	\n	   \n    \
        unlink $unique_prefix;\n	    sleep $pause_\
time;\n	    $id++;\n	    }\n\n    close (OUT);\n  \
  close (F);\n\n    print STDERR \"Your $id Jobs H\
ave Been Submited (NAME=$unique_prefix)\\n\";\n   \
 monitor_process (0, @job_log);\n    foreach $file\
(@job_log) {if (-e $file) {unlink($file);}}\n    \\
n    }\n\nsub monitor_process ( @job_list)\n    {\\
n    my (@job_list)=@_;\n    my $min_sub_jobs=shif\
t (@job_list);\n    my $n_sub_jobs;\n    my $finis\
hed;\n    my $n=0;\n\n    $n_sub_jobs=-1;\n    $fi\
nished=0;\n    print STDERR \"\\nMonitor Batch: [$\
min_sub_jobs]\";\n       \n    while (!$finished &\
& (($n_sub_jobs>$min_sub_jobs)|| $n_sub_jobs==-1) \
) \n	{\n	$finished=1;\n	$n_sub_jobs=0;\n	$n=0;\n	f\
oreach $file (@job_list)\n	        {\n	\n		if (-e \
$file){;}\n		else \n		    {\n		    $finished=0; $n\
_sub_jobs++;\n	            }\n	        }\n	system \
\"sleep 1\";\n        }\n    \n    return $n_sub_j\
obs;\n    }\n    \n    \nif ($tmp_file_name){unlin\
k($tmp_file_name);}\n","\n\nforeach ($np=0; $np<=$\
#ARGV; $np++)\n    {\n    $value=$ARGV[$np];\n\n  \
  if ($value eq \"-file\")\n      {\n      $file= \
$ARGV[++$np];\n      }\n    elsif ($value eq \"-ty\
pe\")\n      {\n        $type= $ARGV[++$np];\n    \
  }\n    elsif ($value eq \"-institute\")\n      {\
\n        $institute= $ARGV[++$np];\n      }\n    \
elsif ($value eq \"-author\")\n      {\n        $a\
uthor= $ARGV[++$np];\n      }\n    elsif ($value e\
q \"-date\")\n      {\n        $date= $ARGV[++$np]\
;\n      }\n     elsif ($value eq \"-program\")\n \
     {\n        $program= $ARGV[++$np];\n      }\n\
    elsif ($value eq \"-email\")\n      {\n       \
 $email= $ARGV[++$np];\n      }\n    else\n      {\
\n	print \"$value is an unkown argument[FATAL]\\n\\
";\n	exit (1);\n      }\n  }\n\n\n\nopen F, $file \
|| die;\nprint $INSTITUTE;\nif ( $type eq \"c\"){p\
rint \"/******************************COPYRIGHT NO\
TICE*******************************/\\n\";}\nif ( \
$type eq \"perl\"){print \"#######################\
#######COPYRIGHT NOTICE###########################\
###/\\n\";}\nif ( $type eq \"txt\"){print \"------\
-------------------------COPYRIGHT NOTICE---------\
---------------------/\\n\";}\n\n\nwhile (<F>)\n  \
{\n  s/\\$INSTITUTE/$institute/g;\n  s/\\$AUTHOR/$\
author/g;\n  s/\\$DATE/$date/g;\n  s/\\$PROGRAM/$p\
rogram/g;  \n  s/\\$EMAIL/$email/g;  \n  if ( $typ\
e eq \"txt\"){print $_;}\n  elsif ($type eq \"c\")\
{chop $_; print \"\\/*$_*\\/\\n\";}\n  elsif ($typ\
e eq \"perl\"){print \"\\#$_\";}\n}\nclose (F);\ni\
f ( $type eq \"c\"){print \"/*********************\
*********COPYRIGHT NOTICE*************************\
******/\\n\";}\nif ( $type eq \"perl\"){print \"##\
############################COPYRIGHT NOTICE######\
########################/\\n\";}\nif ( $type eq \"\
txt\"){print \"-------------------------------COPY\
RIGHT NOTICE------------------------------/\\n\";}\
\n\n","\nwhile (<>)	\n	{\n	s/\\=cc/123456789/g;\n	\
s/\\bcc/\\$\\(CC\\)/g;\n	s/123456789/\\=cc/g;\n	pr\
int $_;\n	}\n\n","$version=\"1.00\";\n$rseed= int(\
rand(100000))+1;\n\n\nif ( $#ARGV==-1)\n  {\n    p\
rint \"msa2bootstrap -i <input_file> -input <seq|m\
sa|matrix|tree> -n <N-Boostrap> -o <outtree> -tmod\
e <nj|upgma|parsimony|ml> -dmode <kimura> -alignpg\
 <t_coffee | muscle | clustalw> -rtree <file> -sty\
pe <prot|cdna|dna> -recompute -system <cygwin|unix\
>\";\n    print \"\\n\\t-i: input file, can be seq\
uneces, msa, matrix, trees, type is specified via \
-input\";\n    print \"\\n\\t-input: Type of input\
 data\";\n    print \"\\n\\t\\tmsa: msa in fasta f\
ormat\";\n    print \"\\n\\t\\tseq: compute an msa\
 with -alignpg\";\n    print \"\\n\\t\\tmatrix: ph\
ylipp distance matrix fed directly to method -tmod\
e [caveat: tmode=nj or upgma]\";\n    print \"\\n\\
\t\\ttree: list of newick trees directly fed to co\
nsence in order to generate a bootstraped tree\";\\
n    \n    print \"\\n\\t-n: number of bootstrap r\
eplicates\";\n    print \"\\n\\t-o: name of the ou\
tput tree. Files are not overwritten. Use -recompu\
te to overwrite existing file\";\n    print \"\\n\\
\t-tmode: tree mode: nj|upgma|parsimony|ml\";\n   \
 print \"\\n\\t-dmode: distance mode\";\n    print\
 \"\\n\\t-alignpg: program for aligning sequences \
(t_coffee=default)\";\n    print \"\\n\\t-rtree: r\
eplicate tree file (default: no file)\";\n    prin\
t \"\\n\\t-rmsa: replicate msa file (default: no f\
ile)\";\n    print \"\\n\\t-rmat: replicate matrix\
 file (default: no file)\";\n    print \"\\n\\t-st\
ype: sequence type: protein, dna or cdna\";\n    p\
rint \"\\n\\t-recompute: force files to be overwri\
tten\";\n    print \"\\n\\t-system: cygwin|unix\";\
\n      \n\n    \n    &my_exit (EXIT_FAILURE);\n  \
}\nforeach $arg (@ARGV){$command.=\"$arg \";}\n\np\
rint \"CLINE: $command\\n\";\n$threshold=100;\n$tr\
im_msa=0;\n$stype=\"prot\";\nprint \"msa2bootstrap\
 \";\n\n$system=\"cygwin\";\nif(($command=~/\\-sys\
tem (\\S+)/))\n  {\n    $system=$1;\n    if ( $sys\
tem eq \"cygwin\")\n      {\n	$exec_extension=\".e\
xe\";\n      }\n    elsif ( $system eq \"unix\")\n\
      {\n	$exec_extension=\"\";\n	print \"system=U\
nix\";die;\n      }\n    else\n      {\n	print \"m\
sa2boostrap: -system=$system is an unknown mode [F\
ATAL]\\n\"; die;\n      }\n    \n    print \"-syst\
em $system \";\n  }\nif(($command=~/\\-stype (\\S+\
)/))\n  {\n    $stype=$1;\n  }\nprint \"-stype=$st\
ype \";\n\n\n\nif(($command=~/\\-i (\\S+)/))\n  {\\
n    $msa=$1;\n    print \"-i $msa \";\n  }\n\nif(\
($command=~/\\-rtree (\\S+)/))\n  {\n    $rtree=$1\
;\n    print \"-rtree=$rtree \";\n  }\n\nif(($comm\
and=~/\\-rmsa (\\S+)/))\n  {\n    $rmsa=$1;\n  }\n\
if(($command=~/\\-rmat (\\S+)/))\n  {\n    $rmat=$\
1;\n  }\n$input=\"seq\";\nif(($command=~/\\-input \
(\\S+)/))\n  {\n    $input=$1;\n  }\nprint \"-inpu\
t=$input \";\n\n$dmode=\"kimura\";\nif(($command=~\
/\\-dmode (\\S+)/))\n  {\n    $dmode=$1;\n  }\npri\
nt \"-dmode=$dmode \";\n$alignpg=\"muscle\";\nif((\
$command=~/\\-alignpg (\\S+)/))\n  {\n    $alignpg\
=$1;\n  }\nprint \"-alignpg=$dmode \";\n\n$tmode=\\
"nj\";\nif(($command=~/\\-tmode (\\S+)/))\n  {\n  \
  $tmode=$1;\n  }\nprint \"-tmode=$tmode \";\n$rec\
ompute=0;\nif(($command=~/\\-recompute/))\n  {\n  \
  $recompute=1;\n    print \"-recompute \";\n  }\n\
\n$out=$msa;\n$out=~s/\\..*//;\n$out.=\".bph\";\ni\
f(($command=~/\\-o (\\S+)/))\n  {\n    $out=$1;\n \
   \n  }\nprint \"-out=$out \";\nif (-e $out && !$\
recompute)\n  {\n    print \"\\nNo Computation Req\
uired $out already exists\\n\";\n    &my_exit (EXI\
T_SUCCESS);\n    \n  }\n\n$n=100;\nif(($command=~/\
\\-n (\\d+)/))\n  {\n    $n=$1;\n  }\nprint \"-n=$\
n \";\n$seed=3;\nif(($command=~/\\-s (\\d+)/))\n  \
{\n    $seed=$1;\n  }\nprint \"-s=$seed \";\n\nif(\
($command=~/\\-run_name (\\d+)/))\n  {\n    $suffi\
x=$1;\n  }\nelse\n  {\n    $msa=~/([^.]+)/;\n    $\
suffix=$1;\n  }\nprint \"-run_name=$suffix\\n\";\n\
\n\nif ( $input eq \"seq\")\n  {\n    $seq=$msa;\n\
    $msa=\"$suffix.prot_msa\";\n    \n    if ($sty\
pe eq \"cdna\")\n      {\n	$cdna_seq=$seq;\n	$clea\
n_cdna_seq=&vtmpnam();\n	$seq=&vtmpnam();\n	`t_cof\
fee -other_pg seq_reformat -in $cdna_seq -action +\
clean_cdna >$clean_cdna_seq`;\n	`t_coffee -other_p\
g seq_reformat -in $clean_cdna_seq -action +transl\
ate >$seq`;\n	\n      }\n\n    if (!-e $msa || $re\
compute)\n      {\n	print \"\\n#####   Compute an \
MSA With $alignpg\\n\";\n	\n	if ( $alignpg eq \"t_\
coffee\")\n	  {`$alignpg $seq -outfile=$msa >/dev/\
null 2>/dev/null`;}\n	elsif ( $alignpg eq \"muscle\
\")\n	  {\n	    `$alignpg -in $seq > $msa 2>/dev/n\
ull`;\n	  }\n	elsif ( $alignpg eq \"clustalw\")\n	\
  {\n	    `$alignpg -infile=$seq -outfile=$msa -qu\
icktree >/dev/null 2>/dev/null`;\n	  }\n	elsif ( $\
align eq \"mafft\")\n	  {\n	    `$alignpg $seq > $\
msa >/dev/null 2>/dev/null`;\n	  }\n	else\n	  {\n	\
    `$alignpg -in=$seq -outfile=$msa`;\n	  }\n    \
  }\n    if (!-e $msa)\n      {\n	print \"\\nError\
: $alignpg Could Not produce the MSA $msa [FATAL]\\
\n\";\n      }\n\n    if ($stype eq \"cdna\")\n   \
   {\n	$msa2=\"$suffix.cdna_msa\";\n	`t_coffee -ot\
her_pg seq_reformat -in $clean_cdna_seq -in2 $msa \
-action +thread_dna_on_prot_aln -output fasta_aln \
 >$msa2`;\n	$msa=$msa2;\n      }\n    \n    $input\
=\"msa\";\n  }\n\n\n\n$seqboot_o=&vtmpnam();\n$seq\
boot_c=&vtmpnam();\n\n$protdist_o=&vtmpnam();\n$pr\
otdist_c=&vtmpnam();\nif ( $input eq \"msa\")\n  {\
\n    if ($tmode eq \"nj\" || $tmode eq \"upgma\")\
{$input=\"matrix\";}\n    \n    $lmsa= &vtmpnam ()\
;\n    `t_coffee -other_pg seq_reformat -in $msa -\
output phylip_aln > $lmsa`;\n    \n    if ( -e \"o\
utfile\"){unlink (\"outfile\");}\n    # run seqboo\
t\n  \n    if ( $n>1)\n      {\n	print \"Run SeqBo\
ot .....\";\n	open (F, \">$seqboot_c\");\n	print F\
 \"$lmsa\\nR\\n$n\\nY\\n$seed\\n\";\n	close (F);\n\
	`seqboot$exec_extension  < $seqboot_c`;\n	if ( -e\
 \"outfile\"){ print \"[OK]\\n\";}\n	else { print \
\"[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n	`mv ou\
tfile $seqboot_o`;\n      }\n    else\n      {\n	`\
cp $lmsa $seqboot_o`;\n      }\n\n    if ($rmsa){`\
cp $seqboot_o $rmsa`;}\n    \n    if ($tmode eq \"\
nj\" || $tmode eq \"upgma\")\n      {\n	if ( $styp\
e eq \"prot\")\n	  {\n	    # run protdist\n	    pr\
int \"Run Protdist [dmode=$dmode]\";\n	    if ($dm\
ode eq \"kimura\")\n	      {\n		$dmode=\"P\\nP\\nP\
\";\n	      }\n	    else\n	      {\n		print \"\\n$\
dmode is an unknown mode for Protdist [FATAL:msa2b\
ootstrap.pl]\\n\";\n		&my_exit (EXIT_FAILURE);\n	 \
     }\n	    open (F, \">$protdist_c\");\n	    if \
($n>1){print F \"$seqboot_o\\n$dmode\\nM\\nD\\n$n\\
\nY\\n\";}\n	    else {printf F \"$seqboot_o\\n$dm\
ode\\nY\\n\";}\n	    close (F);\n	    `protdist$ex\
ec_extension  < $protdist_c`;\n	    if ( -e \"outf\
ile\"){ print \"[OK]\\n\";}\n	    else { print \"[\
FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n	    `mv o\
utfile $protdist_o`;\n	 \n	  }\n	elsif ( $stype eq\
 \"cdna\" || $stype eq \"dna\")\n	  {\n	    print \
\"Run dnadist [dmode=default\";\n	    open (F, \">\
$protdist_c\");\n	    if ($n>1){print F \"$seqboot\
_o\\nM\\nD\\n$n\\nY\\n\";}\n	    else {printf F \"\
$seqboot_o\\nY\\n\";}\n	    close (F);\n	    `prot\
dist$exec_extension  < $protdist_c`;\n	    if ( -e\
 \"outfile\"){ print \"[OK]\\n\";}\n	    else { pr\
int \"[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n	  \
  `mv outfile $protdist_o`;\n	  }\n      }\n  }\ne\
lsif ( $input eq \"matrix\")\n  {\n    $protdist_o\
=&vtmpnam();\n    print \"MSA: $msa\\n\";\n    `cp\
 $msa $protdist_o`;\n    $n=1;\n  }\n\n\n\n\n\n$nb\
_o=&vtmpnam();\n$nb_c=&vtmpnam();\nif ($input eq \\
"matrix\" && $tmode ne \"parsimony\" && $tmode ne \
\"ml\")\n  {\n    print \"Run neighbor [tmode=$tmo\
de]\";\n\n    if ($tmode eq \"nj\")\n      {\n	$tm\
ode=\"\\nN\\nN\";\n      }\n    elsif ( $tmode eq \
\"upgma\")\n      {\n	$tmode = \"\\nN\";\n      }\\
n    else\n      {\n	print \"\\n ERROR: $tmode is \
an unknown tree computation mode\\n\";\n	&my_exit \
(EXIT_FAILURE);\n      }\n\n    open (F, \">$nb_c\\
");\n    if ($n>1){print F \"$protdist_o$tmode\\nM\
\\n$n\\n$seed\\nY\\n\";}\n    else {print F \"$pro\
tdist_o$tmode\\nY\\n\";}\n    close (F);\n\n    `n\
eighbor$exec_extension  < $nb_c`;\n    if ( -e \"o\
uttree\"){ print \"[Neighbor OK]\\n\";}\n    else \
{ print \"[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\\
n    `mv outtree $nb_o`;\n    unlink (\"outfile\")\
;\n  }\nelsif ($input eq \"msa\" && $tmode eq \"pa\
rsimony\")\n  {\n    if ( -e \"outfile\"){unlink (\
\"outfile\");}\n    if ( -e \"outtree\"){unlink (\\
"outtree\");}\n    \n    if ($stype eq \"prot\")\n\
      {\n	print \"Run protpars [tmode=$tmode]\";\n\
	open (F, \">$nb_c\");\n	if ($n>1){print F \"$seqb\
oot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\\n\";}\n	else \
{print F \"$seqboot_o\\nY\\n\";}\n	close (F);\n	`p\
rotpars$exec_extension  < $nb_c`;\n      }\n    el\
sif ( $stype eq \"dna\" || $stype eq \"cdna\")\n  \
    {\n	print \"Run dnapars [tmode=$tmode]\";\n	op\
en (F, \">$nb_c\");\n	if ($n>1){print F \"$seqboot\
_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\\n\";}\n	else {pr\
int F \"$seqboot_o\\nY\\n\";}\n	close (F);\n	`dnap\
ars$exec_extension  < $nb_c`;\n      }\n    if ( -\
e \"outtree\"){ print \"[OK]\\n\";}\n    else { pr\
int \"[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n   \
 `mv outtree $nb_o`;\n   unlink (\"outfile\");\n  \
}\nelsif ($input eq \"msa\" && $tmode eq \"ml\")\n\
  {\n    if ( -e \"outfile\"){unlink (\"outfile\")\
;}\n    if ( -e \"outtree\"){unlink (\"outtree\");\
}\n    \n    if ($stype eq \"prot\")\n      {\n	pr\
int \"Error: ML impossible with Protein Sequences \
[ERROR]\";\n	&my_exit (EXIT_FAILURE);\n      }\n  \
  elsif ( $stype eq \"dna\" || $stype eq \"cdna\")\
\n      {\n	print \"Run dnaml [tmode=$tmode]\";\n	\
open (F, \">$nb_c\");\n	if ($n>1){print F \"$seqbo\
ot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\\n\";}\n	else {\
print F \"$seqboot_o\\nY\\n\";}\n	close (F);\n	`dn\
aml$exec_extension  < $nb_c`;\n      }\n    if ( -\
e \"outtree\"){ print \"[OK]\\n\";}\n    else { pr\
int \"[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n   \
 `mv outtree $nb_o`;\n   unlink (\"outfile\");\n  \
}\n\n\nelse\n  {\n    `cp $msa $nb_o`;\n    $n=2;\\
n  }\n\nif ($rmsa && -e $seqboot_o){print \"\\nOut\
put List of $n Replicate MSA: $rmsa\\n\";`cp $seqb\
oot_o $rmsa`;}\nif ($rmat && -e $protdist_o){print\
 \"\\nOutput List of $n Replicate MATRICES: $rmat\\
\n\";`cp $protdist_o $rmat`;}\nif ($rtree && -e $n\
b_o){print \"\\nOutput List of $n Replicate TREES:\
 $rtree\\n\";`cp $nb_o $rtree`;}\n\n\n\n$con_o=&vt\
mpnam();\n$con_c=&vtmpnam();\nif ($n >1)\n  {\n   \
 print \"Run Consense.....\";\n    open (F, \">$co\
n_c\");\n    print F \"$nb_o\\nY\\n\";\n    close \
(F);\n    `consense$exec_extension  < $con_c`;\n  \
  if ( -s \"outtree\"  > 0) { print \"[OK]\\n\";}\\
n    else { print \"[FAILED]\\n\";&my_exit (EXIT_F\
AILURE);}\n    `mv outtree $con_o`;\n    unlink (\\
"outfile\");\n  }\nelse\n  {\n    `cp $nb_o $con_o\
`;\n  }\n\n\n`cp $con_o $out`;\nif ( !-e $out)\n  \
{\n    print \"Tree Computation failed [FAILED]\\n\
\";\n    &my_exit (EXIT_FAILURE);\n  }\nelsif ($n>\
1)\n  {\n    print \"\\nOutput Bootstrapped Tree: \
$out\\n\";\n    $avg=`t_coffee -other_pg seq_refor\
mat -in $out -action +avg_bootstrap`;\n    $avg=~s\
/\\n//g;\n    print \"$avg\\n\";\n  }\nelse\n  {\n\
    print \"\\nOutput Tree: $out\\n\";\n  }\n\nope\
n (F, \"$out\");\nwhile (<F>)\n  {\n    \n    $tre\
e.=$_;\n  }\nclose (F);\n$tree=~s/\\n//g;\nprint \\
"BPH: $tree\\n\";\n\n\n&my_exit (EXIT_SUCCESS);\n\\
nsub my_exit \n  {\n    my $m=@_[0];\n    &clean_v\
tmpnam();\n    exit ($m);\n  }\nsub vtmpnam \n  {\\
n    my $file;\n\n\n    $ntmp++;\n    $file=\"tmp4\
msa2bootstrap.$rseed.$$.$ntmp\";\n    \n    push (\
@tmpfile, $file);\n    return $file;\n  }\nsub cle\
an_vtmpnam \n  {\n    my $t;\n    foreach $t (@tmp\
file)\n      {\n	if ( -e $t){unlink ($t)};\n      \
}\n  }\n","use Env;\nuse FileHandle;\nuse Cwd;\nus\
e File::Path;\nuse Sys::Hostname;\n\n\nour $PIDCHI\
LD;\nour $ERROR_DONE;\nour @TMPFILE_LIST;\nour $EX\
IT_FAILURE=1;\nour $EXIT_SUCCESS=0;\n\nour $REFDIR\
=getcwd;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=\
1;\n\nour $PROGRAM=\"tc_generic_method.pl\";\nour \
$CL=$PROGRAM;\n\nour $CLEAN_EXIT_STARTED;\nour $de\
bug_lock=$ENV{\"DEBUG_LOCK\"};\nour $debug_generic\
_method=$ENV{\"DEBUG_GENERIC_METHOD\"};\nour $LOCK\
DIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif (!$LOCKDIR){$\
LOCKDIR=getcwd();}\nour $ERRORDIR=$ENV{\"ERRORDIR_\
4_TCOFFEE\"};\nour $ERRORFILE=$ENV{\"ERRORFILE_4_T\
COFFEE\"};\n&set_lock ($$);\nif (isshellpid(getppi\
d())){lock4tc(getppid(), \"LLOCK\", \"LSET\", \"$$\
\\n\");}\nour %RECODE;\nour $RECODE_N;\n\n\n\n\nou\
r $BLAST_MAX_NRUNS=2;\nour $COMMAND;\nour $PIDCHIL\
D;\n\n$REF_EMAIL=\"\";\n$tmp_dir=\"\";\n$init_dir=\
\"\";\n\n\n$test=0;\nif ($test==1)\n  {\n    $SERV\
ER=\"NCBI\";\n    $query=$ARGV[0];\n    $hitf=$ARG\
V[1];\n    %s=read_fasta_seq($query);\n    @sl=key\
s(%s);\n    &blast_xml2profile (\"xx\", $s{$sl[0]}\
{seq},$maxid,$minid,$mincov, $hitf);\n    myexit (\
$EXIT_FAILURE);\n  }\n\nforeach $v(@ARGV){$cl.=\"$\
v \";}\n$COMMAND=$cl;\n($mode)=&my_get_opt ( $cl, \
\"-mode=\",1,0);\n\n($A)=(&my_get_opt ( $cl, \"-na\
me1=\",0,0));\n($B)=(&my_get_opt ( $cl, \"-name2=\\
",0,0));\n($TMPDIR)=(&my_get_opt ( $cl, \"-tmpdir=\
\",0,0));\n($CACHE)=(&my_get_opt ( $cl, \"-cache=\\
",0,0));\n($SERVER)=((&my_get_opt ( $cl, \"-server\
=\",0,0)));\n($EMAIL)=((&my_get_opt ( $cl, \"-emai\
l=\",0,0)));\n\nif (!$A){$A=\"A\";}\nif (!$B){$B=\\
"B\";}\n\n\nif (!$TMPDIR)\n  {\n    $HOME=$ENV{HOM\
E};\n    if ($ENV{TMP_4_TCOFFEE}){$TMPDIR=$ENV{TMP\
_4_TCOFFEE};}\n    else{$TMPDIR=\"$HOME/.t_coffee/\
tmp/\";}\n  }\nif ( ! -d $TMPDIR)\n  {\n    mkdir \
$TMPDIR;\n  }\nif ( ! -d $TMPDIR)\n  {\n    print \
\"ERROR: Could not create temporary dir: $TMPDIR\\\
n\";\n    myexit ($EXIT_FAILURE);\n  }\n\n$EMAIL=~\
s/XEMAILX/\\@/g;\nif (!$EMAIL)\n  {\n    if ($ENV{\
EMAIL_4_TCOFFEE}){$EMAIL=$ENV{EMAIL_4_TCOFFEE};}\n\
    elsif ($ENV{EMAIL}){$EMAIL=$ENV{EMAIL};}\n    \
else {$EMAIL=$REF_EMAIL;}\n  }\n\n($maxid,$minid,$\
mincov,$trim)=(&my_get_opt ( $cl, \"-maxid=\",0,0,\
 \"-minid=\",0,0,\"-mincov=\",0,0, \"-trim=\",0,0)\
);\n\nif (!$cl=~/\\-maxid\\=/){$maxid=95;}\nif (!$\
cl=~/\\-minid\\=/){$minid=35;}\nif (!$cl=~/\\-minc\
ov\\=/){$mincov=80;}\nif (!$cl=~/\\-trim\\=/){$tri\
m;}\n\n\n\n\nif ($mode eq \"seq_msa\")\n  {\n    &\
seq2msa($mode,&my_get_opt ( $cl, \"-infile=\",1,1,\
 \"-method=\",1,2, \"-param=\",0,0,\"-outfile=\",1\
,0, \"-database=\",0,0));\n  }\nelsif ($mode eq \"\
blast2prf\")\n  {\n\n    blast2prf (&my_get_opt ( \
$cl, \"-infile=\",0,0,\"-seqfile=\",0,0,\"-outfile\
=\",0,0));\n  }\nelsif ( $mode eq \"tblastx_msa\")\
\n  {\n    &seq2tblastx_lib ($mode,&my_get_opt ( $\
cl, \"-infile=\",1,1, \"-outfile=\",1,0));\n  }\ne\
lsif ( $mode eq \"tblastpx_msa\")\n  {\n    &seq2t\
blastpx_lib ($mode,&my_get_opt ( $cl, \"-infile=\"\
,1,1, \"-outfile=\",1,0));\n  }\nelsif ( $mode eq \
\"thread_pair\")\n  {\n    &seq2thread_pair($mode,\
&my_get_opt ( $cl, \"-infile=\",1,1, \"-pdbfile1=\\
",1,1, \"-method=\",1,2,\"-param=\",0,0, \"-outfil\
e=\",1,0, ));\n  }\nelsif ( $mode eq \"pdbid_pair\\
")\n  {\n    &seq2pdbid_pair($mode,&my_get_opt ( $\
cl, \"-pdbfile1=\",1,0, \"-pdbfile2=\",1,0, \"-met\
hod=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0, ));\
\n  }\nelsif ( $mode eq \"pdb_pair\")\n  {\n    &s\
eq2pdb_pair($mode,&my_get_opt ( $cl, \"-pdbfile1=\\
",1,1, \"-pdbfile2=\",1,1, \"-method=\",1,2,\"-par\
am=\",0,0, \"-outfile=\",1,0, ));\n  }\nelsif ( $m\
ode eq \"rnapdb_pair\")\n{\n    &seq2rnapdb_pair($\
mode,&my_get_opt ( $cl, \"-pdbfile1=\",1,1, \"-pdb\
file2=\",1,1, \"-method=\",1,2,\"-param=\",0,0, \"\
-outfile=\",1,0, ));\n}\nelsif ( $mode eq \"profil\
e_pair\")\n  {\n     &seq2profile_pair($mode,&my_g\
et_opt ( $cl, \"-profile1=\",1,1, \"-profile2=\",1\
,1, \"-method=\",1,2,\"-param=\",0,0, \"-outfile=\\
",1,0 ));\n  }\nelsif ($mode eq \"pdb_template_tes\
t\")\n  {\n    &blast2pdb_template_test ($mode,&my\
_get_opt ( $cl, \"-infile=\",1,1));\n\n  }\nelsif \
($mode eq \"psi_template_test\")\n  {\n    &psibla\
st2profile_template_test ($mode,&my_get_opt ( $cl,\
 \"-seq=\",1,1,\"-blast=\",1,1));\n\n  }\n\nelsif \
( $mode eq \"pdb_template\")\n  {\n    &blast2pdb_\
template ($mode,&my_get_opt ( $cl, \"-infile=\",1,\
1, \"-database=\",1,0, \"-method=\",1,0, \"-outfil\
e=\",1,0,\"-pdb_type=\",1,0));\n  }\n\nelsif ( $mo\
de eq \"profile_template\")\n  {\n\n    &seq2profi\
le_template ($mode,&my_get_opt ( $cl, \"-infile=\"\
,1,1, \"-database=\",1,0, \"-method=\",1,0, \"-out\
file=\",1,0));\n  }\nelsif ( $mode eq \"psiprofile\
_template\")\n  {\n    &seq2profile_template ($mod\
e,&my_get_opt ( $cl, \"-infile=\",1,1, \"-database\
=\",1,0, \"-method=\",1,0, \"-outfile=\",1,0));\n \
 }\nelsif ( $mode eq \"RNA_template\")\n  {\n    &\
seq2RNA_template ($mode,&my_get_opt ( $cl, \"-infi\
le=\",1,1,\"-pdbfile=\",1,1,\"-outfile=\",1,0));\n\
  }\nelsif ( $mode eq \"tm_template\")\n  {\n    &\
seq2tm_template ($mode,&my_get_opt ( $cl, \"-infil\
e=\",1,1,\"-arch=\",1,1,\"-psv=\",1,1, \"-outfile=\
\",1,0));\n  }\nelsif ( $mode eq \"psitm_template\\
")\n  {\n    &seq2tm_template ($mode,&my_get_opt (\
 $cl, \"-infile=\",1,1, \"-arch=\",1,1,\"-psv=\",1\
,1, \"-outfile=\",1,0,\"-database=\",1,0));\n  }\n\
elsif ( $mode eq \"ssp_template\")\n  {\n    &seq2\
ssp_template ($mode,&my_get_opt ( $cl, \"-infile=\\
",1,1,\"-seq=\",1,1,\"-obs=\",1,1, \"-outfile=\",1\
,0));\n  }\nelsif ( $mode eq \"psissp_template\")\\
n  {\n    &seq2ssp_template ($mode,&my_get_opt ( $\
cl, \"-infile=\",1,1,\"-seq=\",1,1,\"-obs=\",1,1, \
\"-outfile=\",1,0));\n  }\n\n\n\nelse\n  {\n    my\
exit(flush_error( \"$mode is an unknown mode of tc\
_generic_method.pl\"));\n  }\nmyexit ($EXIT_SUCCES\
S);\n\n\nsub seq2ssp_template\n  {\n  my ($mode, $\
infile,$gor_seq,$gor_obs,$outfile)=@_;\n  my %s, %\
h;\n  my $result;\n  my (@profiles);\n  &set_tempo\
rary_dir (\"set\",$infile,\"seq.pep\");\n  %s=read\
_fasta_seq (\"seq.pep\");\n\n\n  open (R, \">resul\
t.aln\");\n\n  #print stdout \"\\n\";\n  foreach $\
seq (keys(%s))\n    {\n\n      open (F, \">seqfile\
\");\n      $s{$seq}{seq}=uc$s{$seq}{seq};\n      \
print (F \">$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\\
n      close (F);\n      $lib_name=\"$s{$seq}{name\
}.ssp\";\n      $lib_name=&clean_file_name ($lib_n\
ame);\n\n      if ($mode eq \"ssp_template\"){&seq\
2gor_prediction ($s{$seq}{name},$s{$seq}{seq}, \"s\
eqfile\", $lib_name,$gor_seq, $gor_obs);}\n      e\
lsif ($mode eq \"psissp_template\")\n	{\n	  &seq2m\
sa_gor_prediction ($s{$seq}{name},$s{$seq}{seq},\"\
seqfile\", $lib_name,$gor_seq, $gor_obs);\n	}\n\n \
     if ( !-e $lib_name)\n	{\n	  myexit(flush_erro\
r(\"GORIV failed to compute the secondary structur\
e of $s{$seq}{name}\"));\n	  myexit ($EXIT_FAILURE\
);\n	}\n      else\n	{\n	  print stdout \"!\\tProc\
ess: >$s{$seq}{name} _E_ $lib_name \\n\";\n	  prin\
t R \">$s{$seq}{name} _E_ $lib_name\\n\";\n	}\n   \
   unshift (@profiles, $lib_name);\n    }\n  close\
 (R);\n  &set_temporary_dir (\"unset\",$mode, $met\
hod,\"result.aln\",$outfile, @profiles);\n}\n\nsub\
 seq2tm_template\n  {\n  my ($mode,$infile,$arch,$\
psv,$outfile,$db)=@_;\n  my %s, %h;\n  my $result;\
\n  my (@profiles);\n  &set_temporary_dir (\"set\"\
,$infile,\"seq.pep\");\n  %s=read_fasta_seq (\"seq\
.pep\");\n\n\n  open (R, \">result.aln\");\n\n  #p\
rint stdout \"\\n\";\n  foreach $seq (keys(%s))\n \
   {\n      open (F, \">seqfile\");\n      print (\
F \">$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\n      \
close (F);\n      $lib_name=\"$s{$seq}{name}.tmp\"\
;\n      $lib_name=&clean_file_name ($lib_name);\n\
\n      if ($mode eq \"tm_template\")\n	{\n	  &saf\
e_system (\"t_coffee -other_pg fasta_seq2hmmtop_fa\
sta.pl -in=seqfile -out=$lib_name -arch=$arch -psv\
=$psv\");\n	}\n      elsif ( $mode eq \"psitm_temp\
late\")\n	{\n	  &seq2msa_tm_prediction ($s{$seq}{n\
ame},$s{$seq}{seq}, $db, \"seqfile\", $lib_name,$a\
rch, $psv);\n	}\n      if ( !-e $lib_name)\n	{\n	 \
 myexit(flush_error(\"hmmtop failed to compute the\
 secondary structure of $s{$seq}{name}\"));\n	  my\
exit ($EXIT_FAILURE);\n	}\n      else\n	{\n	  prin\
t stdout \"!\\tProcess: >$s{$seq}{name} _T_ $lib_n\
ame\\n\";\n	  print R \">$s{$seq}{name} _T_ $lib_n\
ame\\n\";\n	}\n      unshift (@profiles, $lib_name\
);\n    }\n  close (R);\n  &set_temporary_dir (\"u\
nset\",$mode, $method,\"result.aln\",$outfile, @pr\
ofiles);\n}\n\n\n\nsub seq2RNA_template\n  {\n    \
\n    my ($mode, $infile, $pdbfile, $outfile)=@_;\\
n    my %s, %h ;\n    my $result;\n    my (@profil\
es);\n    my ($seq_mode, $pdb_mode, $pwd);\n    \n\
    #use $seq_mode to estimate the template of seq\
uences WITHOUT a PDB\n    #use $pdb_mode to estima\
te the template of sequences WITH    a PDB\n\n    \
$seq_mode=$ENV{\"SEQ2TEMPLATE4_F_\"};\n    $pdb_mo\
de=$ENV{\"PDB2TEMPLATE4_F_\"};\n    \n    if (!$pd\
b_mode){$pdb_mode=\"find_pair-p\";}\n    if (!$seq\
_mode){$seq_mode=\"RNAplfold\";}\n    \n    my $cw\
d = cwd();\n    &set_temporary_dir (\"set\",$infil\
e,\"seq.pep\");\n    %s=read_fasta_seq (\"seq.pep\\
");\n    %pdb_template_h = &read_template_file($pd\
bfile);\n    my $pdb_chain;\n    \n       \n    op\
en (R, \">result.aln\");\n    #print stdout \"\\n\\
";\n    foreach $seq (keys(%s))\n      {\n	\n	open\
 (F, \">seqfile\");\n	print (F \">$s{$seq}{name}\\\
n$s{$seq}{seq}\\n\");\n	close (F);\n	$pdb_chain = \
$pdb_template_h{$seq};\n	$lib_name=\"$s{$seq}{name\
}.rfold\";\n	$lib_name=&clean_file_name ($lib_name\
);\n	if ($pdb_template_h{$seq} eq \"\")\n	  {\n	  \
  if    ($seq_mode eq \"RNAplfold\"){RNAplfold2lib\
 (\"seqfile\", \"$lib_name\");}\n	    elsif ($seq_\
mode eq \"no\"){$lib_name=0;}\n	    else\n	      {\
\n		myexit(add_error (EXIT_FAILURE,$$,$$,getppid()\
, \"seq2RNA_template failure::method $seq_mode not\
 available for sequences without PDB structures\")\
);\n	      }\n	  }\n	elsif ($pdb_template_h{$seq} \
ne \"\")\n	  {\n	    my $pdbf;\n	    if    ( -e \"\
$cwd/$pdb_chain\"   ){$pdbf=\"$cwd/$pdb_chain\"; }\
\n	    elsif ( -e  $pdb_chain         ){$pdbf=\"$p\
db_chain\";      }\n	    elsif ( -e  \"$CACHE$pdb_\
chain\" ){$pdbf=\"$CACHE$pdb_chain\";}\n	    elsif\
 ( -e  \"$CACHE/$pdb_chain\"){$pdbf=\"$CACHE/$pdb_\
chain\";}\n	    else\n	      {\n		myexit(flush_err\
or(\"Could not read $pdb_chain \"));\n	      }\n\n\
	    if($pdb_mode eq \"x3dna-ssr\")\n	      {\n		x\
3dnassr2lib (\"seqfile\", \"$pdbf\", \"$lib_name\"\
);\n	      }\n	    elsif ($pdb_mode eq \"find_pair\
-p\")\n	      {\n		x3dna_find_pair2lib (\"seqfile\\
", \"$pdbf\", \"$lib_name\", \"find_pair -p\");\n	\
      }\n	    elsif ($pdb_mode eq \"find_pair\")\n\
	      {\n		x3dna_find_pair2lib (\"seqfile\", \"$p\
dbf\", \"$lib_name\", \"find_pair\");\n	      }\n	\
    elsif ($pdb_mode eq \"RNAplfold\")\n	      {\n\
		RNAplfold2lib (\"seqfile\", \"$lib_name\");\n	  \
    }\n	    elsif ($pdb_mode eq \"no\"){$lib_name=\
0;}\n	    else\n	      {\n		myexit(add_error (EXIT\
_FAILURE,$$,$$,getppid(), \"seq2RNA_template failu\
re::Could not find method $pdb_mode\"));\n	      }\
\n	  }\n	if ($lib_name)\n	  {\n	    print stdout \\
"!\\tProcess: >$s{$seq}{name} _F_ $lib_name\\n\";\\
n	    print R \">$s{$seq}{name} _F_ $lib_name\\n\"\
;\n	    unshift (@profiles, $lib_name);\n	  }\n   \
   }\n    close (R);\n    &set_temporary_dir (\"un\
set\",$mode, $method,\"result.aln\",$outfile, @pro\
files);\n  }\n\n\n\nsub psiblast2profile_template_\
test\n  {\n  my ($mode, $seq,$blast)=@_;\n  my %s,\
 %h, ;\n  my ($result,$psiblast_output,$profile_na\
me,@profiles);\n  my $trim=0;\n  my $maxid=100;\n \
 my $minid=0;\n  my $mincov=0;\n  my $maxcov=100;\\
n\n  %s=read_fasta_seq ($seq);\n  open (R, \">resu\
lt.aln\");\n\n  #print stdout \"\\n\";\n  foreach \
$seq (keys(%s))\n    {\n\n      open (F, \">seqfil\
e\");\n      print (F \">$A\\n$s{$seq}{seq}\\n\");\
\n      close (F);\n      $psiblast_output=$blast;\
\n      if ( -e $psiblast_output)\n	{\n	  %profile\
=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$\
maxid, $minid,$mincov,$psiblast_output);\n\n\n\n	 \
 $profile_name=\"$s{$seq}{name}.prf\";\n	  $profil\
e_name=&clean_file_name ($profile_name);\n	  unshi\
ft (@profiles, $profile_name);\n	  output_profile \
($profile_name, \\%profile, $trim);\n	  print stdo\
ut \"!\\tProcess: >$s{$seq}{name} _R_ $profile_nam\
e [$profile{n} Seq.] [$SERVER/blast/$db][$CACHE_ST\
ATUS]\\n\";\n	  print R \">$s{$seq}{name} _R_ $pro\
file_name\\n\";\n	}\n    }\n  close (R);\n\n  die;\
\n}\nsub seq2profile_template\n    {\n      my ($m\
ode, $infile, $db, $method, $outfile)=@_;\n      i\
f    ($method eq \"psiblast\"){return psiblast2pro\
file_template ($mode, $infile, $db, $method, $outf\
ile);}\n      elsif ($method eq \"blastp\")   {ret\
urn psiblast2profile_template ($mode, $infile, $db\
, $method, $outfile);}\n      elsif ($method eq \"\
hh\")      {return hh2profile_template ($mode, $in\
file, $db, $method, $outfile);}\n    }\n\nsub psib\
last2profile_template\n  {\n  my ($mode, $infile, \
$db, $method, $outfile)=@_;\n  my %s, %h, ;\n  my \
($result,$psiblast_output,$profile_name,@profiles)\
;\n  &set_temporary_dir (\"set\",$infile,\"seq.pep\
\");\n  %s=read_fasta_seq (\"seq.pep\");\n  open (\
R, \">result.aln\");\n\n  #print stdout \"\\n\";\n\
  foreach $seq (keys(%s))\n    {\n      open (F, \\
">seqfile\");\n      print (F \">$A\\n$s{$seq}{seq\
}\\n\");\n      close (F);\n      $psiblast_output\
=&run_blast ($s{$seq}{name},$method, $db, \"seqfil\
e\",\"outfile\");\n\n      if ( -e $psiblast_outpu\
t)\n	{\n	  my %profile=blast_xml2profile($s{$seq}{\
name}, $s{$seq}{seq},$maxid, $minid,$mincov,$psibl\
ast_output);\n	  unlink ($psiblast_output);\n	  \n\
	  $profile_name=\"$s{$seq}{name}.prf\";\n	  $prof\
ile_name=&clean_file_name ($profile_name);\n	  uns\
hift (@profiles, $profile_name);\n	  output_profil\
e ($profile_name, \\%profile, $trim);\n	  \n	  pri\
nt stdout \"!\\tProcess: >$s{$seq}{name} _R_ $prof\
ile_name [$profile{n} Seq.] [$SERVER/blast/$db][$C\
ACHE_STATUS]\\n\";\n	  print R \">$s{$seq}{name} _\
R_ $profile_name\\n\";\n	  \n	  \n	}\n      \n    \
}\n  close (R);\n  \n  \n\n  &set_temporary_dir (\\
"unset\",$mode, $method,\"result.aln\",$outfile, @\
profiles);\n}\n\nsub hh2profile_template\n  {\n\n \
 #for each sequence, build a profile, in FASTA, wi\
th ungapped querry on top  \n  my ($mode, $infile,\
 $db, $method, $outfile)=@_;\n  my %s, %h, ;\n  my\
 ($result,$psiblast_output,$profile_name,@profiles\
);\n  &set_temporary_dir (\"set\",$infile,\"seq.pe\
p\");\n  %s=read_fasta_seq (\"seq.pep\");\n  open \
(R, \">result.aln\");\n  \n  my $hh=$ENV{\"HHSEARC\
H_4_TCOFFEE\"};\n  if (!$hh)\n    {\n      print \\
"ERROR: HHSEARCH_4_TCOFFEE is not set\\n\";\n     \
 myexit ($EXIT_FAILURE);\n    }\n  \n  #print stdo\
ut \"\\n\";\n  foreach $seq (keys(%s))\n    {\n   \
   my ($profile_name, $nseq);\n      open (F, \">s\
eqfile\");\n      print (F \">$A\\n$s{$seq}{seq}\\\
n\");\n      close (F);\n      \n      #This funct\
ion should input a querry and a database and retur\
n as output a fasta MSA with quesry on top\n      \
$profile_name=\"$s{$seq}{name}.prf\";\n      $prof\
ile_name=&clean_file_name ($profile_name);\n      \
unshift (@profiles, $profile_name);\n      \n     \
 \n      safe_system  (\"$hh -name=$s{$seq}{name} \
-method=search -db=$db -seq=seqfile -outfile=$prof\
ile_name\");\n      if (-e $profile_name){$nseq=fa\
sta2nseq($profile_name);}\n      \n      print std\
out \"!\\tProcess: >$s{$seq}{name} _R_ $profile_na\
me [$nseq Seq.] [$method/$db][$CACHE_STATUS]\\n\";\
\n      print R \">$s{$seq}{name} _R_ $profile_nam\
e\\n\";\n    }\n  close (R);\n  &set_temporary_dir\
 (\"unset\",$mode, $method,\"result.aln\",$outfile\
, @profiles);\n}\n\nsub blast2pdb_template_test\n \
   {\n      my ($mode,$infile)=@_;\n      my ($max\
id,$minid,$mincov);\n      $maxid=100;\n      $min\
id=0;\n      $mincov=0;\n\n      print \"$infile\\\
n\";\n\n      %p=blast_xml2profile($s{$seq}{name},\
 $s{$seq}{seq},$maxid, $minid,$mincov,$infile);\n \
     $c=1;\n      print stdout \"!\\tProcess: >$s{\
$seq}{name} [$SERVER/blast/$db][$CACHE_STATUS]\\n\\
";\n      while (!$found && $c<$p{n})\n	{\n	  $pdb\
id=&id2pdbid($p{$c}{identifyer});\n	  if ( length \
($pdbid)>5){$pdbid=id2pdbid($p{$c}{definition});}\\
n\n	  if ( length ($pdbid)>5)\n	    {\n	      myex\
it(add_error (EXIT_FAILURE,$$,$$,getppid(), \"BLAS\
T_FAILURE::Could Not Parse PDBID ($p{$c}{identifye\
r},$p{$c}{definition})\"));\n	    }\n\n\n	  if (!&\
pdb_is_released($pdbid))\n	    {\n	      print std\
out \"\\t\\t**$pdbid [WARNIG: PDB NOT RELEASED or \
WITHDRAWN]\\n\";\n	      $c++;\n	    }\n	  elsif (\
!&pdb_has_right_type ($pdbid,$type))\n	    {\n	   \
   my $ptype=&pdb2type ($pdbid);\n	      my $etype\
=&type2etype($type);\n\n	      print stdout \"\\t\\
\t**$pdbid [$ptype cannot be used (expected: $etyp\
e)]\\n\";\n	      $c++;\n	    }\n	  else\n	    {\n\
	      $found=1;\n	    }\n	}\n\n      if ($found)\\
n	{\n	  print stdout \"\\t\\t >$s{$seq}{name} _P_ \
$pdbid\\n\";\n	}\n      else\n	{\n	  print stdout \
\"\\t\\t >$s{$seq}{name} No Template Selected\\n\"\
;\n	}\n      die;\n    }\nsub blast2pdb_template\n\
  {\n  my ($mode, $infile, $db, $method, $outfile,\
$type)=@_;\n  my %s, %h, ;\n  my ($result,$blast_o\
utput);\n  &set_temporary_dir (\"set\",$infile,\"s\
eq.pep\");\n  %s=read_fasta_seq (\"seq.pep\");\n  \
open (R, \">result.aln\");\n\n\n  #print stdout \"\
\\n\";\n  foreach $seq (keys(%s))\n    {\n      my\
 $c;\n      my $found;\n\n      open (F, \">seqfil\
e\");\n      print (F \">$A\\n$s{$seq}{seq}\\n\");\
\n      close (F);\n\n      $blast_output=&run_bla\
st ($s{$seq}{name},$method, $db, \"seqfile\",\"out\
file\");\n\n      %p=blast_xml2profile($s{$seq}{na\
me}, $s{$seq}{seq},$maxid, $minid,$mincov,$blast_o\
utput);\n      unlink ($blast_output);\n\n      $c\
=1;\n      print stdout \"!\\tProcess: >$s{$seq}{n\
ame} [$SERVER/blast/$db][$CACHE_STATUS]\\n\";\n   \
   while (!$found && $c<$p{n})\n	{\n	  $pdbid=&id2\
pdbid($p{$c}{identifyer});\n	  if ( length ($pdbid\
)>5){$pdbid=id2pdbid($p{$c}{definition});}\n\n	  i\
f ( length ($pdbid)>5)\n	    {\n	      myexit(add_\
error (EXIT_FAILURE,$$,$$,getppid(), \"BLAST_FAILU\
RE::Could Not Parse PDBID ($p{$c}{identifyer},$p{$\
c}{definition})\"));\n	    }\n\n\n	  if (!&pdb_is_\
released($pdbid))\n	    {\n	      print stdout \"\\
\t\\t**$pdbid [PDB NOT RELEASED or WITHDRAWN]\\n\"\
;\n	      $c++;\n	    }\n	  elsif (!&pdb_has_right\
_type ($pdbid,$type))\n	    {\n	      my $ptype=&p\
db2type ($pdbid);\n	      my $etype=&type2etype($t\
ype);\n\n	      print stdout \"\\t\\t**$pdbid [$pt\
ype cannot be used (expected: $etype)]\\n\";\n	   \
   $c++;\n	    }\n	  else\n	    {\n	      $found=1\
;\n	    }\n	}\n\n      if ($found)\n	{\n	  print R\
 \">$s{$seq}{name} _P_ $pdbid\\n\";\n	  print stdo\
ut \"\\t\\t >$s{$seq}{name} _P_ $pdbid\\n\";\n	}\n\
      else\n	{\n	  print R \">$s{$seq}{name}\\n\";\
\n	  print stdout \"\\t\\t >$s{$seq}{name} No Temp\
late Selected\\n\";\n	}\n    }\n  close (R);\n  &s\
et_temporary_dir (\"unset\",$mode, $method,\"resul\
t.aln\",$outfile);\n}\nsub type2etype\n  {\n    my\
 $type=shift;\n    my $etype;\n\n    if ( $type=~/\
n/){$etype.=\"NMR \";}\n    if ( $type=~/d/){$etyp\
e.=\"diffraction \";}\n    if ( $type=~/m/){$etype\
.=\"model \";}\n    return $etype;\n  }\nsub pdb2t\
ype\n  {\n     my $pdb=shift;\n     my $f=vtmpnam(\
);\n\n     my $value= &safe_system (\"t_coffee -ot\
her_pg extract_from_pdb -model_type $pdb > $f\");\\
n     my $r=&file2string ($f);\n     chomp($r);\n \
    return $r;\n   }\nsub pdb_has_right_type\n  {\\
n    my $pdb=shift;\n    my $type=shift;\n\n    my\
 $f=vtmpnam();\n\n    my $value= &safe_system (\"t\
_coffee -other_pg extract_from_pdb -model_type $pd\
b > $f\");\n    my $r=&file2string ($f);\n    chom\
p($r);\n\n\n    if ( $r eq \"NMR\" && $type=~/n/){\
return 1;}\n    elsif ( $r eq \"diffraction\" && $\
type=~/d/){return 1;}\n    elsif ( $r eq \"model\"\
 && $type=~/m/){return 1;}\n    else {return 0;}\n\
  }\nsub pdb_is_released\n  {\n    my $pdb=shift;\\
n    my $f=vtmpnam();\n\n    $value= &safe_system \
(\"t_coffee -other_pg extract_from_pdb -is_release\
d_pdb_name $pdb > $f\");\n    my $r=&file2string (\
$f);\n    chomp($r);\n    return $r;\n  }\nsub bla\
st_msa\n  {\n    my ($blast,$infile,$db,$outfile)=\
@_;\n    my ($a, %s1, %s, %qs, %qs1);\n    my $seq\
file;\n    my $SEQ=new FileHandle;\n    my $seqfil\
e=\"seqfile\";\n    my @txt;\n\n\n    %s1=&read_fa\
sta_seq ($db);\n    %s=&fasta_hash2index_hash(%s1)\
;\n    %qs1=&read_fasta_seq ($infile);\n    %qs=&f\
asta_hash2index_hash(%qs1);\n\n\n    #&safe_system\
 (\"formatdb -i $db\");\n    if ($blast eq \"blast\
p\"){&safe_system  (\"blastall -i $infile -d $db -\
m7 -p blastp -o io\");}\n    elsif ($blast eq \"bl\
astn\"){&safe_system  (\"blastn -query $infile -db\
 $db -outfmt 5 -word_size 4 -out io\");}\n\n    &s\
et_blast_type (\"io\");\n\n\n    my %FB=&xml2tag_l\
ist (\"io\", \"Iteration\");\n    open (F, \">$out\
file\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\n\
    print F \"$s{n}\\n\";\n    for ( my $a=0; $a<$\
s{n}; $a++)\n      {\n	print F \"$s{$a}{name} $s{$\
a}{len} $s{$a}{seq}\\n\";\n      }\n\n\n    for ( \
my $a=0; $a<$FB{n}; $a++)\n      {\n	my %p=blast_x\
ml2profile ($qs{$a}{name}, $qs{$a}{seq},100, 0, 0,\
 $FB{$a}{body});\n	my $query=$p{0}{name};\n	my $i=\
 $s1{$query}{order}+1;\n	for (my $b=1; $b<$p{n}; $\
b++)\n	  {\n	    my $l=length ($p{$b}{Qseq});\n	  \
  my $hit=$p{$b}{definition};\n	    my $Qstart=$p{\
$b}{Qstart};\n	    my $Hstart=$p{$b}{Hstart};\n	  \
  my $identity=$p{$b}{identity};\n	    my @lrQ=spl\
it (//,$p{$b}{Qseq});\n	    my @lrH=split (//,$p{$\
b}{Hseq});\n\n	    my $j= $s1{$hit}{order}+1;\n	  \
  #if ( $j==$i){next;}\n	    printf F \"# %d %d\\n\
\", $i, $j;\n	    #  print  F \"\\n$p{$b}{Qseq} ($\
Qstart)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for ($c=\
0; $c<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n		m\
y $rH=$lrH[$c];\n		my $n=0;\n\n		if ($rQ ne \"-\")\
{$n++, $Qstart++;}\n		if ($rH ne \"-\"){$n++; $Hst\
art++;}\n\n		if ( $n==2)\n		  {\n		    printf F \"\
\\t%d %d %d\\n\", $Qstart-1, $Hstart-1,$identity;\\
n		  }\n	      }\n	  }\n      }\n    print F \"! S\
EQ_1_TO_N\\n\";\n    close (F);\n    return $outpu\
t;\n  }\n\nsub blast_msa_old\n  {\n    my ($infile\
,$outfile)=@_;\n    my ($a, %seq);\n    %s1=&read_\
fasta_seq ($infile);\n    foreach $s (keys (%s1))\\
n      {\n	$i=$s1{$s}{order};\n	$s{$i}{name}=$s;\n\
	$s{$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len}=length( $\
s{$i}{seq});\n	$s{n}++;\n      }\n    &safe_system\
 (\"formatdb -i $infile\");\n    &safe_system (\"b\
lastall -i $infile -d $infile -m7 -o io\");\n    &\
set_blast_type (\"io\");\n\n    %FB=&xml2tag_list \
(\"io\", \"Iteration\");\n\n    open (F, \">$outfi\
le\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\n  \
  print F \"$s{n}\\n\";\n    for ( $a=0; $a<$s{n};\
 $a++)\n      {\n	print F \"$s{$a}{name} $s{$a}{le\
n} $s{$a}{seq}\\n\";\n      }\n    for ( $a=0; $a<\
$FB{n}; $a++)\n      {\n	%p=blast_xml2profile ($s{\
$a}{name}, $s{$a}{seq},100, 0, 0, $FB{$a}{body});\\
n	for ($b=1; $b<$p{n}; $b++)\n	  {\n	    my $l=len\
gth ($p{$b}{Qseq});\n	    my $hit=$p{$b}{definitio\
n};\n	    my $Qstart=$p{$b}{Qstart};\n	    my $Hst\
art=$p{$b}{Hstart};\n	    my $identity=$p{$b}{iden\
tity};\n	    my @lrQ=split (//,$p{$b}{Qseq});\n	  \
  my @lrH=split (//,$p{$b}{Hseq});\n	    my $i= $s\
1{$s{$a}{name}}{order}+1;\n	    my $j= $s1{$hit}{o\
rder}+1;\n	    #if ( $j==$i){next;}\n	    printf F\
 \"# %d %d\\n\", $i, $j;\n	    #  print  F \"\\n$p\
{$b}{Qseq} ($Qstart)\\n$p{$b}{Hseq} ($Hstart)\";\n\
	    for ($c=0; $c<$l; $c++)\n	      {\n		my $rQ=$\
lrQ[$c];\n		my $rH=$lrH[$c];\n		my $n=0;\n\n		if (\
$rQ ne \"-\"){$n++, $Qstart++;}\n		if ($rH ne \"-\\
"){$n++; $Hstart++;}\n\n		if ( $n==2)\n		  {\n		  \
  printf F \"\\t%d %d %d\\n\", $Qstart-1, $Hstart-\
1,$identity;\n		  }\n	      }\n	  }\n      }\n    \
print F \"! SEQ_1_TO_N\\n\";\n    close (F);\n    \
return $output;\n\n  }\n\nsub seq2msa\n  {\n    my\
 ($mode, $infile, $method, $param, $outfile,$datab\
ase)=@_;\n    &set_temporary_dir (\"set\",$infile,\
\"seq.pep\", $database, \"db.pep\");\n    $param.=\
\" >/dev/null 2>&1 \";\n\n\n    #make sure test.pe\
p is in FASTA\n    &safe_system (\"t_coffee -other\
_pg seq_reformat -in seq.pep -output fasta_seq > x\
\");\n    `mv x seq.pep`;\n\n    if ( $method eq \\
"blastp\")\n      {\n	&blast_msa (\"blastp\",\"seq\
.pep\",$database,\"result.aln\");\n      }\n    el\
sif ( $method eq \"blastn\")\n      {\n	&blast_msa\
 (\"blastn\",\"seq.pep\",$database,\"result.aln\")\
;\n      }\n\n    elsif ( $method eq \"muscle\")\n\
      {\n	`muscle -in seq.pep -out result.aln $par\
am`;\n      }\n    elsif ( $method eq \"probcons\"\
)\n      {\n	`probcons seq.pep >result.aln 2>/dev/\
null`;\n      }\n    elsif ( $method eq \"mafft\")\
\n      {\n	`mafft --quiet --localpair --maxiterat\
e 1000 seq.pep> result.aln  2>/dev/null`\n      }\\
n    elsif ( $method=~/prank/)\n      {\n	`$method\
 -d=seq.pep -o=result.aln -quiet 2>/dev/null`;\n	`\
mv result.aln.1.fas result.aln`;\n      }\n    els\
if ($method eq \"clustalo\")\n      {\n	`clustalo \
-i seq.pep > result.aln`;\n      }\n\n    else\n  \
    {\n	`$method -infile=seq.pep -outfile=result.a\
ln`;\n      }\n\n    &set_temporary_dir (\"unset\"\
,$mode, $method,\"result.aln\",$outfile);\n    mye\
xit ($EXIT_SUCCESS);\n  }\n\nsub seq2thread_pair\n\
  {\n    my ($mode, $infile, $pdbfile1, $method, $\
param, $outfile)=@_;\n    &set_temporary_dir (\"se\
t\",$infile,\"seq.pep\",$pdbfile1,\"struc.pdb\");\\
n    if ($method eq \"fugueali\")\n      {\n	#Env \
Variable that need to be defined for Fugue\n	if (!\
$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_LIB_LIST}=\"DUMMY\
\";}\n	if (!$ENV{HOMSTRAD_PATH})  {$ENV{HOMSTRAD_P\
ATH}=\"DUMMY\";}\n	if (!$ENV{HOMS_PATH}){$ENV{HOMS\
_PATH}=\"DUMMY\";}\n\n	`joy struc.pdb >x 2>x`;\n	&\
check_file(\"struc.tem\", \"Joy failed [FATAL:$PRO\
GRAM/$method]\");\n	`melody -t struc.tem >x 2>x`;\\
n	&check_file(\"struc.tem\", \"Melody failed [FATA\
L:$PROGRAM/$method]\");\n	`fugueali -seq seq.pep -\
prf struc.fug -print > tmp_result.aln`;\n\n	&check\
_file(\"tmp_result.aln\", \"Fugue failed [FATAL:$P\
ROGRAM/$method]\");\n	&safe_system (\"t_coffee -ot\
her_pg seq_reformat -in tmp_result.aln -output fas\
ta_aln >result.aln\");\n      }\n    elsif ( $meth\
od eq \"t_coffee\")\n      {\n	&safe_system (\"t_c\
offee -in Pstruc.pdb Sseq.pep Mslow_pair -outfile \
result.aln -quiet\");\n      }\n    else\n      {\\
n	&safe_system (\"$method -infile=seq.pep -pdbfile\
1=struc.pdb -outfile=result.aln $param>x 2>x\");\n\
      }\n    &set_temporary_dir (\"unset\",$mode,$\
method,\"result.aln\",$outfile);\n    myexit ($EXI\
T_SUCCESS);\n  }\nsub seq2pdbid_pair\n  {\n    my \
($mode, $pdbfile1, $pdbfile2, $method, $param, $ou\
tfile)=@_;\n    my ($name);\n\n\n    &set_temporar\
y_dir (\"set\");\n    $name=$pdbfile1.\" \".$pdbfi\
le2;\n\n    if (    &cache_file(\"GET\",\"\",\"$na\
me\",\"$method\",\"dali\",$outfile,\"EBI\"))\n    \
  {return $outfile;}\n    else\n      {\n	if ($met\
hod eq \"daliweb\")\n	  {\n	    $pdbfile1=~/(....)\
(.)/;\n	    $id1=$1; $c1=$2;\n\n	    $pdbfile2=~/(\
....)(.)/;\n	    $id2=$1; $c2=$2;\n\n	    $command\
=\"t_coffee -other_pg dalilite.pl --pdb1 $id1 --ch\
ainid1 $c1 --pdb2 $id2 --chainid2 $c2 --email=$EMA\
IL  >dali_stderr 2>dali_stderr\";\n	    $dali=`$co\
mmand`;\n\n	    open (F, \"dali_stderr\");\n	    w\
hile (<F>)\n	      {\n		if ( /JobId: dalilite-(\\S\
+)/)\n		{\n		  $jobid=$1;\n		}\n	      }\n	    clo\
se (F);\n	    unlink (\"dali_stderr\");\n\n	    $o\
utput1=\"dalilite-$jobid.txt\";\n	    if ( -e $out\
put1)\n	      {\n		unlink ($output1);\n		&url2file\
 (\"http://www.ebi.ac.uk/Tools/es/cgi-bin/jobresul\
ts.cgi/dalilite/dalilite-$jobid/aln.html\", \"outp\
ut2\");\n\n		if ( -e \"output2\")\n		  {\n		    my\
 ($seq1, $seq2);\n		    $seq1=$seq2=\"\";\n\n		   \
 open (F, \"output2\");\n		    while (<F>)\n		    \
  {\n			$l=$_;\n			if ( $l=~/Query\\s+(\\S+)/)\n		\
	  {\n			    $seq1.=$1;\n			  }\n			elsif ( $l=~/S\
bjct\\s+(\\S+)/)\n			  {\n			    $seq2.=$1;\n			  \
}\n		      }\n		    close (F);\n		    unlink (\"ou\
tput2\");\n		    if ($seq1 ne \"\" && $seq2 ne \"\\
")\n		      {\n			$output3=\">$A\\n$seq1\\n>$B\\n$\
seq2\\n\";\n			$output3=~s/\\./-/g;\n			open (F, \\
">result.aln\");\n			print F \"$output3\";\n			clo\
se (F);\n		      }\n		  }\n	      }\n	  }\n      }\
\n    &cache_file(\"SET\",\"\",\"$name\",\"$method\
\",\"dali\",\"result.aln\",\"EBI\");\n    &set_tem\
porary_dir (\"unset\",$mode, $method, \"result.aln\
\",$outfile);\n    myexit ($EXIT_SUCCESS);\n  }\ns\
ub seq2pdb_pair\n  {\n    my ($mode, $pdbfile1, $p\
dbfile2, $method, $param, $outfile)=@_;\n\n    &se\
t_temporary_dir (\"set\",$pdbfile1,\"pdb1.pdb\",$p\
dbfile2,\"pdb2.pdb\");\n    if ($method eq \"t_cof\
fee\")\n      {\n	&safe_system (\"t_coffee -in Ppd\
b1.pdb Ppdb2.pdb -quiet -outfile=result.aln\");\n \
     }\n    elsif ( $method eq \"DaliLite\")\n    \
  {\n	if ( &safe_system (\"DaliLite -pairwise pdb1\
.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCESS)\n	  {\n	   \
  my ($seq1, $seq2);\n	     $seq1=$seq2=\"\";\n\n	\
     open (F, \"tmp1\");\n	     while (<F>)\n	    \
   {\n		 $l=$_;\n		 if ( $l=~/Query\\s+(\\S+)/)\n	\
	   {\n		     $seq1.=$1;\n		   }\n		 elsif ( $l=~/\
Sbjct\\s+(\\S+)/)\n		   {\n		     $seq2.=$1;\n		  \
 }\n	       }\n	     close (F);\n	     unlink (\"t\
mp1\");\n	     if ($seq1 ne \"\" && $seq2 ne \"\")\
\n	       {\n		 my $output3=\">$A\\n$seq1\\n>$B\\n\
$seq2\\n\";\n		 $output3=~s/\\./-/g;\n		 open (F, \
\">result.aln\");\n		 print F \"$output3\";\n		 cl\
ose (F);\n	       }\n	   }\n	else\n	  {\n	    prin\
t \"ERROR: DalLite failed to align the considered \
structures[tc_generic_method.pl]\\n\";\n	  }\n    \
  }\n    elsif ( $method eq \"TMalign\")\n      {\\
n	if ( &safe_system (\"TMalign pdb1.pdb pdb2.pdb >\
tmp1\")==$EXIT_SUCCESS)\n	  {\n	    `tail -4 tmp1 \
> tmp2`;\n\n	    open (F, \"tmp2\");\n	    while (\
<F>)\n	      {\n		unshift(@l, $_);\n	      }\n	   \
 close (F);\n	    open (F, \">result.aln\");\n	   \
 $l[3]=~s/[^a-zA-Z0-9-]/\\-/g;\n	    $l[1]=~s/[^a-\
zA-Z0-9-]/\\-/g;\n	    print F \">$A\\n$l[3]\\n>$B\
\\n$l[1]\\n\";\n	    close (F);\n	  }\n	else\n	  {\
\n	    print \"ERROR: TMalign failed to align the \
considered structures[tc_generic_method.pl]\\n\";\\
n	    `rm result.aln >/dev/null 2>/dev/null`;\n	  \
}\n      }\n    elsif ( $method eq \"mustang\")\n \
     {\n	if ( &safe_system (\"mustang -i pdb1.pdb \
pdb2.pdb -F fasta >/dev/null 2>/dev/null\")==$EXIT\
_SUCCESS)\n	  {\n	    `mv results.afasta result.al\
n`;\n	  }\n	else\n	  {\n	    print \"ERROR: mustan\
g failed to align the considered structures[tc_gen\
eric_method.pl]\\n\";\n	    `rm result.aln >/dev/n\
ull 2>/dev/null`;\n	  }\n      }\n    else\n      \
{\n	if ( &safe_system (\"$method -pdbfile1=pdb1.pe\
p -pdbfile2=pdb2.pdb -outfile=result.aln $param>x \
2>x\")==$EXIT_SUCCESS)\n	  {\n	    `mv results.afa\
sta result.aln`;\n	  }\n	else\n	  {\n	    print \"\
ERROR: $method failed to align the considered stru\
ctures[tc_generic_method.pl]\\n\";\n	    `rm resul\
t.aln >/dev/null 2>/dev/null`;\n	  }\n      }\n   \
 &set_temporary_dir (\"unset\",$mode, $method, \"r\
esult.aln\",$outfile);\n    myexit ($EXIT_SUCCESS)\
;\n  }\n\nsub seq2rnapdb_pair\n  {\n    my ($mode,\
 $pdbfile1, $pdbfile2, $method, $param, $outfile)=\
@_;\n    \n    if ($method eq \"runsara.py\")\n   \
   {\n	my $path=$ENV{PATH};\n	\n	if ($ENV{X3DNA_4_\
SARA}){$ENV{PATH}=\"$ENV{X3DNA_4_SARA}:$path\";}\n\
	\n	open(TMP,\"<$pdbfile1\");\n	my $count = 0;\n	m\
y $line;\n	while (<TMP>)\n	  {\n	    $line = $_;\n\
	    if ($count ==1)\n	      {\n		last;\n	      }\\
n	    $count += 1;\n	  }\n	\n	\n	$chain1 = substr(\
$line,length($line)-3,1);\n	\n	close TMP;\n	open(T\
MP,\"<$pdbfile2\");\n	my $count = 0;\n	while (<TMP\
>)\n	  {\n	    $line = $_;\n	    if ($count ==1)\n\
	      {\n		last;\n	      }\n	    $count += 1;\n	 \
 }\n	$chain2 = substr($line,length($line)-3,1);\n	\
close TMP;\n	\n	$tmp_file=&vtmpnam();\n	\n	safe_sy\
stem(\"runsara.py $pdbfile1 $chain1 $pdbfile2 $cha\
in2 -s -o $tmp_file --limitation 5000 > /dev/null \
2> /dev/null\");\n	if ($ENV{X3DNA_4_SARA}){$ENV{PA\
TH}=$path;}\n	\n	open(TMP,\"<$tmp_file\") or die \\
"cannot open the sara tmp file:$!\\n\";\n	open(OUT\
,\">$outfile\") or die \"cannot open the $outfile \
file:$!\\n\";\n	\n	my $switch = 0;\n	my $seqNum = \
0;\n	foreach my $line (<TMP>)\n	  {\n	    next unl\
ess ($line=~/SARAALI/);\n	    if ($line=~/>/)\n	  \
    {\n		$switch =0;\n		print OUT \">seq$seqNum\\n\
\";\n		$seqNum++;\n	      }\n	    if ($switch < 2)\
{\n	      $switch++;\n	      next;\n	    }\n	    \\
n	    if ($line =~/REMARK\\s+SARAALI\\s+([^\\*]+)\\
\*/)\n	      {\n		my $string = $1;\n		print OUT \"\
$string\\n\";\n	      }\n	  }\n	close TMP;\n	close\
 OUT;\n	unlink($tmp_file);\n      }\n  }\nsub seq2\
profile_pair\n  {\n    my ($mode, $profile1, $prof\
ile2, $method, $param, $outfile)=@_;\n    \n    \n\
    if ($method eq \"clustalw\")\n      {\n	`clust\
alw -profile1=$profile1 -profile2=$profile2 -outfi\
le=$outfile`;\n      }\n    elsif ( $method eq \"c\
lustalo\")\n      {\n	\n	`clustalo --p1 $profile1 \
--p2 $profile2 -o $outfile --force`;\n      }\n   \
 elsif ( $method eq \"hhalign\")\n      {\n	hhalig\
n ( $profile1,$profile2,$outfile,$param);\n      }\
\n    else\n      {\n	`$method -profile1=$profile1\
 -profile2=$profile2 -outfile=$outfile $param> /de\
v/null 2>/dev/null`;\n      }\n    myexit ($EXIT_S\
UCCESS);\n  }\n\nsub pg_is_installed\n  {\n    my \
@ml=@_;\n    my ($r, $p, $m);\n    my $supported=0\
;\n\n    my $p=shift (@ml);\n    if ($p=~/::/)\n  \
    {\n	if (safe_system (\"perl -M$p -e 1\")==$EXI\
T_SUCCESS){return 1;}\n	else {return 0;}\n      }\\
n    else\n      {\n	my $cwhich=vtmpnam();\n	$r=`w\
hich $p >$cwhich 2>/dev/null`;\n	$r=file2string ($\
cwhich);\n	if ($r=~/^\\//){$r=1;}\n	else {$r=0;}\n\
	\n	if ($r==0 && is_blast_package ($p)){return pg_\
is_installed (\"legacy_blast.pl\");}\n	else {retur\
n $r;}\n      }\n  }\n\nsub is_blast_package\n  {\\
n    my $p=shift;\n    if ( $p=~/blastp/){return 1\
;}\n    elsif ($p=~/blastall/){return 1;}\n    els\
if ($p=~/blastn/){return 1;}\n    elsif ($p=~/blas\
tx/){return 1;}\n    elsif ($p=~/formatdb/){return\
 1;}\n    else {return 0;}\n  }\n\nsub check_inter\
net_connection\n  {\n    my $internet;\n    my $tm\
p;\n    &check_configuration ( \"wget\");\n\n    $\
tmp=&vtmpnam ();\n\n    if     (&pg_is_installed  \
  (\"wget\")){`wget www.google.com -O$tmp >/dev/nu\
ll 2>/dev/null`;}\n    elsif  (&pg_is_installed   \
 (\"curl\")){`curl www.google.com -o$tmp >/dev/nul\
l 2>/dev/null`;}\n\n    if ( !-e $tmp || -s $tmp <\
 10){$internet=0;}\n    else {$internet=1;}\n    i\
f (-e $tmp){unlink $tmp;}\n\n    return $internet;\
\n  }\nsub check_pg_is_installed\n  {\n    my @ml=\
@_;\n    my $r=&pg_is_installed (@ml);\n    if (!$\
r && $p=~/::/)\n      {\n	print STDERR \"\\nYou Mu\
st Install the perl package $p on your system.\\nR\
UN:\\n\\tsudo perl -MCPAN -e 'install $pg'\\n\";\n\
      }\n    elsif (!$r)\n      {\n	myexit(flush_e\
rror(\"\\nProgram $p Supported but Not Installed o\
n your system\"));\n      }\n    else\n      {\n	r\
eturn 1;\n      }\n  }\nsub set_temporary_dir\n  {\
\n    my @list=@_;\n    my $dir_mode, $a, $mode, $\
method;\n\n    $dir_mode=shift (@list);\n\n\n    i\
f ( $dir_mode eq \"set\")\n      {\n	$initial_dir=\
cwd();\n	if ( !$tmp_dir)\n	  {\n	    $rand=rand (1\
00000);\n	    $tmp_dir=\"$TMPDIR/tmp4tcoffee_profi\
le_pair_dir_$$\\_P_$rand\";\n	  }\n	if ( !-d $tmp_\
dir)\n	  {\n	    push (@TMPDIR_LIST, $tmp_dir);\n	\
    `mkdir $tmp_dir`;\n	  }\n\n	for ( $a=0; $a<=$#\
list; $a+=2)\n	      {\n		if (-e $list[$a]){ `cp $\
list[$a] $tmp_dir/$list[$a+1]`;}\n	      }\n	chdir\
 $tmp_dir;\n      }\n    elsif ( $dir_mode eq \"un\
set\")\n      {\n	$mode=shift (@list);\n	$method=s\
hift (@list);\n\n	if (!-e $list[0])\n	  {\n	   mye\
xit(flush_error(\"Program $method failed to produc\
e $list[1]\" ));\n	    myexit ($EXIT_FAILURE);\n	 \
 }\n	else\n	  {\n	    chdir $initial_dir;\n	    # \
`t_coffee -other_pg seq_reformat -in $tmp_dir/$lis\
t[0] -output fasta_aln -out $tmp_dir/result2.aln`;\
\n	    `cp $tmp_dir/$list[0] $tmp_dir/result2.aln`\
;\n	    if ( $list[1] eq \"stdout\")\n	      {\n		\
open (F, \"$tmp_dir/result2.aln\");\n		while (<F>)\
{print $_;}close(F);\n	      }\n	    else\n	      \
{\n		`mv $tmp_dir/result2.aln $list[1]`;\n	      }\
\n	    shift (@list); shift (@list);\n	    foreach\
 $f (@list)\n	      {\n		if (-e (\"$tmp_dir/$f\"))\
{`mv $tmp_dir/$f .`;}\n	      }\n	  }\n      }\n  \
}\n\n\n\n\nsub my_get_opt\n  {\n    my @list=@_;\n\
    my ($cl, $a, $argv, @argl);\n\n    \n    @argl\
=();\n    $cl=shift @list;\n    for ( my $a=0; $a<\
=$#list; $a+=3)\n      {\n	my $option=$list[$a];\n\
	my $optional=$list[$a+1];\n	my $status=$list[$a+2\
];\n	my $argv=\"\";\n	if ($cl=~/$option(\\S+)/){$a\
rgv=$1;}\n	@argl=(@argl,$argv);\n\n\n	#$optional:0\
=>optional\n	#$optional:1=>must be set\n	#$status:\
 0=>no requirement\n	#$status: 1=>must be an exist\
ing file\n	#$status: 2=>must be an installed packa\
ge\n	\n\n	if ($optional==0){;}\n	elsif ( $optional\
==1 && $argv eq \"\")\n	  {\n	    myexit(flush_err\
or( \"ERROR: Option $option must be set\"));\n	   \
 myexit ($EXIT_FAILURE);\n	  }\n	if ($status==0){;\
}\n	elsif ($status ==1 && $argv ne \"\" && !-e $ar\
gv)\n	  {\n	    myexit(flush_error( \"File [$argv]\
 must exist\"));\n	    myexit ($EXIT_FAILURE);\n	 \
 }\n	elsif ( $status==2 && $argv ne \"\" && &check\
_pg_is_installed ($argv)==0)\n	  {\n	    myexit(fl\
ush_error( \" $argv is not installed\"));\n	    my\
exit ($EXIT_FAILURE);\n	  }\n      }\n    return @\
argl;\n    }\n\nsub check_file\n  {\n    my ($file\
, $msg)=@_;\n\n    if ( !-e $file)\n      {\n	myex\
it(flush_error(\"$msg\"));\n      }\n    }\nsub hh\
align\n  {\n    my ($aln1, $aln2, $outfile, $param\
)=@_;\n    my $hh=$ENV{\"HHALIGN_4_TCOFFEE\"};\n  \
  \n    \n    if ($hh)\n      {\n	\n	#external_hha\
lign\n	# set via HHALIGN_4_TCOFFEE\n	#<pg> -profil\
e1 <fasta_prf with seq1 top> -profile2 <fasta prof\
ile with seq2 top> -outfile < fasta alignmentof se\
q1 and 2 | tc_lib of seq 1 and 2>\n	\n	safe_system\
 (\"$hh -method=align -profile1=$aln1 -profile2=$a\
ln2 -outfile=$outfile\");\n      }\n    else\n    \
  {\n	&local_hhalign ($aln1, $aln2, $outfile, $par\
am);\n      }\n  }\n\n    \n    \nsub local_hhalig\
n\n  {\n    my ($aln1, $aln2, $outfile, $param)=@_\
;\n    my $h1, $h2;\n\n    $h{0}{index}=0;\n    $h\
{1}{index}=1;\n\n    $h{0}{aln}=$aln1;\n    $h{1}{\
aln}=$aln2;\n\n\n\n    %{$h{0}}=aln2psi_profile (%\
{$h{0}});\n    %{$h{1}}=aln2psi_profile (%{$h{1}})\
;\n\n    $param=~s/#S/ /g;\n    $param=~s/#M/\\-/g\
;\n    $param=~s/#E/\\=/g;\n\n\n\n    $command=\"h\
halign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tm\
p -rank 1 -mapt 0 $param\";\n    `$command`;\n\n  \
#  `hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfi\
le.tmp -rank 1 -mapt 0 -gapf 0.8 -gapg 0.8`;\n\n\n\
    # To run global use the following\n\n    open \
(I, \"$outfile.tmp\");\n    open (O, \">$outfile\"\
);\n    $h{0}{cons}=s/\\./x/g;\n    $h{1}{cons}=s/\
\\./x/g;\n\n    print O \"! TC_LIB_FORMAT_01\\n2\\\
n$h{0}{name} $h{0}{len} $h{0}{seq}\\n$h{1}{name} $\
h{1}{len} $h{1}{seq}\\n#1 2\\n\";\n\n    while (<I\
>)\n      {\n	if (/(\\d+)\\s+(\\d+)\\s+(\\d+)/)\n	\
  {\n	    print O \"\\t$h{0}{$1}\\t$h{1}{$2}\\t$3\\
\n\";\n	  }\n      }\n    print O \"! SEQ_1_TO_N\\\
n\";\n\n    close (O);\n    close (I);\n  }\n\nsub\
 aln2psi_profile\n  {\n    my (%h)=@_;\n    my ($a\
ln,$i,$hv, $a, @c, $n);\n\n\n    $i=$h{index};\n  \
  $aln=$h{aln};\n\n    `cp $aln $$.hhh_aln`;\n    \
$command=\"t_coffee -other_pg seq_reformat -in $al\
n -output hasch\";\n    $hv=`$command`;chomp ($hv)\
;\n\n    $h{a2m}=\"$tmp/$hv.tmp4hhpred.a2m\";\n   \
 $h{a3m}=\"$tmp/$hv.tmp4hhpred.a3m\";\n    if ( -e\
 $h{a3m}){;}\n    else\n      {\n	$x=`which hhcons\
ensus`;\n	`hhconsensus  -M 50 -i $h{aln} -oa2m $h{\
a2m}`;\n	if (!-e $h{a2m})\n	  {\n	    print STDERR\
 \"Program tc_generic_method.pl FAILED to run:\\n\\
\thhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}\";\n\
	    myexit ($EXIT_FAILURE);\n	  }\n\n	`hhconsensu\
s  -M 50 -i $h{aln} -oa3m $h{a3m}`;\n	if (!-e $h{a\
3m})\n	  {\n	    print STDERR \"Program tc_generic\
_method.pl FAILED to run:\\n\\thhconsensus  -M 50 \
-i $h{aln} -oa3m $h{a3m}\";\n	    myexit ($EXIT_FA\
ILURE);\n	  }\n       `buildali.pl $h{a3m} -n 1`;\\
n      }\n\n\n    $h{a2m_seq}=`head -n 2 $h{a2m} |\
 grep -v \">\"`;chomp ($h{a2m_seq});\n    $h{a3m_s\
eq}=`head -n 2 $h{a3m} | grep -v \">\"`;chomp ($h{\
a3m_seq});\n    $h{cons}=$h{a2m_seq};\n    $h{seq}\
=`head -n 2 $h{aln} | grep -v \">\"`;chomp ($h{seq\
});\n\n\n\n    @c=split (//, $h{cons});\n    $h{le\
n}=$#c+1;\n    for ($n=0,$a=0, $b=0; $a<$h{len};$a\
++)\n      {\n	if ( $c[$a]=~/[A-Z]/)\n	  {\n	    $\
h{++$n}=++$b;\n\n	  }\n	elsif ( $c[$a]=~/[a-z\\.]/\
)\n	  {\n	    ++$b;\n	  }\n      }\n\n    $name=`h\
ead -n 2 $h{aln} | grep \">\"`;\n    $name=~/\\>(\\
\S+)/;\n    $h{name}=$1;\n\n    `cp $h{a2m} $i.a2m\
`;\n    `cp $h{a3m} $i.a3m`;\n    `cp $h{aln} $i.h\
h_aln`;\n\n    return %h;\n  }\nsub read_fasta_seq\
_index\n  {\n    my $f=@_[0];\n    my %hseq;\n    \
my (@seq, @com, @name);\n    my ($a, $s,$nseq);\n\\
n    open (F, $f);\n    while (<F>)\n      {\n	$s.\
=$_;\n      }\n    close (F);\n\n\n    @name=($s=~\
/>(\\S*).*\\n[^>]*/g);\n\n    @seq =($s=~/>.*.*\\n\
([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)/g\
);\n\n\n    $nseq=$#name+1;\n\n    for ($a=0; $a<$\
nseq; $a++)\n      {\n	my $s;\n	my $n=$name[$a];\n\
	$hseq{$a}{name}=$n;\n	$seq[$a]=~s/[^A-Za-z]//g;\n\
	$hseq{$a}{order}=$a;\n	$hseq{$a}{seq}=$seq[$a];\n\
	$hseq{$a}{com}=$com[$a];\n\n      }\n    return %\
hseq;\n  }\nsub read_fasta_seq\n  {\n    my $f=@_[\
0];\n    my %hseq;\n    my (@seq, @com, @name);\n \
   my ($a, $s,$nseq);\n\n    open (F, $f);\n    wh\
ile (<F>)\n      {\n	$s.=$_;\n      }\n    close (\
F);\n\n\n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n\n\
    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s\
=~/>\\S*(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name+1\
;\n\n    for ($a=0; $a<$nseq; $a++)\n      {\n	my \
$s;\n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$se\
q[$a]=~s/[^A-Za-z]//g;\n	$hseq{$n}{order}=$a;\n	$h\
seq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\\
n\n      }\n    return %hseq;\n  }\n\n\nsub read_f\
asta_aln\n  {\n    my $f=@_[0];\n    my %hseq;\n  \
  my (@seq, @com, @name);\n    my ($a, $s,$nseq);\\
n\n    open (F, $f);\n    while (<F>)\n      {\n	$\
s.=$_;\n      }\n    close (F);\n\n\n    @name=($s\
=~/>(\\S*).*\\n[^>]*/g);\n\n    @seq =($s=~/>.*.*\\
\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)\
/g);\n\n\n    $nseq=$#name+1;\n\n    for ($a=0; $a\
<$nseq; $a++)\n      {\n	my $s;\n	my $n=$name[$a];\
\n	$hseq{$n}{name}=$n;\n	$seq[$a]=~s/[^A-Za-z-.()[\
\\]]//g;\n	$hseq{$n}{order}=$a;\n	$hseq{$n}{seq}=$\
seq[$a];\n	$hseq{$n}{com}=$com[$a];\n\n      }\n  \
  return %hseq;\n  }\n\nsub recode_name2\n{\n	my (\
$in)=shift;\n	my $mode=shift;\n\n	my %seq;\n	my $n\
ew_name;\n\n	if (! -e $in){return;}\n\n	#needed by\
 ClustalOmega to avoid very long names\n	open (INF\
ILE, \"+<$in\");\n\n	my $line;\n\n	if ($mode eq \"\
code\")\n	{\n		chomp($line = <INFILE>);\n		my $lin\
e_length = length($line);\n		$new_name=++$RECODE_N\
;\n		$new_name=\">$new_name\";\n		my $new_length =\
 length($new_name);\n		$RECODE {$new_name}=$line;\\
n		for ($count = $new_length; $count < $line_lengt\
h; $count++)\n		{\n			$new_name .= \" \";\n		}\n		\
$new_name=\"$new_name\\n\";\n		seek INFILE, 0, 0\n\
			or die \"could not seek: $!\";\n		print INFILE \
\"$new_name\";\n	}\n	else\n	{\n		my $n_found = 0;\\
n		my $file_pos=0;\n		$file_pos=tell INFILE;\n		wh\
ile (<INFILE>)\n		{\n			$line=$_;\n			$line =~ s/\\
\s*$//;\n\n			$old_name= $RECODE{$line};\n			if ($\
old_name ne \"\")\n			{\n				seek INFILE, $file_po\
s, 0\n					or die \"could not seek: $!\";\n				pri\
nt INFILE \"$old_name\\n\";\n				$file_pos++;\n			\
	if ($file_pos == 2)\n				{\n					print \"stop\\n\\
";\n					break;\n				}\n			}\n			$file_pos=tell IN\
FILE;\n		}\n\n	}\n\n\n	close INFILE;\n}\n\n\nsub r\
ecode_name\n{\n	my ($in)=shift;\n	my $mode=shift;\\
n	my $f=new FileHandle;\n	my %seq;\n	my $new_name;\
\n\n	if (! -e $in){return;}\n\n	#needed by Clustal\
Omega to avoid very long names\n	%seq=read_fasta_a\
ln ($in);\n\n	open ($f, \">$in\");\n	foreach my $s\
 (keys(%seq))\n	{\n		if ($mode eq \"code\")\n		{\n\
			$new_name=++$RECODE_N;\n			$RECODE {$new_name}=\
$seq{$s}{name};\n		}\n		else\n		{\n			$new_name=$R\
ECODE{$seq{$s}{name}};\n		}\n		print $f \">$new_na\
me\\n$seq{$s}{seq}\\n\";\n	}\n	close $f;\n}\n\nsub\
 fasta_hash2index_hash\n  {\n    my %s1=@_;\n    m\
y %s;\n    foreach my $s (keys (%s1))\n      {\n	m\
y $i=$s1{$s}{order};\n	$s{$i}{name}=$s;\n	$s{$i}{s\
eq}=$s1{$s}{seq};\n	$s{$i}{len}=length( $s{$i}{seq\
});\n	$s{n}++;\n      }\n    return %s;\n  }\nsub \
file_contains\n  {\n    my ($file, $tag, $max)=(@_\
);\n    my ($n);\n    $n=0;\n\n    if ( !-e $file \
&& ($file =~/$tag/)) {return 1;}\n    elsif ( !-e \
$file){return 0;}\n    else\n      {\n	open (FC, \\
"$file\");\n	while ( <FC>)\n	  {\n	    if ( ($_=~/\
$tag/))\n	      {\n		close (FC);\n		return 1;\n	  \
    }\n	    elsif ($max && $n>$max)\n	      {\n		c\
lose (FC);\n		return 0;\n	      }\n	    $n++;\n	  \
}\n      }\n    close (FC);\n    return 0;\n  }\n\\
n\nsub file2string\n  {\n    my $f=@_[0];\n    my \
$string, $l;\n    open (F,\"$f\");\n    while (<F>\
)\n      {\n\n	$l=$_;\n	#chomp ($l);\n	$string.=$l\
;\n      }\n    close (F);\n    $string=~s/\\r\\n/\
/g;\n    $string=~s/\\n//g;\n    return $string;\n\
  }\n\n\nsub tag2value\n  {\n\n    my $tag=(@_[0])\
;\n    my $word=(@_[1]);\n    my $return;\n\n    $\
tag=~/$word=\"([^\"]+)\"/;\n    $return=$1;\n    r\
eturn $return;\n  }\n\nsub hit_tag2pdbid\n  {\n   \
 my $tag=(@_[0]);\n    my $pdbid;\n\n    $tag=~/id\
=\"(\\S+)\"/;\n    $pdbid=$1;\n    $pdbid=~s/_//;\\
n    return $pdbid;\n  }\nsub id2pdbid\n  {\n    m\
y $in=@_[0];\n    my $id;\n\n    $in=~/(\\S+)/;\n \
   $id=$in;\n    $id=~s/PDB/pdb/g;\n\n    if ($id \
=~/pdb(.*)/){$id=$1;}\n    elsif ( $id=~/(\\S+)\\s\
+mol:protein/){$id=$1;}\n    $id=~s/[:|��_]//g\
;\n    return $id;\n  }\nsub set_blast_type\n  {\n\
    my $file =@_[0];\n    if (&file_contains ($fil\
e,\"EBIApplicationResult\",100)){$BLAST_TYPE=\"EBI\
\";}\n    elsif (&file_contains ($file,\"NCBI_Blas\
tOutput\",100)) {$BLAST_TYPE=\"NCBI\";}\n    else\\
n      {\n	$BLAST_TYPE=\"\";\n      }\n    return \
$BLAST_TYPE;\n  }\nsub is_valid_blast_xml\n    {\n\
      my $file=shift;\n      my $line;\n\n\n      \
if ( !-e $file) {return 0;}\n      $line=&file2tai\
l ($file,100);\n\n      if ( $line=~/<\\/EBIApplic\
ationResult/ || $line=~/<\\/NCBI_BlastOutput/ || $\
line=~/<\\/BlastOutput/ ){return 1;}\n      return\
 0;\n    }\nsub file2blast_flavor\n      {\n	my $f\
ile=shift; \n	if (&file_contains ($file,\"EBIAppli\
cationResult\",100)){return \"EBI\";}\n	elsif (&fi\
le_contains ($file,\"NCBI_BlastOutput\",100)){retu\
rn \"NCBI\";}\n	else {return \"UNKNOWN\";}\n      \
}\nsub blast2prf\n	{\n	  my ($blastF, $seqF,$outfi\
le)=@_;\n	  my (%s, %profile);\n	  my ($result,$ps\
iblast_output,$profile_name,@profiles);\n	  %s=rea\
d_fasta_seq_index ($seqF);\n	  my ($z1,$z1m)=uncom\
press($blastF);\n	  %profile=blast_xml2profile($s{\
0}{name}, $s{0}{seq},$maxid, $minid,$mincov,$blast\
F);\n	  output_profile ($outfile, \\%profile, $tri\
m);\n	  compress($z1,$z1m);\n	}\nsub blast_xml2pro\
file\n  {\n    my ($name,$seq,$maxid, $minid, $min\
cov, $file)=(@_);\n    my (%p, $a, $string, $n);\n\
\n\n\n    if ($BLAST_TYPE eq \"EBI\" || &file_cont\
ains ($file,\"EBIApplicationResult\",100)){%p=ebi_\
blast_xml2profile(@_);}\n    elsif ($BLAST_TYPE eq\
 \"NCBI\" || &file_contains ($file,\"NCBI_BlastOut\
put\",100)){%p=ncbi_blast_xml2profile(@_);}\n    e\
lse\n      {\n	myexit(add_error ( $$,$$,getppid(),\
 \"BLAST_FAILURE::unkown XML\",$CL));\n      }\n  \
  for ($a=0; $a<$p{n}; $a++)\n      {\n	my $name=$\
p{$a}{name};\n	$p{$name}{seq}=$p{$a}{seq};\n	$p{$n\
ame}{index}=$a;\n      }\n    return %p;\n  }\nsub\
 ncbi_tblastx_xml2lib_file\n  {\n    my  ($outlib,\
$string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$i,$nhi\
ts,@identifyerL);\n    my (%ITERATION);\n\n    ope\
n (F, \">>$outlib\");\n\n    $seq=~s/[^a-zA-Z]//g;\
\n    $L=length ($seq);\n\n    %ITERATION=xml2tag_\
list ($string, \"Iteration\");\n    for ($i=0; $i<\
$ITERATION{n};$i++)\n      {\n	my ($qindex, $qlen,\
 %hit, $string);\n	$string=$ITERATION{$i}{body};\n\
\n	$qindex=xmltag2value($string,\"Iteration_iter-n\
um\");\n	$qlen  =xmltag2value($string,\"Iteration_\
query-len\");\n	%hit=&xml2tag_list  ($string, \"Hi\
t\");\n\n	for ($a=0; $a<$hit{n}; $a++)\n	  {\n	   \
 my ($string);\n	    $string=$hit{$a}{body};\n\n	 \
   $hindex=xmltag2value($string,\"Hit_accession\")\
+1;\n	    if ($hindex<=$qindex){next;}\n	    else \
 {print F  \"# $qindex $hindex\\n\";}\n\n\n	    $h\
len=xmltag2value  ($string,\"Hit_len\");\n	    %HS\
P=&xml2tag_list  ($string, \"Hsp\");\n\n	    for (\
$b=0; $b<$HSP{n}; $b++)\n	      {\n		my ($string, \
$qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e);\n		$string=$H\
SP{$b}{body};\n\n		$qs=xmltag2value  ($string,\"Hs\
p_query-from\");\n		$qe=xmltag2value  ($string,\"H\
sp_query-to\");\n		$qf=xmltag2value  ($string,\"Hs\
p_query-frame\");\n\n		$hs=xmltag2value  ($string,\
\"Hsp_hit-from\");\n		$he=xmltag2value  ($string,\\
"Hsp_hit-to\");\n		$hf=xmltag2value  ($string,\"Hs\
p_hit-frame\");\n\n		$s=xmltag2value  ($string,\"H\
sp_identity\");\n		$l=xmltag2value  ($string,\"Hsp\
_align-len\");\n		$s=int(($s*100)/$l);\n\n		if ($q\
f>0)\n		  {$rqs=$qs; $rqe=$qe;}\n		else\n		  {\n		\
    $rqe=($qlen-$qs)+1;\n		    $rqs=($qlen-$qe)+1;\
\n		  }\n\n		if ($hf>0)\n		  {$rhs=$hs; $rhe=$he;}\
\n		else\n		  {\n		    $rhe=($hlen-$hs)+1;\n		    \
$rhs=($hlen-$he)+1;\n		  }\n		for ($d=0,$e=$rqs; $\
e<$rqe; $e++,$d++)\n		  {\n		    my ($r1,$r2);\n		\
    $r1=$e;\n		    $r2=$rhs+$d;\n		    print F \" \
$r1 $r2 $s 0\\n\";\n		  }\n	      }\n	  }\n      }\
\n    print F \"! SEQ_1_TO_N\\n\";\n\n    close (F\
);\n    return %lib;\n  }\n\nsub ncbi_tblastpx_xml\
2lib_file\n  {\n    my  ($outlib,$string,%s)=(@_);\
\n    my ($L,$l, $a,$b,$c,$d,$i,$nhits,@identifyer\
L);\n    my (%ITERATION,%hdes, %qdes);\n\n    open\
 (F, \">>$outlib\");\n\n    $seq=~s/[^a-zA-Z]//g;\\
n    $L=length ($seq);\n\n    %ITERATION=xml2tag_l\
ist ($string, \"Iteration\");\n    for ($i=0; $i<$\
ITERATION{n};$i++)\n      {\n	my ($qindex, $qlen, \
%hit, $string);\n	$string=$ITERATION{$i}{body};\n\\
n	$qdef=xmltag2value($string,\"Iteration_query-def\
\");\n	%qdes=&tblastpx_name2description($qdef,%s);\
\n	$qlen  =xmltag2value($string,\"Iteration_query-\
len\");\n	%hit=&xml2tag_list  ($string, \"Hit\");\\
n\n	for ($a=0; $a<$hit{n}; $a++)\n	  {\n	    my ($\
string);\n	    $string=$hit{$a}{body};\n	    $hdef\
=xmltag2value($string,\"Hit_def\");\n	    %hdes=&t\
blastpx_name2description($hdef,%s);\n	    if ($hde\
s{index}<=$qdes{index}){next;}\n	    else  {print \
F  \"# $qdes{index} $hdes{index}\\n\";}\n\n\n	    \
$hlen=xmltag2value  ($string,\"Hit_len\");\n	    %\
HSP=&xml2tag_list  ($string, \"Hsp\");\n\n	    for\
 ($b=0; $b<$HSP{n}; $b++)\n	      {\n		my ($string\
, $l,$qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e, @s1, @s2)\
;\n		$string=$HSP{$b}{body};\n\n		$qs=xmltag2value\
  ($string,\"Hsp_query-from\");\n		$qe=xmltag2valu\
e  ($string,\"Hsp_query-to\");\n		$qf=$qdes{frame}\
;\n		$qseq=xmltag2value  ($string,\"Hsp_qseq\");\n\
\n		$hs=xmltag2value  ($string,\"Hsp_hit-from\");\\
n		$he=xmltag2value  ($string,\"Hsp_hit-to\");\n		\
$hf=$hdes{frame};\n		$hseq=xmltag2value  ($string,\
\"Hsp_hseq\");\n\n		$s=xmltag2value  ($string,\"Hs\
p_identity\");\n		$l=xmltag2value  ($string,\"Hsp_\
align-len\");\n		$s=int(($s*100)/$l);\n		@s1=tblas\
tpx_hsp2coordinates($qseq,$qs,$qe,%qdes);\n		@s2=t\
blastpx_hsp2coordinates($hseq,$hs,$he,%hdes);\n\n\\
n		for ($f=0; $f<=$#s1; $f++)\n		  {\n		    if ($s\
1[$f]==-1 || $s2[$f]==-1){next;}\n		    else\n		  \
    {\n			print F \" $s1[$f] $s2[$f] $s 0\\n\";\n	\
	      }\n		  }\n	      }\n	  }\n      }\n    prin\
t F \"! SEQ_1_TO_N\\n\";\n\n    close (F);\n    re\
turn %lib;\n  }\nsub tblastpx_hsp2coordinates\n  {\
\n    my ($seq, $s, $e, %des)=@_;\n    my @list;\n\
    my @sa;\n    my @gap=(-1,-1,-1);\n\n    $s=$de\
s{start}+3*($s-1);\n\n    if ($des{strand} eq \"d\\
"){;}\n    else {$s=($des{length}-$s)+1;}\n\n    f\
oreach $c (split (//,$seq))\n      {\n	if ( $c eq \
'-'){push (@list,@gap);}\n	elsif ($des{strand} eq \
\"d\")\n	  {\n	    push(@list,$s++,$s++,$s++);\n	 \
 }\n	else\n	  {\n	    push(@list, $s--,$s--,$s--);\
\n	  }\n      }\n    return @list;\n  }\n\nsub tbl\
astpx_name2description\n  {\n    my ($name, %s)=@_\
;\n    my @at=split(\"__\", $name);\n    my %des;\\
n\n    $des{name}=$at[0];\n    $des{strand}=$at[1]\
;\n\n    $des{start}=$at[2];\n    $des{end}=$at[3]\
;\n    $des{length}=$at[4];\n    $des{index}=$s{$a\
t[0]}{order}+1;\n    return %des;\n  }\nsub ncbi_b\
last_xml2profile\n  {\n    my ($name,$seq,$maxid, \
$minid, $mincov, $string)=(@_);\n    my ($L,$l, $a\
,$b,$c,$d,$nhits,@identifyerL);\n\n    $seq=~s/[^a\
-zA-Z]//g;\n    $L=length ($seq);\n    \n    #This\
 is causing the NCBI parser to fail when Iteration\
_query-def is missing\n    #%query=&xml2tag_list (\
$string, \"Iteration_query-def\");\n    #$name=$qu\
ery{0}{body};\n\n    %hit=&xml2tag_list ($string, \
\"Hit\");\n\n    \n    for ($nhits=0,$a=0; $a<$hit\
{n}; $a++)\n      {\n	my ($ldb,$id, $identity, $ex\
pectation, $start, $end, $coverage, $r);\n	my (%ID\
,%DE,%HSP);\n\n	$ldb=\"\";\n\n	%ID=&xml2tag_list (\
$hit{$a}{body}, \"Hit_id\");\n	$identifyer=$ID{0}{\
body};\n\n	%DE=&xml2tag_list ($hit{$a}{body}, \"Hi\
t_def\");\n	$definition=$DE{0}{body};\n\n	%HSP=&xm\
l2tag_list ($hit{$a}{body}, \"Hsp\");\n	for ($b=0;\
 $b<$HSP{n}; $b++)\n	  {\n	    my (%START,%END,%E,\
%I,%Q,%M);\n\n\n	    %START=&xml2tag_list ($HSP{$b\
}{body}, \"Hsp_query-from\");\n	    %HSTART=&xml2t\
ag_list ($HSP{$b}{body}, \"Hsp_hit-from\");\n\n	  \
  %LEN=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_alig\
n-len\");\n	    %END=  &xml2tag_list ($HSP{$b}{bod\
y}, \"Hsp_query-to\");\n	    %HEND=  &xml2tag_list\
 ($HSP{$b}{body}, \"Hsp_hit-to\");\n	    %E=&xml2t\
ag_list     ($HSP{$b}{body}, \"Hsp_evalue\");\n	  \
  %I=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_iden\
tity\");\n	    %Q=&xml2tag_list     ($HSP{$b}{body\
}, \"Hsp_qseq\");\n	    %M=&xml2tag_list     ($HSP\
{$b}{body}, \"Hsp_hseq\");\n\n	    for ($e=0; $e<$\
Q{n}; $e++)\n\n	      {\n		$qs=$Q{$e}{body};\n		$m\
s=$M{$e}{body};\n\n		$expectation=$E{$e}{body};\n	\
	$identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN\
{$e}{body}*100;\n		$start=$START{$e}{body};\n		$en\
d=$END{$e}{body};\n		$Hstart=$HSTART{$e}{body};\n	\
	$Hend=$HEND{$e}{body};\n\n		$coverage=($L)?(($end\
-$start)*100)/$L:0;\n		if ($identity>$maxid || $id\
entity<$minid || $coverage<$mincov)\n		  {\n		    \
next;\n		  }\n		@lr1=(split (//,$qs));\n		@lr2=(sp\
lit (//,$ms));\n		$l=$#lr1+1;\n		for ($c=0;$c<$L;$\
c++){$p[$nhits][$c]=\"-\";}\n		for ($d=0,$c=0; $c<\
$l; $c++)\n		  {\n		    $r=$lr1[$c];\n		    if ( $\
r=~/[A-Za-z]/)\n		      {\n\n			$p[$nhits][$d + $s\
tart-1]=$lr2[$c];\n			$d++;\n		      }\n		  }\n		$\
Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$ms;\n		$Qstart\
L[$nhits]=$start;\n		$HstartL[$nhits]=$Hstart;\n		\
$identityL[$nhits]=$identity;\n		$endL[$nhits]=$en\
d;\n		$definitionL[$nhits]=$definition;\n		$identi\
fyerL[$nhits]=$identifyer;\n		$comment[$nhits]=\"$\
ldb|$identifyer [Eval=$expectation][id=$identity%]\
[start=$Hstart end=$Hend]\";\n		$nhits++;\n	      \
}\n	  }\n      }\n\n\n    $profile{n}=0;\n    $pro\
file{$profile{n}}{name}=$name;\n    $profile{$prof\
ile{n}}{seq}=$seq;\n    $profile {n}++;\n\n    for\
 ($a=0; $a<$nhits; $a++)\n      {\n	$n=$a+1;\n\n	$\
profile{$n}{name}=\"$name\\_$a\";\n	$profile{$n}{s\
eq}=\"\";\n	$profile{$n}{Qseq}=$Qseq[$a];\n	$profi\
le{$n}{Hseq}=$Hseq[$a];\n	$profile{$n}{Qstart}=$Qs\
tartL[$a];\n	$profile{$n}{Hstart}=$HstartL[$a];\n	\
$profile{$n}{identity}=$identityL[$a];\n	$profile{\
$n}{definition}=$definitionL[$a];\n	$profile{$n}{i\
dentifyer}=$identifyerL[$a];\n	$profile{$n}{commen\
t}=$comment[$a];\n\n	for ($b=0; $b<$L; $b++)\n	  {\
\n	    if ($p[$a][$b])\n	      {\n		$profile{$n}{s\
eq}.=$p[$a][$b];\n	      }\n	    else\n	      {\n	\
	$profile{$n}{seq}.=\"-\";\n	      }\n	  }\n      \
}\n\n    $profile{n}=$nhits+1;\n    return %profil\
e;\n  }\nsub ebi_blast_xml2profile\n  {\n    my ($\
name,$seq,$maxid, $minid, $mincov, $string)=(@_);\\
n    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,$i\
dentifyer);\n\n\n\n    $seq=~s/[^a-zA-Z]//g;\n    \
$L=length ($seq);\n    %hit=&xml2tag_list ($string\
, \"hit\");\n\n    for ($nhits=0,$a=0; $a<$hit{n};\
 $a++)\n      {\n	my ($ldb,$id, $identity, $expect\
ation, $start, $end, $coverage, $r);\n	my (%Q,%M,%\
E,%I);\n\n	$ldb=&tag2value ($hit{$a}{open}, \"data\
base\");\n	$identifyer=&tag2value ($hit{$a}{open},\
 \"id\");\n\n	$description=&tag2value ($hit{$a}{op\
en}, \"description\");\n\n	%Q=&xml2tag_list ($hit{\
$a}{body}, \"querySeq\");\n	%M=&xml2tag_list ($hit\
{$a}{body}, \"matchSeq\");\n	%E=&xml2tag_list ($hi\
t{$a}{body}, \"expectation\");\n	%I=&xml2tag_list \
($hit{$a}{body}, \"identity\");\n\n\n	for ($b=0; $\
b<$Q{n}; $b++)\n	  {\n\n	    $qs=$Q{$b}{body};\n	 \
   $ms=$M{$b}{body};\n\n	    $expectation=$E{$b}{b\
ody};\n	    $identity=$I{$b}{body};\n\n\n	    $sta\
rt=&tag2value ($Q{$b}{open}, \"start\");\n	    $en\
d=&tag2value ($Q{$b}{open}, \"end\");\n	    $start\
M=&tag2value ($M{$b}{open}, \"start\");\n	    $end\
M=&tag2value ($M{$b}{open}, \"end\");\n	    $cover\
age=(($end-$start)*100)/$L;\n\n	   # print \"$id: \
ID: $identity COV: $coverage [$start $end]\\n\";\n\
\n	    if ($identity>$maxid || $identity<$minid ||\
 $coverage<$mincov){next;}\n	    # print \"KEEP\\n\
\";\n\n\n	    @lr1=(split (//,$qs));\n	    @lr2=(s\
plit (//,$ms));\n	    $l=$#lr1+1;\n	    for ($c=0;\
$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n	    for ($d=0\
,$c=0; $c<$l; $c++)\n	      {\n		$r=$lr1[$c];\n		i\
f ( $r=~/[A-Za-z]/)\n		  {\n\n		    $p[$nhits][$d \
+ $start-1]=$lr2[$c];\n		    $d++;\n		  }\n	      \
}\n\n	    $Qseq[$nhits]=$qs;\n	    $Hseq[$nhits]=$\
ms;\n	    $QstartL[$nhits]=$start;\n	    $HstartL[\
$nhits]=$Hstart;\n	    $identityL[$nhits]=$identit\
y;\n	    $endL[$nhits]=$end;\n	    $definitionL[$n\
hits]=$definition;\n	    $identifyerL[$nhits]=$ide\
ntifyer;\n	    $comment[$nhits]=\"$ldb|$identifyer\
 [Eval=$expectation][id=$identity%][start=$startM \
end=$endM]\";\n	    $nhits++;\n	  }\n      }\n\n  \
  $profile{n}=0;\n    $profile{$profile{n}}{name}=\
$name;\n    $profile{$profile{n}}{seq}=$seq;\n    \
$profile {n}++;\n\n    for ($a=0; $a<$nhits; $a++)\
\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\"$name\\
\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{Q\
seq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	\
$profile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}{\
Hstart}=$HstartL[$a];\n	$profile{$n}{identity}=$id\
entityL[$a];\n	$profile{$n}{definition}=$definitio\
nL[$a];\n	$profile{$n}{identifyer}=$identifyerL[$a\
];\n	$profile{$n}{comment}=$comment[$a];\n\n	for (\
$b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	 \
     {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      }\
\n	    else\n	      {\n		$profile{$n}{seq}.=\"-\";\
\n	      }\n	  }\n      }\n    $profile{n}=$nhits+\
1;\n\n    return %profile;\n  }\nsub output_profil\
e\n  {\n    my ($outfile,$profileR, $trim)=(@_);\n\
    my ($a);\n    my %profile=%$profileR;\n    my \
$P= new FileHandle;\n    my $tmp=vtmpnam();\n\n   \
 open ($P, \">$tmp\");\n    for ($a=0; $a<$profile\
{n}; $a++)\n      {\n	print $P \">$profile{$a}{nam\
e} $profile{$a}{comment}\\n$profile{$a}{seq}\\n\";\
\n      }\n    close ($P);\n\n    if ( $trim)\n   \
   {\n	if ($ENV{psitrim_mode_4_TCOFFEE} eq \"trim\\
")# old trimming\n	  {\n	    if ($trim>0)\n	      \
{\n		&safe_system (\"t_coffee -other_pg seq_reform\
at -in $tmp -action +trim _aln_n$trim\\_K1 -output\
 fasta_aln -out $outfile\");\n	      }\n	    else\\
n	      {\n		&safe_system (\"t_coffee -other_pg se\
q_reformat -in $tmp -action +trim _aln_%%$trim\\_K\
1 -output fasta_aln -out $outfile\");\n	      }\n	\
  }\n	else # newtrimming\n	  {\n	    my $tm;\n	   \
 if ($ENV{psitrim_tree_4_TCOFFEE}){$tm=\"-treemode\
=\".$ENV{psitrim_tree_4_TCOFFEE};}\n	    else {$tm\
=\"-treemode=codnd\";}\n	    if ($trim>0)\n	      \
{\n		&safe_system (\"t_coffee -other_pg seq_reform\
at -in $tmp $tm -keep 1 -action +regtrim $trim -ou\
tput fasta_aln -out $outfile\");\n	      }\n	    e\
lse\n	      {\n		&safe_system (\"t_coffee -other_p\
g seq_reformat -in $tmp $tm -keep 1 -action +regtr\
im $trim\\% -output fasta_aln -out $outfile\");\n	\
      }\n	   \n	  }\n	\n      }\n    else\n      {\
\n	&safe_system (\"mv $tmp $outfile\");\n      }\n\
    return;\n  }\nsub blast_xml2hit_list\n  {\n   \
 my $string=(@_[0]);\n    return &xml2tag_list ($s\
tring, \"hit\");\n  }\nsub xmltag2value\n  {\n    \
my ($string_in, $tag)=@_;\n    my %TAG;\n    %TAG=\
xml2tag_list ($string_in, $tag);\n    return $TAG{\
0}{body};\n  }\n\nsub xml2tag_list\n  {\n    my ($\
string_in,$tag)=@_;\n    my $tag_in, $tag_out;\n  \
  my %tag;\n\n    if (-e $string_in)\n      {\n	$s\
tring=&file2string ($string_in);\n      }\n    els\
e\n      {\n	$string=$string_in;\n      }\n    $ta\
g_in1=\"<$tag \";\n    $tag_in2=\"<$tag>\";\n    $\
tag_out=\"/$tag>\";\n    $string=~s/>/>##1/g;\n   \
 $string=~s/</##2</g;\n    $string=~s/##1/<#/g;\n \
   $string=~s/##2/#>/g;\n    @l=($string=~/(\\<[^>\
]+\\>)/g);\n    $tag{n}=0;\n    $in=0;$n=-1;\n\n\n\
\n    foreach $t (@l)\n      {\n\n	$t=~s/<#//;\n	$\
t=~s/#>//;\n\n	if ( $t=~/$tag_in1/ || $t=~/$tag_in\
2/)\n	  {\n\n	    $in=1;\n	    $tag{$tag{n}}{open}\
=$t;\n	    $n++;\n\n	  }\n	elsif ($t=~/$tag_out/)\\
n	  {\n\n\n	    $tag{$tag{n}}{close}=$t;\n	    $ta\
g{n}++;\n	    $in=0;\n	  }\n	elsif ($in)\n	  {\n\n\
	    $tag{$tag{n}}{body}.=$t;\n	  }\n      }\n\n  \
  return %tag;\n  }\n\n\nsub seq2gor_prediction\n \
 {\n    my ($name, $seq,$infile, $outfile, $gor_se\
q, $gor_obs)=(@_);\n    my ($l);\n\n    `gorIV -pr\
d $infile -seq $gor_seq -obs $gor_obs > gor_tmp`;\\
n    open (GR, \">$outfile\");\n    open (OG, \"go\
r_tmp\");\n\n    while (<OG>)\n      {\n\n	$l=$_;\\
n	if ($l=~/\\>/){print GR \"$l\";}\n	elsif ( $l=~/\
Predicted Sec. Struct./)\n	  {\n	    $l=~s/Predict\
ed Sec. Struct\\.//;\n	    print GR \"$l\";\n	  }\\
n      }\n    close (GR);\n    close (OG);\n    re\
turn;\n  }\nsub seq2msa_tm_prediction\n  {\n    my\
 ($name, $seq, $db, $infile, $outfile, $arch, $psv\
)=(@_);\n    my (%p,%gseq,%R, $blast_output, %s, $\
l);\n    my $R2=new FileHandle;\n    my $method=\"\
psitm\";\n    my $SERVER=\"EBI\";\n\n    $blast_ou\
tput=&run_blast ($name,\"blastp\", $db, $infile, \\
"outfile\");\n\n    if (&cache_file(\"GET\",$infil\
e,$name,$method,$db,$outfile,$SERVER))\n      {\n	\
print \"\\tPSITM: USE Cache\\n\";\n	return $outfil\
e;\n      }\n    else\n      {\n	$CACHE_STATUS=\"C\
OMPUTE CACHE\";\n	%p=blast_xml2profile($name,$seq,\
$maxid, $minid,$mincov,$blast_output);\n\n\n	open \
(F, \">tm_input\");\n	for (my $a=0; $a<$p{n}; $a++\
)\n	  {\n	    my $s;\n\n	    $s=$p{$a}{seq};\n	   \
 $s=uc($s);\n	    print F \">$p{$a}{name}\\n$s\\n\\
";\n	    #print stdout \">$p{$a}{name}\\n$s\\n\";\\
n	  }\n	close (F);\n	print \"\\tPSITM: kept  $p{n}\
 Homologues for Sequence $p{0}{name}\\n\";\n	&safe\
_system (\"t_coffee -other_pg fasta_seq2hmmtop_fas\
ta.pl -in=tm_input -out=$outfile -output=cons -cov\
=70 -trim=95 -arch=$arch -psv=$psv\");\n	unlink (\\
"tm_input\");\n	&cache_file(\"SET\",$infile,$name,\
$method,$db,$outfile,$SERVER);\n	return;\n      }\\
n  }\n\n\nsub seq2msa_gor_prediction\n  {\n    my \
($name, $seq,$infile, $outfile, $gor_seq, $gor_obs\
)=(@_);\n    my (%p,%gseq,%R, $blast_output, %s, $\
l);\n    my $R2=new FileHandle;\n    my $db=\"unip\
rot\";\n    my $method=\"psigor\";\n    my $SERVER\
=\"EBI\";\n\n    $blast_output=&run_blast ($name,\\
"blastp\", \"uniprot\", $infile, \"outfile\");\n\n\
    if (&cache_file(\"GET\",$infile,$name,$method,\
$db,$outfile,$SERVER, $psiJ))\n      {\n	print \"\\
\tPSIGOR: USE Cache\\n\";\n	return $outfile;\n    \
  }\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE C\
ACHE\";\n	%p=blast_xml2profile($name,$seq,$maxid, \
$minid,$mincov,$blast_output);\n\n\n	open (F, \">g\
or_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n	  {\
\n	    my $s;\n\n	    $s=$p{$a}{seq};\n	    $s=uc(\
$s);\n	    print F \">$p{$a}{name}\\n$s\\n\";\n	  \
  #print stdout \">$p{$a}{name}\\n$s\\n\";\n	  }\n\
	close (F);\n	print \"\\tGORTM: kept  $p{n} Homolo\
gues for Sequence $p{0}{name}\\n\";\n	&safe_system\
 (\"t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -\
in=gor_input -out=$outfile -output=cons -cov=70 -t\
rim=95 -gor_seq=$gor_seq -gor_obs=$gor_obs -mode=g\
or\");\n	unlink (\"tm_input\");\n	&cache_file(\"SE\
T\",$infile,$name,$method,$db,$outfile,$SERVER);\n\
	return;\n      }\n  }\n\n\n\nsub run_blast\n  {\n\
    my ($name, $method, $db, $infile, $outfile, $r\
un)=(@_);\n    if (!$run){$run=1;}\n    my $error_\
log=vtmpnam();\n    my $cl_db;\n    my $psiJ=($ENV\
{psiJ_4_TCOFFEE})?$ENV{psiJ_4_TCOFFEE}:1;\n    \n \
   my $psiJFlag=\"-j$psiJ\";\n   \n    if (&cache_\
file(\"GET\",$infile,$name,$method,$db,$outfile,$S\
ERVER, $psiJ) && is_valid_blast_xml ($outfile))\n \
     {return $outfile;}\n    else\n      {\n	$CACH\
E_STATUS=\"COMPUTE CACHE\";\n	if ( $SERVER eq \"EB\
I_SOAP\")\n	  {\n	    &check_configuration (\"EMAI\
L\",\"SOAP::Light\",\"INTERNET\");\n\n	    $cl_met\
hod=$method;\n	    if ($cl_method =~/wu/)\n	      \
{\n		if ( $method eq \"psiblast\" || $psiJ>1){myex\
it(add_error (EXIT_FAILURE,$$,$$,getppid(), \"BLAS\
T_FAILURE::$SERVER does not Support psiblast mode \
($psiJFlag)\",$CL));}\n		$cl_method=~s/wu//;\n		if\
 ( $cl_method eq \"psiblast\" || $psiJ>1)\n		  {\n\
		    add_warning($$,$$,\"PSI BLAST cannot be used\
 with the wuBLAST Client. Use server=EBI Or server\
=LOCAL. blastp will be used instead\");\n		    $cl\
_method=\"blastp\";\n		  }\n\n		$command=\"t_coffe\
e -other_pg wublast.pl --email $EMAIL $infile -D $\
db -p $cl_method --outfile $outfile -o xml>/dev/nu\
ll 2>$error_log\";\n		&safe_system ( $command);\n	\
	if (-e \"$outfile.xml\") {`mv $outfile.xml $outfi\
le`;}\n	      }\n	    else\n	      {\n		if ($cl_me\
thod eq \"psiblast\"){$cl_method =\"blastp $psiJFl\
ag\";}\n\n		$command=\"t_coffee -other_pg blastpgp\
.pl --email $EMAIL $infile -d $db --outfile $outfi\
le -p $cl_method --mode PSI-Blast>/dev/null 2>$err\
or_log\";\n		&safe_system ( $command);\n\n		if (-e\
 \"$outfile.xml\") {`mv $outfile.xml $outfile`;}\n\
	      }\n	  }\n	elsif ($SERVER eq \"EBI_REST\" ||\
 $SERVER eq \"EBI\")\n	  {\n	    $cl_method=$metho\
d;\n	    &check_configuration(\"EMAIL\",\"XML::Sim\
ple\", \"INTERNET\");\n	    if ($db eq \"uniprot\"\
){$db1=\"uniprotkb\";}\n	    else {$db1=$db;}\n\n	\
    \n	    if ($cl_method =~/wu/)\n	      {\n		$cl\
_method=~s/wu//;\n\n		\n		if ( $cl_method eq \"psi\
blast\"){$cl_method=\"blastp\";}\n\n		$command=\"t\
_coffee -other_pg wublast_lwp.pl --email $EMAIL -D\
 $db1 -p $cl_method --outfile $outfile --align 5 -\
-stype protein $infile>/dev/null 2>error_log\";\n	\
      }\n	    else\n	      {\n		if ( $cl_method =~\
/psiblast/){$cl_method =\"blastp $psiJFlag\";}\n		\
$command=\"t_coffee -other_pg ncbiblast_lwp.pl --e\
mail $EMAIL --database $db1 --program $cl_method -\
-outfile $outfile --alignments 5 --stype protein $\
infile>/dev/null 2>$error_log\";\n		#DEBUG\n		#$co\
mmand=\"t_coffee -other_pg ncbiblast_lwp.pl --emai\
l $EMAIL -D $db1 -p $cl_method --outfile $outfile \
--align 5 --stype protein $infile\";\n		\n		my $ma\
xrun=5;#number of crashes accepetd\n		my $nrun;\n	\
	my $keep_going=1;\n		while ($keep_going)\n		  {\n\
		    \n		    #print \"----- $command [$nrun]\\n\"\
;\n		    $nrun++;\n		    $keep_going=0;\n		    &sa\
fe_system ( $command,5);\n		    \n		    my $succes\
s=0;\n		    $success =$success || -e \"$outfile.ou\
t.xml\";\n		    $success =$success || -e \"$outfil\
e.xml.xml\";\n		    $success =$success || -e \"$ou\
tfile.out..xml\";\n		    $success =$success || -e \
\"$outfile.xml..xml\";\n		    \n		    if (!$succes\
s && ($nrun<$maxrun || -e \"$outfile.out.txt\"))\n\
		      {\n			$keep_going=1;\n			add_warning($$,$$\
,\"[ncbiblast_lwp.pl] [$command] failed to produce\
 xml output -- will ne tried again [$nrun]\");\n		\
      }\n		  }\n		\n		if (-e \"$outfile.out.xml\")\
 {`mv $outfile.out.xml $outfile`;}\n		elsif (-e \"\
$outfile.xml.xml\"){`mv $outfile.xml.xml $outfile`\
;}\n		elsif (-e \"$outfile.out..xml\") {`mv $outfi\
le.out..xml $outfile`;}\n		elsif (-e \"$outfile.xm\
l..xml\"){`mv $outfile.xml..xml $outfile`;}\n		els\
e\n		  {\n		    add_warning($$,$$,\"[ncbiblast_lwp\
.pl] [$command] failed to produce xml output\");\n\
		  }\n	      }\n	  }\n	elsif ($SERVER eq \"NCBI\"\
)\n	  {\n	    &check_configuration (\"INTERNET\");\
\n	    if ($db eq \"uniprot\"){$cl_db=\"swissprot\\
";}\n	    else {$cl_db=$db;}\n\n	    if ( $method \
eq \"psiblast\" || $psiJ>1){myexit(add_error (EXIT\
_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::$SERVER\
 does not Support psiblast mode ($psiJFlag)\",$CL)\
);}\n	    my $cl_method=$method;\n	    \n	    &che\
ck_configuration ($cl_method);  \n	    $command=\"\
$cl_method -db $cl_db -query $infile -out $outfile\
 -outfmt 5 -remote\";\n	    &safe_system ($command\
);\n	  }\n	elsif ($SERVER =~/CLIENT_(.*)/)\n	  {\n\
	    my $client=$1;\n	    if ( $method eq \"psibla\
st\" || $psiJ>1){myexit(add_error (EXIT_FAILURE,$$\
,$$,getppid(), \"BLAST_FAILURE::$SERVER does not S\
upport psiblast mode ($psiJFlag)\",$CL));}\n	    $\
command=\"$client -p $method -d $db -i $infile -o \
$outfile -m 7\";\n	    &safe_system ($command);\n	\
  }\n	elsif ( $SERVER eq \"LOCAL_blastall\")\n	  {\
\n	    &check_configuration (\"blastall\");\n	    \
if ( $method eq \"psiblast\" || $psiJ>1){myexit(ad\
d_error (EXIT_FAILURE,$$,$$,getppid(), \"BLAST_FAI\
LURE::$SERVER does not Support psiblast mode ($psi\
JFlag)\",$CL));}\n	    $command=\"blastall -d $db \
-i $infile -o $outfile -m7 -p blastp\";\n	    &saf\
e_system ($command);\n	  }\n	elsif ( $SERVER eq \"\
LOCAL\")\n	  {\n	    my $legacy=0;\n	    if ($ENV{\
\"BLAST_DB_DIR\"}) \n	      {\n	    	$x=$ENV{\"BLA\
ST_DB_DIR\"};\n		$cl_db=\"$x/$db\";\n	      }\n	  \
  else\n	      {\n		$cl_db=$db;\n	      }\n	    \n\
	    ##\n	    ## BLAST+ provide different binaries\
 names and CLI options\n	    ## Use the 'legacy_bl\
ast.pl' to keep compatibility with old blast comma\
nds\n	    ##\n	    $path=`which legacy_blast.pl 2>\
/dev/null`;  \n	    $path=`dirname $path`; \n	    \
chomp($path);\n	    \n	    if    (!$legacy && ($me\
thod eq \"blastp\" || $method eq \"psiblast\"))\n	\
      {\n		\n		&check_configuration(\"psiblast\");\
\n		$command=\"psiblast -db $cl_db -query $infile \
-num_iterations $psiJ -out $outfile -outfmt 5\";\n\
	      }\n	    elsif ($legacy && $method eq \"blas\
tp\")\n	     {\n	       &check_configuration(\"leg\
acy_blast.pl\");\n	       $command=\"legacy_blast.\
pl blastpgp --path $path -d $cl_db -i $infile -o $\
outfile -m7 $psiJFlag\";		\n	     }\n	    elsif ($\
legacy && $method eq \"psiblast\")\n	      {\n		&c\
heck_configuration(\"legacy_blast.pl\");\n		$comma\
nd=\"legacy_blast.pl blastpgp --path $path -d $cl_\
db -i $infile -o $outfile -m7 $psiJFlag\";\n	     \
 }\n	    elsif ($method eq \"blastn\")\n	      {\n\
		&check_configuration(\"legacy_blast.pl\");\n		$c\
ommand=\"legacy_blast.pl blastall --path $path -p \
blastn -d $cl_db -i $infile -o $outfile -m7 -W6\";\
\n	      }\n	    &safe_system ($command);\n	  }\n	\
else\n	  {\n\n	    myexit(add_error (EXIT_FAILURE,\
$$,$$,getppid(), \"BLAST_FAILURE::UnknownServer\",\
$CL));\n	  }\n\n\n	#Check that everything went wel\
l\n\n	if ( !-e $outfile || !&is_valid_blast_xml($o\
utfile))\n	  {\n\n	    if ( -e $outfile)\n	      {\
\n		add_warning ($$,$$,\"Corrupted Blast Output (R\
un $run)\");\n		unlink($outfile);\n	      }\n	    \
if ( -e $error_log)\n	      {\n\n		my $error_msg=f\
ile2string ($error_log);\n\n		if ( $error_msg =~/e\
nter a valid email/)\n		  {\n		    myexit(add_erro\
r (EXIT_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::\
Invalid_or_rejected_email::$EMAIL\", \"$command\")\
);\n		  }\n	      }\n	    if ( $run==$BLAST_MAX_NR\
UNS)\n	      {\n\n		myexit(add_error (EXIT_FAILURE\
,$$,$$,getppid(), \"BLAST_FAILURE::UnknownReason\"\
, \"$command\"));\n	      }\n	    else\n	      {\n\
		my $out;\n		if ($SERVER eq \"NCBI\") {$SERVER=\"\
EBI\"; }\n		elsif ($SERVER eq \"EBI\"){$SERVER=\"N\
CBI\";}\n		add_warning ($$,$$,\"Blast for $name fa\
iled (Run: $run out of $BLAST_MAX_NRUNS. Use $SERV\
ER)\");\n		$out=&run_blast ($name, $method, $db,$i\
nfile, $outfile, $run+1);\n		if ($SERVER eq \"NCBI\
\") {$SERVER=\"EBI\"; }\n		elsif ($SERVER eq \"EBI\
\"){$SERVER=\"NCBI\";}\n		return $out;\n	      }\n\
	  }\n\n	&cache_file(\"SET\",$infile,$name,$method\
,$db,$outfile,$SERVER, $psiJ);\n	return $outfile;\\
n      }\n  }\n\nsub cache_file\n  {\n    my ($cac\
he_mode,$infile,$name,$method,$db, $outfile,$serve\
r,$it)=(@_);\n    my $cache_file;\n    #Protect na\
mes so that they can be turned into legal filename\
s\n    $name=&clean_file_name ($name);\n    if (!$\
it){$it=1;}\n    if ($db=~/\\//)\n      {\n	$db=~/\
([^\\/]+)$/;\n	$db=$1;\n      }\n\n    $cache_file\
=\"$CACHE/$name.$method.$db.$server.$it.tmp\";\n  \
  #print \"Look for $cache_file [$cache_mode][$CAC\
HE] \\n\";\n    if ($infile ne \"\"){$cache_file_i\
nfile=\"$CACHE/$name.$method.$db.$server.$it.infil\
e.tmp\";}\n\n    if ($cache_mode eq \"GET\")\n    \
  {\n	if ($CACHE eq \"\" || $CACHE eq \"no\" || $C\
ACHE eq \"ignore\"  || $CACHE eq \"local\" || $CAC\
HE eq \"update\"){return 0;}\n	elsif ( !-d $CACHE)\
\n	  {\n	    print STDERR \"ERROR: Cache Dir: $CAC\
HE Does not Exist\";\n	    return 0;\n	  }\n	else\\
n	  {\n	    my ($z1,$z1m)=uncompress($cache_file_i\
nfile);\n	    my ($z2,$z2m)=uncompress($cache_file\
);\n	    \n	    if ( -e $cache_file && &fasta_file\
1_eq_fasta_file2($infile,$cache_file_infile)==1)\n\
	      {\n		`cp $cache_file $outfile`;\n		$CACHE_S\
TATUS=\"READ CACHE\";\n		compress($z1,$z1m);\n		co\
mpress($z2,$z2m);\n		\n		return 1;\n	      }\n	  }\
\n      }\n    elsif ($cache_mode eq \"SET\")\n   \
   {\n	if ($CACHE eq \"\" || $CACHE eq \"no\" || $\
CACHE eq \"ignore\"  || $CACHE eq \"local\" || $CA\
CHE eq \"update\"){return 0;}\n	elsif ( !-d $CACHE\
)\n	  {\n	    print STDERR \"ERROR: Cache Dir: $CA\
CHE Does not Exist\";\n	    return 0;\n	  }\n	elsi\
f (-e $outfile)\n	  {\n	    `cp $outfile $cache_fi\
le`;\n	    if ($cache_file_infile ne \"\"){ `cp $i\
nfile $cache_file_infile`;}\n	    return 1;\n	  }\\
n      }\n    $CACHE_STATUS=\"COMPUTE CACHE\";\n  \
  return 0;\n  }\nsub file1_eq_file2\n  {\n    my \
($f1, $f2)=@_;\n    if ( $f1 eq \"\"){return 1;}\n\
    elsif ( $f2 eq \"\"){return 1;}\n    elsif ( !\
-e $f1){return 0;}\n    elsif ( !-e $f2){return 0;\
}\n    elsif ($f1 eq \"\" || $f2 eq \"\" || `diff \
$f1 $f2` eq \"\"){return 1;}\n\n    return 0;\n  }\
\nsub clean_file_name\n  {\n    my $name=@_[0];\n\\
n    $name=~s/[^A-Za-z1-9.-]/_/g;\n    return $nam\
e;\n  }\nsub url2file\n  {\n    my ($address, $out\
)=(@_);\n\n    if (&pg_is_installed (\"wget\"))\n	\
{\n	  return &safe_system (\"wget $address -O$out \
>/dev/null 2>/dev/null\");\n	}\n    elsif (&pg_is_\
installed (\"curl\"))\n      {\n	return &safe_syst\
em (\"curl $address -o$out >/dev/null 2>/dev/null\\
");\n      }\n    else\n      {\n	myexit(flus_erro\
r(\"neither curl nor wget are installed. Imnpossib\
le to fectch remote file\"));\n	exit ($EXIT_FAILUR\
E);\n      }\n  }\nsub fasta_file1_eq_fasta_file2\\
n  {\n    my ($f1, $f2)=@_;\n    my (%s1, %s2);\n \
   my @names;\n    %s1=read_fasta_seq ($f1);\n    \
%s2=read_fasta_seq ($f2);\n\n    @names=(keys (%s1\
));\n\n    foreach $n (keys(%s1))\n      {\n	my $s\
s1=lc($s1{$n}{seq});\n	my $ss2=lc($s2{$n}{seq});\n\
	if ($ss1 ne $ss2){return 0;}\n      }\n    foreac\
h $n (keys(%s2))\n      {\n	my $ss1=lc($s1{$n}{seq\
});\n	my $ss2=lc($s2{$n}{seq});\n	if ($ss1 ne $ss2\
){return 0;}\n      }\n    \n    return 1;\n  }\n\\
n\n\nsub read_template_file\n  {\n    my $pdb_temp\
lates = @_[0];\n    my $tmp=new FileHandle;\n    o\
pen ($tmp, \"<$pdb_templates\");\n    my %temp_h;\\
n    my ($skip,$seq, $temp);\n\n    #supports both\
 a simple [seq pdb] format and the regular fasta l\
ike template format\n    while (<$tmp>)\n      {\n\
	\n	$line = $_;\n	if ($line=~/\\>(\\S+)\\s_._\\s(\\
\S+)/){$temp_h{$1}= $2;}\n	elsif ($line =~/(\\S+)\\
\s(\\S+)/){$temp_h{$1}= $2;}\n      }\n    close($\
tmp);\n    return %temp_h;\n  }\n\n\n\n\n\n\nsub s\
eq2tblastx_lib\n  {\n    my ($mode, $infile, $outf\
ile)=@_;\n    my (%s, $method,$nseq);\n\n    $meth\
od=$mode;\n    &set_temporary_dir (\"set\",$infile\
,\"infile\");\n    %s=read_fasta_seq(\"infile\");\\
n\n\n    foreach $seq (keys(%s))\n      {\n	$slist\
[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$s{$seq}\
{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{order}]=\
length ($s{$seq}{seq});\n      }\n    $nseq=$#snam\
e+1;\n    open (F, \">outfile\");\n    print F \"!\
 TC_LIB_FORMAT_01\\n\";\n    print F \"$nseq\\n\";\
\n    for ($a=0; $a<$nseq;$a++)\n      {\n	print F\
 \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n      }\\
n    close (F);\n    &safe_system (\"formatdb -i i\
nfile -p F\");\n    &safe_system (\"blastall -p tb\
lastx -i infile -d infile -m 7 -S1>blast.output\")\
;\n\n    ncbi_tblastx_xml2lib_file (\"outfile\", f\
ile2string (\"blast.output\"));\n    &set_temporar\
y_dir (\"unset\",$mode, $method, \"outfile\",$outf\
ile);\n    myexit ($EXIT_SUCCESS);\n    }\nsub seq\
2tblastpx_lib\n  {\n    my ($mode, $infile, $outfi\
le)=@_;\n    my (%s, $method,$nseq);\n    $method=\
$mode;\n    &set_temporary_dir (\"set\",$infile,\"\
infile\");\n    %s=read_fasta_seq(\"infile\");\n\n\
    foreach $seq (keys(%s))\n      {\n	$slist[$s{$\
seq}{order}]=$s{$seq}{seq};\n	$sname[$s{$seq}{orde\
r}]=$s{$seq}{name};\n	$slen[$s{$seq}{order}]=lengt\
h ($s{$seq}{seq});\n      }\n    $nseq=$#sname+1;\\
n    open (F, \">outfile\");\n    print F \"! TC_L\
IB_FORMAT_01\\n\";\n    print F \"$nseq\\n\";\n   \
 for ($a=0; $a<$nseq;$a++)\n      {\n	print F \"$s\
name[$a] $slen[$a]  $slist[$a]\\n\"\n      }\n    \
close (F);\n    &safe_system(\"t_coffee -other_pg \
seq_reformat -in infile -output tblastx_db1 > tbla\
stxdb\");\n    &safe_system (\"formatdb -i tblastx\
db -p T\");\n    &safe_system (\"blastall -p blast\
p -i tblastxdb -d tblastxdb -m7 >blast.output\");\\
n    ncbi_tblastpx_xml2lib_file (\"outfile\", file\
2string (\"blast.output\"), %s);\n    &set_tempora\
ry_dir (\"unset\",$mode, $method, \"outfile\",$out\
file);\n    myexit ($EXIT_SUCCESS);\n    }\n\nsub \
x3dna_find_pair2lib\n      {\n      my ($seq, $pdb\
, $lib, $pg)=@_;\n      my $outfile1=\"dssr-2ndstr\
s.dbn\";\n      my $outfile2=\"simple.output\";\n \
     my $f= new FileHandle;\n      my ($rnaSS,$pdb\
SS);\n      my $command;\n      my %s_pdb;\n      \
my %s_seq;\n      my $pdbseq=vtmpnam(NULL);\n     \
 \n      #$pg: \"find_pair\" OR \"find_pair -p\"\n\
      \n      if (!pg_is_installed (\"find_pair\")\
)\n	{\n	  add_warning ($$,$$, \"x3dna/find_pair co\
uld not be used to extract RNA secondary structure\
s. Secondary structures will be extracted by x3dna\
-ssr instead -- Install the find-pair module of x3\
dna  [http://x3dna.org/]\");\n	  return x3dnassr2l\
ib ($seq, $pdb, $lib);\n	}\n      \n\n      #get P\
DB sequence\n      safe_system (\"t_coffee -other_\
pg extract_from_pdb $pdb -seq >$pdbseq\");\n      \
\n      #get find_pair contacts\n      $command=\"\
$pg $pdb simple.output > /dev/null 2>/dev/null\";\\
n      safe_system ($command);\n\n      if (($comm\
and=~/find_pair -p/)){$outfile2=\"allpairs.ana\";}\
\n      else {$outfile2=\"simple.output\";}\n     \
 \n      if ( !-e $outfile2)\n	{\n	  myexit(flush_\
error(\"x3dna failed to compute the secondary stru\
cture file $outfile2 for $pdb\"));\n	  myexit ($EX\
IT_FAILURE);\n	}\n      \n\n      #Handle situatio\
ns when the pdb sequence differs from the RNA sequ\
ence\n      #my @out=file2array($outfile1);\n     \
 %s_pdb=read_fasta_seq_index ($pdbseq);\n      %s_\
seq=read_fasta_seq_index ($seq);\n      my $rnaS=u\
c($s_seq{0}{seq});\n      my $pdbS=uc($s_pdb{0}{se\
q});\n      my $vienna;\n      my @lu;\n    \n    \
  if ($rnaS ne $pdbS)\n	{\n	  \n	  my ($rna,$pdb);\
\n	  $rnaSS=$rnaS;\n	  $pdbSS=$pdbS;\n	  $rnaSS=~s\
/T/U/g;\n	  $pdbSS=~s/T/U/g;\n	  ($rnaSS,$pdbSS)=n\
w ($rnaS, $pdbS);\n	  \n	  my @rnaA =split (//,$rn\
aSS);\n	  my @pdbA=split (//,$pdbSS);\n	  my $l=@r\
naA;\n	  \n	  #print \"\\n--- $s_seq{0}{name} $rna\
SS\\n--- $s_seq{0}{name} $pdbSS\\n\\n\";\n	  \n	  \
for (my $b=0,my $a=0; $a<$l; $a++)\n	    {\n	     \
 if   ($rnaA[$a] ne '-' && $pdbA[$a] ne '-'){$lu[+\
+$pdb]=++$rna;}\n	      elsif($rnaA[$a] eq '-'){$l\
u[++$pdb]=-1;}\n	      elsif($pdbA[$a] eq '-'){++$\
rna;}\n	    }\n	}\n      else\n	{\n	  for (my $a=0\
; $a<=length ($rnaS); $a++)\n	    {\n	      $lu[$a\
]=$a;\n	    }\n	}\n      my $l=length ($rnaS);\n  \
    open ($f, \">$lib\");\n      print $f \"! TC_L\
IB_FORMAT_01\\n\";\n      print $f \"1\\n\";\n    \
  print $f \"$s_seq{0}{name} $l $rnaS\\n\";\n     \
 print $f \"!CMT: [SOURCE] >$s_seq{0}{name} 3D con\
tact library Generated by $pg (x3dna)\\n\";\n     \
 print $f \"#1 1\\n\";\n      \n      my $ne;\n   \
   my @array=file2array($outfile2);\n      for (my\
 $a=0; $a<5; $a++){shift (@array);}\n      while (\
!($array[0]=~/####/))\n	{\n	  my $line= shift (@ar\
ray);\n	  my @l=($line=~/(\\d+)/g);\n	  \n	 \n	  m\
y $f1=$lu[$l[0]];\n	  my $s1=$lu[$l[1]];\n\n	  #pr\
int \"\\n$line\\n$l[0] --> $f1\\n$l[1] --> $s1\\n\\
\n\"; \n	  \n	  if (!$f1 || !$s1)\n	    {\n	      \
print \"\\n1---- $rnaSS\\n2---- $pdbSS\\n$line\\n[\
$l[0] --- $l[1]]<---->[$f1 --- $s1]\\n\";\n	      \
myexit(flush_error(\"Error while parsing s3dna::fi\
nd_pair output\"));\n	    }\n	  elsif ($f1==-1 || \
$s1==-1){;}\n	  elsif ($f1<$s1){print $f \"$f1 $s1\
 100\\n\";}\n	  else {print $f \"$s1 $f1 100\\n\";\
$ne++;}\n	}\n      print $f \"! SEQ_1_TO_N\\n\";\n\
      close ($f);\n      return;\n    }\nsub RNApl\
fold2lib\n  {\n    my ($seq, $lib)=@_;\n    my $f=\
 new FileHandle;\n    \n    &safe_system (\"t_coff\
ee -other_pg RNAplfold2tclib.pl -in=$seq -out=$lib\
\");\n    \n    if ( !-e $lib)\n	{\n	 myexit(flush\
_error(\"RNAplfold failed to compute the secondary\
 structure of $s{$seq}{name}\"));\n	 myexit ($EXIT\
_FAILURE);\n       }\n    open ($f, \">>$lib\");\n\
    print $f \"!CMT: [SOURCE] 2D contact library G\
enerated by RNAPlfold (Vienna Package)\\n\";\n    \
close $f;\n    return;\n  }\nsub x3dnassr2lib\n   \
 {\n      my ($seq, $pdb, $lib)=@_;\n      my $out\
file=\"dssr-2ndstrs.dbn\";\n      my $f= new FileH\
andle;\n      \n\n      if (!pg_is_installed (\"x3\
dna-ssr\"))\n	{\n	  add_warning ($$,$$, \"x3dna-ss\
r could not be used to extract RNA secondary struc\
tures. Secondary structures will be predicted ab-i\
nitio instead with RNAPlfold -- Install s3dna [htt\
p://x3dna.org/] \");\n	  return RNAplfold2lib ($se\
q,$lib);\n	}\n      \n      safe_system (\"x3dna-s\
sr -i=$pdb >/dev/null 2>/dev/null\");\n      if ( \
!-e $outfile)\n	{\n	  myexit(flush_error(\"x3dna-s\
sr failed to compute the secondary structure file \
\"));\n	  myexit ($EXIT_FAILURE);\n	}\n\n      #Ha\
ndle situations when the pdb sequence differs from\
 the RNA sequence\n      @out=file2array($outfile)\
;\n      my %s=read_fasta_seq ($seq);\n      my @n\
ames=keys (%s);\n      my $rnaS=uc($s{$names[0]}{s\
eq});\n      my $pdbS=uc($out[1]);\n      my $vien\
na;\n      \n      #x3dna returns non legitimate n\
ucleotides\n       $pdbS=~s/[^AGCTU]//g;\n      \n\
      if ($rnaS ne $pdbS)\n	{\n	  my ($rna,$pdb);\\
n	  my $rnaSS=$rnaS;\n	  my $pdbSS=$pdbS;\n	  $rna\
SS=~s/T/U/g;\n	  $pdbSS=~s/T/U/g;\n	  ($rnaSS,$pdb\
SS)=nw ($rnaSS, $pdbSS);\n	  my @rnaA =split (//,$\
rnaSS);\n	  my @pdbA=split (//,$pdbSS);\n	  my @SS\
=split (//, $out[2]);\n	  \n	  my $l=@rnaA;\n	  fo\
r (my $b=0,my $a=0; $a<$l; $a++)\n	    {\n	      i\
f   ($rnaA[$a] ne '-' && $pdbA[$a] ne '-'){$vienna\
.=$SS[$b++];}\n	      elsif($rnaA[$a] eq '-'){$b++\
;}\n	      elsif($pdbA[$a] eq '-'){$vienna.='.';}\\
n	    }\n	}\n      else\n	{\n	  $vienna=$out[2];\n\
	}\n    \n\n      open ($f, \">seq\");\n      prin\
t $f \">$names[0]\\n$rnaS\\n\";\n      close $f;\n\
      \n      open ($f, \">str\");\n      print $f\
 \">$names[0]\\n$vienna\\n\";\n      close $f;\n  \
    \n      safe_system (\"t_coffee -other_pg seq_\
reformat -in seq -in2 str -output vienna2tc_lib >$\
lib\");\n      if ( !-e $lib)\n	    {\n	      myex\
it(flush_error(\"seq_reformat failed to convert yo\
ur secondary structure\"));\n	      myexit ($EXIT_\
FAILURE);\n	    }\n      \n      open ($f, \">>$li\
b\");\n      print $f \"!CMT: [SOURCE] >$names[0] \
2D Contact library generated by x3dna-ssr\\n\";\n \
     #print $f \"! Vienna_Format: >$names[0]\\n\";\
\n      #print $f \"! Vienna_Format: $vienna\\n\";\
\n      \n      close $f;\n      return;\n    }\n\\
n\nsub file2head\n      {\n	my $file = shift;\n	my\
 $size = shift;\n	my $f= new FileHandle;\n	my $lin\
e;\n	open ($f,$file);\n	read ($f,$line, $size);\n	\
close ($f);\n	return $line;\n      }\nsub file2tai\
l\n      {\n	my $file = shift;\n	my $size = shift;\
\n	my $f= new FileHandle;\n	my $line;\n\n	open ($f\
,$file);\n	seek ($f,$size*-1, 2);\n	read ($f,$line\
, $size);\n	close ($f);\n	return $line;\n      }\n\
\n\nsub vtmpnam\n      {\n	my $r=rand(100000);\n	m\
y $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n	    $\
f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f);\n	r\
eturn $f;\n      }\n\nsub myexit\n  {\n    my $cod\
e=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){return;}\
\n    else {$CLEAN_EXIT_STARTED=1;}\n    ### ONLY \
BARE EXIT\n    exit ($code);\n  }\nsub set_error_l\
ock\n    {\n      my $name = shift;\n      my $pid\
=$$;\n\n\n      &lock4tc ($$,\"LERROR\", \"LSET\",\
 \"$$ -- ERROR: $name $PROGRAM\\n\");\n      retur\
n;\n    }\nsub set_lock\n  {\n    my $pid=shift;\n\
    my $msg= shift;\n    my $p=getppid();\n    &lo\
ck4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$msg\\n\");\n\
  }\nsub unset_lock\n   {\n\n    my $pid=shift;\n \
   &lock4tc ($pid,\"LLOCK\",\"LRELEASE\",\"\");\n \
 }\nsub shift_lock\n  {\n    my $from=shift;\n    \
my $to=shift;\n    my $from_type=shift;\n    my $t\
o_type=shift;\n    my $action=shift;\n    my $msg;\
\n\n    if (!&lock4tc($from, $from_type, \"LCHECK\\
", \"\")){return 0;}\n    $msg=&lock4tc ($from, $f\
rom_type, \"LREAD\", \"\");\n    &lock4tc ($from, \
$from_type,\"LRELEASE\", $msg);\n    &lock4tc ($to\
, $to_type, $action, $msg);\n    return;\n  }\nsub\
 isshellpid\n  {\n    my $p=shift;\n    if (!lock4\
tc ($p, \"LLOCK\", \"LCHECK\")){return 0;}\n    el\
se\n      {\n	my $c=lock4tc($p, \"LLOCK\", \"LREAD\
\");\n	if ( $c=~/-SHELL-/){return 1;}\n      }\n  \
  return 0;\n  }\nsub isrootpid\n  {\n    if(lock4\
tc (getppid(), \"LLOCK\", \"LCHECK\")){return 0;}\\
n    else {return 1;}\n  }\nsub lock4tc\n	{\n	  my\
 ($pid,$type,$action,$value)=@_;\n	  my $fname;\n	\
  my $host=hostname;\n\n	  if ($type eq \"LLOCK\")\
{$fname=\"$LOCKDIR/.$pid.$host.lock4tcoffee\";}\n	\
  elsif ( $type eq \"LERROR\"){ $fname=\"$LOCKDIR/\
.$pid.$host.error4tcoffee\";}\n	  elsif ( $type eq\
 \"LWARNING\"){ $fname=\"$LOCKDIR/.$pid.$host.warn\
ing4tcoffee\";}\n\n	  if ($debug_lock)\n	    {\n	 \
     print STDERR \"\\n\\t---lock4tc(tcg): $action\
 => $fname =>$value (RD: $LOCKDIR)\\n\";\n	    }\n\
\n	  if    ($action eq \"LCHECK\") {return -e $fna\
me;}\n	  elsif ($action eq \"LREAD\"){return file2\
string($fname);}\n	  elsif ($action eq \"LSET\") {\
return string2file ($value, $fname, \">>\");}\n	  \
elsif ($action eq \"LRESET\") {return string2file \
($value, $fname, \">\");}\n	  elsif ($action eq \"\
LRELEASE\")\n	    {\n	      if ( $debug_lock)\n		{\
\n		  my $g=new FileHandle;\n		  open ($g, \">>$fn\
ame\");\n		  print $g \"\\nDestroyed by $$\\n\";\n\
		  close ($g);\n		  safe_system (\"mv $fname $fna\
me.old\");\n		}\n	      else\n		{\n		  unlink ($fn\
ame);\n		}\n	    }\n	  return \"\";\n	}\n\nsub fil\
e2string\n	{\n	  my $file=@_[0];\n	  my $f=new Fil\
eHandle;\n	  my $r;\n	  open ($f, \"$file\");\n	  \
while (<$f>){$r.=$_;}\n	  close ($f);\n	  return $\
r;\n	}\nsub file2array\n	{\n	  my $file=@_[0];\n	 \
 my $f=new FileHandle;\n	  my @r;\n	  open ($f, \"\
$file\");\n	  while (<$f>){@r=(@r,$_);}\n	  close \
($f);\n	  return @r;\n	}\nsub string2file\n    {\n\
    my ($s,$file,$mode)=@_;\n    my $f=new FileHan\
dle;\n\n    open ($f, \"$mode$file\");\n    print \
$f  \"$s\";\n    close ($f);\n  }\n\nBEGIN\n    {\\
n      srand;\n\n      $SIG{'SIGUP'}='signal_clean\
up';\n      $SIG{'SIGINT'}='signal_cleanup';\n    \
  $SIG{'SIGQUIT'}='signal_cleanup';\n      $SIG{'S\
IGILL'}='signal_cleanup';\n      $SIG{'SIGTRAP'}='\
signal_cleanup';\n      $SIG{'SIGABRT'}='signal_cl\
eanup';\n      $SIG{'SIGEMT'}='signal_cleanup';\n \
     $SIG{'SIGFPE'}='signal_cleanup';\n\n      $SI\
G{'SIGKILL'}='signal_cleanup';\n      $SIG{'SIGPIP\
E'}='signal_cleanup';\n      $SIG{'SIGSTOP'}='sign\
al_cleanup';\n      $SIG{'SIGTTIN'}='signal_cleanu\
p';\n      $SIG{'SIGXFSZ'}='signal_cleanup';\n    \
  $SIG{'SIGINFO'}='signal_cleanup';\n\n      $SIG{\
'SIGBUS'}='signal_cleanup';\n      $SIG{'SIGALRM'}\
='signal_cleanup';\n      $SIG{'SIGTSTP'}='signal_\
cleanup';\n      $SIG{'SIGTTOU'}='signal_cleanup';\
\n      $SIG{'SIGVTALRM'}='signal_cleanup';\n     \
 $SIG{'SIGUSR1'}='signal_cleanup';\n\n\n      $SIG\
{'SIGSEGV'}='signal_cleanup';\n      $SIG{'SIGTERM\
'}='signal_cleanup';\n      $SIG{'SIGCONT'}='signa\
l_cleanup';\n      $SIG{'SIGIO'}='signal_cleanup';\
\n      $SIG{'SIGPROF'}='signal_cleanup';\n      $\
SIG{'SIGUSR2'}='signal_cleanup';\n\n      $SIG{'SI\
GSYS'}='signal_cleanup';\n      $SIG{'SIGURG'}='si\
gnal_cleanup';\n      $SIG{'SIGCHLD'}='signal_clea\
nup';\n      $SIG{'SIGXCPU'}='signal_cleanup';\n  \
    $SIG{'SIGWINCH'}='signal_cleanup';\n\n      $S\
IG{'INT'}='signal_cleanup';\n      $SIG{'TERM'}='s\
ignal_cleanup';\n      $SIG{'KILL'}='signal_cleanu\
p';\n      $SIG{'QUIT'}='signal_cleanup';\n\n     \
 our $debug_lock=$ENV{\"DEBUG_LOCK\"};\n\n\n\n\n  \
    foreach my $a (@ARGV){$CL.=\" $a\";}\n      if\
 ( $debug_lock ){print STDERR \"\\n\\n\\n*********\
* START PG: $PROGRAM *************\\n\";}\n      i\
f ( $debug_lock ){print STDERR \"\\n\\n\\n********\
**(tcg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\\
n      if ( $debug_lock ){print STDERR \"\\n --- $\
$ -- $CL\\n\";}\n\n\n\n\n    }\nsub flush_error\n \
 {\n    my $msg=shift;\n    $msg.=\" [tc_generic_m\
ethod.pl/FATAL]\";\n    return add_error ($EXIT_FA\
ILURE,$$, $$,getppid(), $msg, $CL);\n  }\nsub add_\
error\n  {\n    my $code=shift;\n    my $rpid=shif\
t;\n    my $pid=shift;\n    my $ppid=shift;\n    m\
y $type=shift;\n    my $com=shift;\n\n    $ERROR_D\
ONE=1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\"\
$pid -- ERROR: $type\\n\");\n    lock4tc ($$, \"LE\
RROR\",\"LSET\", \"$pid -- COM: $com\\n\");\n    l\
ock4tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK: \
$ppid -> $pid\\n\");\n\n    return $code;\n  }\nsu\
b add_warning\n  {\n    my $rpid=shift;\n    my $p\
id =shift;\n    my $command=shift;\n    my $msg=\"\
$$ -- WARNING: $command\\n\";\n    print STDERR \"\
$msg\";\n    lock4tc ($$, \"LWARNING\", \"LSET\", \
$msg);\n  }\n\nsub signal_cleanup\n  {\n    print \
dtderr \"\\n**** $$ (tcg) was killed\\n\";\n    &c\
leanup;\n    exit ($EXIT_FAILURE);\n  }\nsub clean\
_dir\n  {\n    my $dir=@_[0];\n    if ( !-d $dir){\
return ;}\n    elsif (!($dir=~/tmp/)){return ;}#sa\
fety check 1\n    elsif (($dir=~/\\*/)){return ;}#\
safety check 2\n    else\n      {\n	`rm -rf $dir`;\
\n      }\n    return;\n  }\nsub cleanup\n  {\n   \
 #print stderr \"\\n----tc: $$ Kills $PIDCHILD\\n\\
";\n    #kill (SIGTERM,$PIDCHILD);\n    my $p=getp\
pid();\n    $CLEAN_EXIT_STARTED=1;\n\n\n\n    if (\
&lock4tc($$,\"LERROR\", \"LCHECK\", \"\"))\n      \
{\n	my $ppid=getppid();\n	if (!$ERROR_DONE)\n	  {\\
n	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- ST\
ACK: $p -> $$\\n\");\n	    &lock4tc($$,\"LERROR\",\
 \"LSET\", \"$$ -- COM: $CL\\n\");\n	  }\n      }\\
n    my $warning=&lock4tc($$, \"LWARNING\", \"LREA\
D\", \"\");\n    my $error=&lock4tc($$,  \"LERROR\\
", \"LREAD\", \"\");\n    #release error and warni\
ng lock if root\n\n    if (isrootpid() && ($warnin\
g || $error) )\n      {\n\n	print STDERR \"*******\
********* Summary *************\\n$error\\n$warnin\
g\\n\";\n\n	&lock4tc($$,\"LERROR\",\"RELEASE\",\"\\
");\n	&lock4tc($$,\"LWARNING\",\"RELEASE\",\"\");\\
n      }\n\n\n    foreach my $f (@TMPFILE_LIST)\n \
     {\n	if (-e $f){unlink ($f);}\n      }\n    fo\
reach my $d (@TMPDIR_LIST)\n      {\n	clean_dir ($\
d);\n      }\n    #No More Lock Release\n    #&loc\
k4tc($$,\"LLOCK\",\"LRELEASE\",\"\"); #release loc\
k\n\n    if ( $debug_lock ){print STDERR \"\\n\\n\\
\n********** END PG: $PROGRAM ($$) *************\\\
n\";}\n    if ( $debug_lock ){print STDERR \"\\n\\\
n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ *********\
****\\n\";}\n  }\nEND\n  {\n\n    &cleanup();\n  }\
\n\nsub blast_com2new_blast_com\n    {\n      my $\
com=shift;\n      if ($com=~/formatdb/)\n	{\n	  $c\
om=~s/formatdb/makeblastdb/;\n	  $com=~s/\\-i/\\-i\
n/;\n	  if ($com =~/pF/){$com=~s/\\-pF/\\-dbtype n\
ucl/;}\n	  if ($com =~/p F/){$com=~s/\\-p F/\\-dbt\
ype nucl/;}\n	  $com=\"$com -logfile /dev/null\";\\
n	  return $com;\n	}\n      else {return $com;}\n \
     \n    }\nsub safe_system\n{\n  my $com=shift;\
\n  my $ntry=shift;\n  my $ctry=shift;\n  my $pid;\
\n  my $status;\n  my $ppid=getppid();\n  if ($com\
 eq \"\"){return 1;}\n \n  if ( ($com=~/^blast/) |\
|($com=~/^formatdb/)){$com=&blast_com2new_blast_co\
m($com);}\n \n  if (($pid = fork ()) < 0){return (\
-1);}\n  if ($pid == 0)\n    {\n      set_lock($$,\
 \" -SHELL- $com (tcg)\");\n      if( $debug_gener\
ic_method ) { printf \"~ exec: %s\\n\", $com; }\n \
     exec ($com);\n      if( $debug_generic_method\
 ) { printf \"~ exitcode: %s\\n\", $?; }\n    }\n \
 else\n    {\n      lock4tc ($$, \"LLOCK\", \"LSET\
\", \"$pid\\n\");#update parent\n      $PIDCHILD=$\
pid;\n    }\n  if ($debug_lock){printf STDERR \"\\\
n\\t .... safe_system (fasta_seq2hmm)  p: $$ c: $p\
id COM: $com\\n\";}\n\n  waitpid ($pid,WTERMSIG);\\
n\n  shift_lock ($pid,$$, \"LWARNING\",\"LWARNING\\
", \"LSET\");\n\n  if ($? == $EXIT_FAILURE || lock\
4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\n    {\n \
     if ($ntry && $ctry <$ntry)\n	{\n\n	  add_warn\
ing ($$,$$,\"$com failed [retry: $ctry out of $ntr\
y]\");\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\"\
, \"\");\n	  #if ($com=~/EBI/){$com=~s/EBI/NCBI/;}\
\n	  #elsif ($com=~/NCBI/){$com=~s/NCBI/EBI/;}\n\n\
	  return safe_system ($com, $ntry, ++$ctry);\n	}\\
n      elsif ($ntry == -1)\n	{\n	  if (!shift_lock\
 ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"))\n\
	    {\n	      add_warning ($$,$$,\"$com failed\")\
;\n	    }\n	  else\n	    {\n	      lock4tc ($pid, \
\"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  retur\
n $?;}\n      else\n	{\n	  if (!shift_lock ($pid,$\
$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n	   \
   myexit(add_error ($EXIT_FAILURE,$$,$pid,getppid\
(), \"UNSPECIFIED system\", $com));\n	    }\n	}\n \
   }\n  return $?;\n}\n\nsub check_configuration\n\
    {\n      my @l=@_;\n      my $v;\n      foreac\
h my $p (@l)\n	{\n\n	  if   ( $p eq \"EMAIL\")\n	 \
   {\n	      if ( !($EMAIL=~/@/))\n		{\n		add_warn\
ing($$,$$,\"Could Not Use EMAIL\");\n		myexit(add_\
error ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\",\"$\
CL\"));\n	      }\n	    }\n	  elsif( $p eq \"INTER\
NET\")\n	    {\n	      if ( !&check_internet_conne\
ction())\n		{\n		  myexit(add_error ($EXIT_FAILURE\
,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n		}\n	  \
  }\n	  elsif( $p eq \"wget\")\n	    {\n	      if \
(!&pg_is_installed (\"wget\") && !&pg_is_installed\
 (\"curl\"))\n		{\n		  myexit(add_error ($EXIT_FAI\
LURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\",\"$\
CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_installed\
 ($p)))\n	    {\n	      myexit(add_error ($EXIT_FA\
ILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\"$C\
L\"));\n	    }\n	}\n      return 1;\n    }\nsub nw\
\n      {\n	my($A,$B)=@_;\n	my ($i,$j, $s);\n	my $\
gep=-2;\n	my $match=+2;\n	my $mismatch=0;\n	my ($s\
ub, $ins, $del);\n\n\n	if ($A eq $B){return ($A,$B\
);}\n	\n	$A=~s/[\\s\\d]//g;	\n	$B=~s/[\\s\\d]//g;	\
\n\n\n	my @rA=split ('',$A);\n	my @rB=split ('',$B\
);\n	\n	#evaluate substitutions\n	my $lenA=@rA;\n	\
my $lenB=@rB;\n	\n	for ($i=0; $i<=$lenA; $i++){$sm\
at[$i][0]=$i*$gep;$tb[$i][0 ]= 1;}\n	for ($j=0; $j\
<=$lenB; $j++){$smat[0][$j]=$j*$gep;$tb[0 ][$j]=-1\
;}\n	\n	for ($i=1; $i<=$lenA; $i++)\n	  {\n	    fo\
r ($j=1; $j<=$lenB; $j++)\n	      {\n		if ($rA[$i-\
1] eq $rB[$j-1]){$s=$match;}\n		else {$s=$mismatch\
;}\n		\n		$sub=$smat[$i-1][$j-1]+$s;\n		$del=$smat\
[$i  ][$j-1]+$gep;\n		$ins=$smat[$i-1][$j  ]+$gep;\
\n		\n		if   ($sub>=$del && $sub>=$ins){$smat[$i][\
$j]=$sub;$tb[$i][$j]=0;}\n		elsif($del>$ins){$smat\
[$i][$j]=$del;$tb[$i][$j]=-1;}\n		else {$smat[$i][\
$j]=$ins;$tb[$i][$j]=1;}\n		}\n	  }\n	#print \"\\n\
---- SCORE=$smat[$lenA][$lenB]\\n\";\n	\n	$i=$lenA\
;\n	$j=$lenB;\n	my $aln_len=0;\n\n	while (!($i==0 \
&& $j==0))\n	  {\n	    if ($tb[$i][$j]==0)\n	    {\
\n	      $aA[$aln_len]=$rA[--$i];\n	      $aB[$aln\
_len]=$rB[--$j];\n	    }\n	  elsif ($tb[$i][$j]==-\
1)\n	    {\n	      $aA[$aln_len]='-';\n	      $aB[\
$aln_len]=$rB[--$j];\n	    }\n	  elsif ($tb[$i][$j\
]==1)\n	    {\n	      $aA[$aln_len]=$rA[--$i];\n	 \
     $aB[$aln_len]='-';\n	    }\n	  $aln_len++;\n	\
  }\n	\n	\n	@aA=reverse (@aA);\n	@aB=reverse (@aB)\
;\n	my $sA=join('',@aA);\n	my $sB=join('',@aB);\n	\
return ($sA,$sB);\n      }\n      \nsub fasta2nseq\
\n	{\n	  \n	  my $f=@_[0];\n	  my $nseq;\n\n	  ope\
n (F, \"$f\") or return 0;\n	  while (<F>)\n	    {\
\n	      if ($_=~/\\>/){$nseq++;}\n	    }\n	  clos\
e (F);\n	  return $nseq;\n	}\n	\nsub compress\n	  \
{\n	    my ($f, $mode)=@_;\n	    if    ($mode eq \\
"gz\"){system (\"gzip $f\");}\n	    elsif ($mode e\
q \"zip\" ){system (\"zip $f\");}\n	    return;\n	\
  }\nsub uncompress \n	  {\n	    my $f=@_[0];\n	  \
  if ( -e $f) {return \"\";}\n\n	    my $gz=$f.\".\
gz\";\n	    if ( -e $gz)\n	      {\n		system (\"gu\
nzip $gz\");\n		return ($f, \"gz\");\n	      }\n	 \
   my $gz=$f.\".zip\";\n	    \n	    if ( -e $gz)\n\
	      {\n		system (\"unzip $gz\");\n		return ($f,\
 \"zip\");\n	      }\n	    return \"\";\n	  }\nmy \
$program=\"T-COFFEE (Version_13.45.15.8d1d86e)\";\\
n","use Env;\nuse strict;\nuse FileHandle;\nuse Cw\
d;\nuse File::Path;\nuse Sys::Hostname;\nuse File:\
:Temp qw/ tempfile tempdir /;\nmy $QUIET=\"2>/dev/\
null\";\nmy $VERBOSE=$ENV{VERBOSE_4_DYNAMIC};\nour\
 $EXIT_FAILURE=1;\nour $EXIT_SUCCESS=0;\nour $LAST\
_COM=\"\";\n\nmy %method;\nmy $method2use;\nmy $tr\
eeF;\nmy $tree=$ENV{\"child_tree_4_TCOFFEE\"};\nmy\
 $thread=$ENV{\"child_thread_4_TCOFFEE\"};\nmy $dy\
namic=$ENV{dynamic_config_4_TCOFFEE};\nmy $clean;\\
nmy $treeFlag;\nmy $blastFlag;\nmy $infile;\nmy $o\
utfile;\nmy $flush;\nmy $do_exit=0;\nmy ($h1, $h2)\
;\nmy @tmpL;\nmy $tmpdir = File::Temp->newdir();\n\
my $stderrF=\"$tmpdir/stderr\";\n$QUIET=\"2>$stder\
rF\";\nmy $CDIR=getcwd();\nmy $threadFlag4tc;\nmy \
$threadFlag4famsa;\nmy $threadFlag;\nmy $tcarg;\n\\
nmy $QUIET_ENV=$ENV{QUIET_ENV};\n\nif ($QUIET_ENV=\
=1){$QUIET=\"\";}\n\nfor ($a=0; $a<=$#ARGV; $a++)\\
n  {\n    if    ($ARGV[$a] eq \"-seq\"){$infile=fi\
le2abs($ARGV[++$a]);}\n    elsif ($ARGV[$a] eq \"-\
outfile\"){$outfile=file2abs($ARGV[++$a], \"new\")\
;}\n    elsif ($ARGV[$a] eq \"-dynamic_config\"){\\
n    	$dynamic=file2abs($ARGV[++$a]);\n    	if ($V\
ERBOSE){print \"\\n![dynamic.pl] --- -dynamic_conf\
ig flag if/else--- $dynamic\\n\";}\n	}\n    \n    \
elsif ($ARGV[$a] eq \"-tree\") {$tree=$ARGV[++$a];\
}\n    elsif ($ARGV[$a] eq \"-method\") {$method2u\
se=$ARGV[++$a];}\n    elsif ($ARGV[$a] eq \"-verbo\
se\"){$VERBOSE=1; $QUIET=\"\";}\n    elsif ($ARGV[\
$a] eq \"-clean\"){$clean=1;}\n  \n    elsif ($ARG\
V[$a] eq \"-thread\"){$thread=$ARGV[++$a]}\n    el\
sif ($ARGV[$a] eq \"-tcarg\") {$tcarg=file2string(\
$ARGV[++$a]);}\n    else \n      {\n	add2tcenv($a+\
+,@ARGV);\n      }\n  }\n\n\n\nif ($tree eq \"list\
\")\n  {\n    my $f=\"$tmpdir/f\";\n    open (F, \\
">$f\");\n    print F \">a\\nxxx\\n>b\\nyyyyy\\n\"\
;\n    close (F);\n    print STDOUT (\"**** Suppor\
ted Guide tree modes:\\n\");\n    my_system (\"t_c\
offee -other_pg seq_reformat -in $f -action +seq2d\
nd list \");\n    $do_exit=1;\n  }\nif ($method2us\
e eq \"list\")\n  {\n    my %ml;\n    my $listfile\
=\"$tmpdir/list\";\n    \n    $ml{tcoffee}=1;\n   \
 $ml{psicoffee}=1;\n    $ml{accurate}=1;\n    $ml{\
'3dcoffee'}=1;\n    $ml{expresso}=1;\n    $ml{clus\
talo}=1;\n    $ml{mafft}=1;\n    $ml{famsa}=1;\n  \
    print STDOUT (\"**** Supported MSA mode:\\n\")\
;\n    my_system (\"t_coffee 2>/dev/null | grep _m\
sa > $listfile\");\n    open (F, $listfile);\n    \
while (<F>)\n      {\n	my $l=$_;\n	$l=~/(.*_msa)\\\
s+(.*)/;\n	my $m=$1;\n	my $i=\"$2\\n\";\n	if ($m=~\
/mafftsparsescore/)\n	  {\n	   printf STDOUT \"%-2\
0s DOES NOT Support [-tree] -- $i\", $m;\n	  }\n	e\
lsif ($m=~/tcoffee/){;}\n	elsif ($m=~/mafft/){;}\n\
	elsif (!$ml{$m})\n	  {\n	    printf STDOUT \"%-20\
s DOES     Support [-tree] -- $i\", $m;\n	  }\n   \
   }\n    $do_exit=1;\n  }\nif ($do_exit){my_exit \
($CDIR,$EXIT_SUCCESS);}\n\nmy $NSEQ=file2nseq($inf\
ile);\n\n\nif ($NSEQ==0)\n  {\n    print \"ERROR -\
 No sequences provided [FATAL:dynamic.pl]\\n\";\n \
   exit ($EXIT_FAILURE);\n  }\nif (!$outfile)\n  {\
\n    ($h1,$outfile)=tempfile();\n    push (@tmpL,\
$outfile);\n    $flush=1;\n  }\n\nif (!($method2us\
e=~/dynamic/)){;}\nelse \n  {\n    if (-e $dynamic\
)\n      {\n       if ($VERBOSE){print \"\\n![dyna\
mic.pl] --- -dynamic_config FILE: \\n\";}\n       \
 my @dynamicFile;\n	open (F, $dynamic);\n	while (<\
F>)\n	  {\n	    my $f=$_;\n	    if ($VERBOSE){prin\
t \"\\n![dynamic.pl] --- FILE content: $f\\n\";}\n\
	    ## $f=~/(\\W)+ (\\d)+/;\n	    @dynamicFile = \
split ' ', $f;\n	    if ($VERBOSE){print \"\\n![dy\
namic.pl] --- -dynamic_config --- $dynamicFile[0] \
:: $dynamicFile[1]\\n\";}\n	    $method{$dynamicFi\
le[0]} = $dynamicFile[1];\n	  }\n	close(F);\n     \
 }\n    else\n      {\n	$method{\"psicoffee_msa\"}\
=50;\n	$method{\"famsa_msa\"}=1000000000;\n      }\
\n    \n    foreach my $name (sort { $method{$a} <\
=> $method{$b} } keys %method) \n      {\n	if ($NS\
EQ<=$method{$name})\n	  {\n	    $method2use=$name;\
\n	    last;\n	  }\n      }\n  }\n\n\nif ($tree)\n\
  {\n    ($h2,$treeF)=tempfile();\n    my ($h2,$tm\
ptree)=tempfile();\n    push (@tmpL,$treeF);\n    \
if ( $tree eq \"default\"){$treeF=0;}\n    elsif (\
 -e $tree)\n      {\n	my_system (\"cp $tree $tmptr\
ee\");\n      }\n    elsif ($tree eq \"master\" ||\
 $tree eq \"main\" || $tree eq \"parent\")\n      \
{\n	if ($ENV{CHILD_TREEF_4_TCOFFEE} && -e $ENV{CHI\
LD_TREEF_4_TCOFFEE})\n	  {\n	    my $ctree=$ENV{CH\
ILD_TREEF_4_TCOFFEE};\n	    my_system (\"mv $ctree\
 $tmptree\");\n	  }\n	else\n	  {\n	    my $master_\
tree=$ENV{CHILD_TREE_4_TCOFFEE};\n	    my_system (\
\"t_coffee -other_pg seq_reformat -in $master_tree\
 -in2 $infile -action +prune_tree -output newick >\
 $tmptree\");\n	  }\n      }\n    else \n      {\n\
	my_system (\"t_coffee -other_pg seq_reformat -in \
$infile -action +seq2dnd $tree -output newick> $tm\
ptree\");\n      }\n    \n    if ($method2use=~/ma\
fft/)\n      {\n	#print \"cp $tmptree /Users/cnotr\
edame/.Trash/$$.tmptree\\n\";\n	#system (\"cp $tmp\
tree /Users/cnotredame/.Trash/$$.tmptree\");\n	\n	\
my_system (\"t_coffee -other_pg seq_reformat -in $\
tmptree -output mafftdndmatrix> $treeF\");\n      \
}\n    else\n      {\n	my_system (\"mv $tmptree $t\
reeF\");\n      }\n  }\nchdir ($tmpdir);\n\nmy $CL\
4tc=get_cl4tc();#will collect from env every CLTCO\
FEE env variable\n\nif (!$treeF || $NSEQ<=2){$tree\
Flag=\"\";}\nelsif ( $method2use=~/coffee/ || $met\
hod2use=~/accurate/){$treeFlag=\"-usetree $treeF \\
";}\nelsif ( $method2use=~/clustalo/){$treeFlag=\"\
--guidetree-in=$treeF \";}\nelsif ( $method2use=~/\
mafftsparsecore/){;}\nelsif ( $method2use=~/mafft/\
){$treeFlag=\"--treein $treeF \";}\nelsif ( $metho\
d2use=~/famsa/){$treeFlag=\"-gt import $treeF \";}\
\n$CL4tc.=\" $treeFlag \";\n\n$threadFlag=($thread\
)?\"--thread $thread \":\"--thread 1 \";\n$threadF\
lag4tc=($thread)?\"-thread $thread \":\"-thread 1 \
\";\n$threadFlag4famsa=($thread)?\"-t $thread \":\\
"-t 1 \";\n$CL4tc.=\" $threadFlag4tc \";\n\nprint \
\"\\n![dynamic.pl] --- CL4tc == $CL4tc\\n\";  \n\n\
my $cmethod=$method2use;\n$cmethod=~s/_pair/_msa/;\
\n$cmethod=~s/_msa//;\n\nprint \"\\n![dynamic.pl] \
--- cmethod == $cmethod\\n\";  \n\nif ($cmethod eq\
 \"tcoffee\"|| $cmethod eq \"t_coffee\" )\n  {\n  \
  my_system (\"t_coffee -seq $infile -outfile $out\
file -output fasta_aln $CL4tc>/dev/null  $QUIET\")\
;    \n  }\nelsif ($cmethod=~/(.*coffee)/ || $cmet\
hod=~/(accurate)/ || $cmethod=~/(expresso)/)\n  {\\
n    my $mode=$1;\n    my_system (\"t_coffee  -mod\
e $mode -seq $infile -outfile $outfile -output fas\
ta_aln $CL4tc >/dev/null  $QUIET\");    \n  }\nels\
if ($cmethod eq \"clustalo\")\n  {\n    my_system \
(\"clustalo -i $infile $treeFlag -o $outfile  --fo\
rce $threadFlag $QUIET\");\n    }\nelsif ($cmethod\
 =~/sparsecore/)\n  {\n    my_system (\"mafft-spar\
secore.rb -i $infile > $outfile $QUIET\");\n  }\ne\
lsif (($cmethod =~/mafft/))\n  {\n    my $mm;\n   \
 my $retree;\n    \n    if ( $cmethod eq \"mafft\"\
 || $cmethod=~/\\-/ )\n      {\n	$mm=$cmethod;\n  \
    }\n    elsif (($cmethod=~/mafft(.*)/))\n      \
{\n	$mm=\"mafft-\".$1;\n      }\n\n    if ($mm =~/\
1/)\n      {\n	$mm=~s/1/i/;\n	$retree=\"--retree 1\
 \"\n      };\n    \n    my_system (\"$mm --anysym\
bol $threadFlag $treeFlag $retree $infile > $outfi\
le $QUIET\");\n  }\n\nelsif ($method2use=~/famsa/)\
\n  {\n    print \"\\n![dynamic.pl] --- FAMSA DEFA\
ULT\\n\";    \n    my_system (\"famsa $treeFlag $t\
hreadFlag4famsa $infile $outfile >/dev/null $QUIET\
\");\n  }\nelsif ($method2use=~/famsaUpgma/)\n  {\\
n    print \"\\n![dynamic.pl] --- FAMSA Upgma\\n\"\
;\n    print \"\\n![dynamic.pl] --- Command: famsa\
 -gt upgma $treeFlag $threadFlag4famsa $infile $ou\
tfile >/dev/null $QUIET\\n\";\n    my_system (\"fa\
msa -gt upgma $treeFlag $threadFlag4famsa $infile \
$outfile >/dev/null $QUIET\");\n  }\nelse\n  {\n  \
  if ($treeF)\n      {\n	printf (STDERR \"WARNING:\
 Method $method2use CANNOT use pre-sepecified guid\
e tree [dynamic.pl]\\n\");\n      }\n    my_system\
 (\"t_coffee -in $infile -method $method2use -outf\
ile $outfile -output fasta_aln $tcarg -quiet $QUIE\
T\");\n  }\n\nif ( ! -e $outfile)\n  {\n    print \
\"ERROR - No MSA computed - $LAST_COM -- [FATAL:dy\
namic.pl]\\n\";\n    my_exit ($CDIR,$EXIT_FAILURE)\
;\n  }\nelsif ( $flush)\n {\n   open (F, $outfile)\
;\n   while (<F>){print $_;}\n   close (F);\n }\nf\
oreach my $f (@tmpL){unlink($f);}\n\n\nif ($VERBOS\
E!=-1)\n  {\n    open (F, \"$stderrF\");\n    whil\
e (<F>)\n      {\n	my $l=$_;\n	if ( $VERBOSE || $l\
=~/WARNING/ || $l=~/ERROR/ || $l=~/INFORNATION/){p\
rint stderr \"$l\";}\n      }\n    close (F);\n  }\
\n\nmy_exit ($CDIR,$EXIT_SUCCESS);\n\n\nsub file2n\
seq\n  {\n    my ($f)=@_;\n    my $n=`grep -c '>' \
$f`; \n    \n    return $n;\n  }\nsub file2abs\n  \
  {\n      my ($f, $mode)=@_;\n      \n      if (!\
$f || $f=~/^\\//){return $f;}\n      elsif (!-e $f\
 && $mode eq \"new\"){return \"$CDIR/$f\";}\n     \
 elsif (!-e $f){return $f;}\n    \n      return \"\
$CDIR/$f\";\n    }\nsub file2string \n    {\n     \
 my ($f)=@_;\n      my $s;\n\n      open (F, $f) |\
| return 0;\n      while (<F>)\n	{\n	  $s.=$_;\n	}\
\n      close (F);\n      chomp($s);\n      return\
 $s;\n    }   \n\nsub get_psicl\n      {\n	my ($ps\
itrim, $psitrim_mode, $pisN);\n	my $cl;\n	\n	if ($\
ENV{psitrim_tree_4_TCOFFEE}){$cl.=\" -psitrim_tree\
=\".$ENV{psitrim_tree_4_TCOFFEE}.\" \";}\n	if ($EN\
V{psitrim_mode_4_TCOFFEE}){$cl.=\" -psitrim_mode=\\
".$ENV{psitrim_mode_4_TCOFFEE}.\" \";}\n	if ($ENV{\
psitrim_4_TCOFFEE}){$cl.=\" -psitrim=\".$ENV{psitr\
im_4_TCOFFEE}.\" \";}\n	if ($ENV{psiJ_4_TCOFFEE}){\
$cl.=\" -psiJ=\".$ENV{psiJ_4_TCOFFEE}.\" \";}\n	\n\
\n	return $cl;\n      }\n      \nsub get_cl4tc\n	{\
\n	  my $cl;\n	  \n	  foreach my $arg (keys(%ENV))\
\n	    {\n	      if ($arg=~/(.*)_4_CLTCOFFEE/)\n		\
{\n		  my $name=$1;\n		  my $val=$ENV{$arg};\n		  \
if (-e $val){$val=file2abs($val);}\n		  \n\n		  if\
 ($val eq \"FLAGSET\"){$val=\"\";}\n		  $cl.=\"-$n\
ame $val \";\n		}\n	    }\n	  if ($VERBOSE){print \
\"\\n![dynamic.pl] --- get_psicl --- $cl\\n\";}\n	\
  return $cl;\n	}\n\nsub add2tcenv\n	    {\n	     \
 my ($p, @argv)=@_;\n\n	      my $flag=$argv[$p];\\
n	      $flag =~s/^-//;\n	      my $val =file2abs(\
$argv[$p+1]);\n	      my $envv=\"$flag\\_4_CLTCOFF\
EE\";\n	      $ENV{$envv}=$val;\n	    }\n	      \n\
sub my_exit\n    {\n      my ($dir,$ec)=@_;\n     \
 my $a;\n      if ($VERBOSE)\n	{\n	  print \"\\n![\
dynamic.pl] --- CDIR: $CDIR\\n\";\n	  print \"\\n!\
[dynamic.pl] --- Processed $NSEQ\\n\";\n	  print \\
"\\n![dynamic.pl] --- \";\n	  foreach my $arg (@AR\
GV)\n	    {\n	      print \"$arg \";\n	    }\n	  \\
n	  print \"\\n![dynamic.pl] --- EXIT: $ec ($EXIT_\
SUCCESS:success, $EXIT_FAILURE:failure)-- Verbose \
mode -- unset VERBOSE_4_DYNAMIC to turn verbose mo\
de off\\n\";\n	}\n      chdir ($dir);\n      exit \
($ec);\n    }\n	      \nsub my_system \n  {\n    m\
y ($com)=@_;\n    $LAST_COM=$com;\n    \n    if ($\
VERBOSE){print \"\\n![dynamic.pl] --- SysCall --- \
$com\\n\";}\n    \n    system ($com);\n  }\nmy $pr\
ogram=\"T-COFFEE (Version_13.45.15.8d1d86e)\";\n",\
"use Env;\nuse FileHandle;\nuse Cwd;\nuse File::Pa\
th;\nuse Sys::Hostname;\nmy $f = new FileHandle;\n\
\nopen ($f, $ARGV[1]);\n$atom=$ARGV[0];\n\n$atom=~\
s/PRIME/\\'/;\nwhile (<$f>)\n  {\n    my $l=$_;\n\\
n    $l=~s/$atom/CA /;\n    \n    \n    $l=~s/  G \
/GLY /g;\n    $l=~s/  C /CYS /g;\n    $l=~s/  T /T\
HR /g;\n    $l=~s/  A /ALA /g;\n    $l=~s/  U /THR\
 /g;\n    \n    $l=~s/ DG /GLY /g;\n    $l=~s/ DC \
/CYS /g;\n    $l=~s/ DT /THR /g;\n    $l=~s/ DA /A\
LA /g;\n    $l=~s/ DU /THR /g;\n    \n    print $l\
;\n  }\n\n\n\n","*TC_METHOD_FORMAT_01\n***********\
*******generic_method.tc_method*************\n*\n*\
       Incorporating new methods in T-Coffee\n*   \
    Cedric Notredame 26/08/08\n*\n****************\
***************************************\n*This fil\
e is a method file\n*Copy it and adapt it to your \
need so that the method \n*you want to use can be \
incorporated within T-Coffee\n********************\
***********************************\n*            \
      USAGE                              *\n******\
*************************************************\\
n*This file is passed to t_coffee via -in:\n*\n*	t\
_coffee -in Mgeneric_method.method\n*\n*	The metho\
d is passed to the shell using the following\n*cal\
l:\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM\
2><OUT_FLAG><outname><PARAM>\n*\n*Conventions:\n*<\
FLAG_NAME> 	<TYPE>		<VALUE>\n*<VALUE>:	no_name 	<=\
> Replaced with a space\n*<VALUE>:	&nbsp	<=> Repla\
ced with a space\n*\n*****************************\
**************************\n*                  ALN\
_MODE                           *\n***************\
****************************************\n*pairwis\
e   ->all Vs all (no self )[(n2-n)/2aln]\n*m_pairw\
ise ->all Vs all (no self)[n^2-n]^2\n*s_pairwise -\
>all Vs all (self): [n^2-n]/2 + n\n*multiple   ->A\
ll the sequences in one go\n*\nALN_MODE		pairwise\\
n*\n**********************************************\
*********\n*                  OUT_MODE            \
               *\n********************************\
***********************\n* mode for the output:\n*\
External methods: \n* aln -> alignmnent File (Fast\
a or ClustalW Format)\n* lib-> Lib file (TC_LIB_FO\
RMAT_01)\n*Internal Methods:\n* fL -> Internal Fun\
ction returning a List (Librairie)\n* fA -> Intern\
al Function returning an Alignmnent\n*\nOUT_MODE		\
aln\n*********************************************\
**********\n*                  SEQ_TYPE           \
                *\n*******************************\
************************\n*G: Genomic, S: Sequence\
, P: PDB, R: Profile\n*Examples:\n*SEQTYPE	S	seque\
nces against sequences (default)\n*SEQTYPE	S_P	seq\
uence against structure\n*SEQTYPE	P_P	structure ag\
ainst structure\n*SEQTYPE	PS	mix of sequences and \
structure	\n*\nSEQ_TYPE	S\n*\n\n******************\
*************************************\n*          \
      COMMAND LINE                         *\n*EXE\
CUTABLE PARAM1 IN_FLAG OUT_FLAG PARAM             \
*\n***********************************************\
********\n****************************************\
***************\n*                  EXECUTABLE    \
                     *\n**************************\
*****************************\n*name of the execut\
able\n*passed to the shell: executable\n*	\nEXECUT\
ABLE	tc_generic_method.pl\n*\n********************\
***********************************\n*            \
      IN_FLAG                             *\n*****\
**************************************************\
\n*IN_FLAG\n*flag indicating the name of the in co\
ming sequences\n*IN_FLAG S no_name ->no flag\n*IN_\
FLAG S &bnsp-in&bnsp -> \" -in \"\n*\nIN_FLAG		-in\
file=\n*\n****************************************\
***************\n*                  OUT_FLAG      \
                     *\n**************************\
*****************************\n*OUT_FLAG\n*flag in\
dicating the name of the out-coming data\n*same co\
nventions as IN_FLAG\n*OUT_FLAG	S no_name ->no fla\
g\n*if you want to redirect, pass the parameters v\
ia PARAM1\n*set OUT_FLAG to >\n*\nOUT_FLAG		-outfi\
le=\n*\n******************************************\
*************\n*                  PARAM_1         \
                     *\n**************************\
*****************************\n*<EXECUTABLE><PARAM\
1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><P\
ARAM>\n*Parameters sent to the EXECUTABLE and spec\
ified *before* IN_FLAG \n*If there is more than 1 \
PARAM line, the lines are\n*concatenated\n*Command\
_line: @EP@PARAM@-gapopen%e10%s-gapext%e20\n*	%s w\
hite space\n*	%e equal sign\n*\n*PARAM1	\n*\n*\n*\\
n*************************************************\
******\n*                  PARAM_2                \
              *\n*********************************\
**********************\n*<EXECUTABLE><PARAM1><IN_F\
LAG><seq_file><PARAM2><OUT_FLAG><outname><PARAM>\n\
*Parameters sent to the EXECUTABLE and specified \\
n*after* IN_FLAG and *before* OUT_FLAG\n*If there \
is more than 1 PARAM line, the lines are\n*concate\
nated\n*\n*PARAM1	\n*\n*\n************************\
*******************************\n*                \
  PARAM                              *\n**********\
*********************************************\n*<E\
XECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_\
FLAG><outname><PARAM>\n*Parameters sent to the EXE\
CUTABLE and specified *after* OUT_FLAG\n*If there \
is more than 1 PARAM line, the lines are\n*concate\
nated\n*\nPARAM	-mode=seq_msa -method=clustalw2\nP\
ARAM   -OUTORDER=INPUT -NEWTREE=core -align -gapop\
en=-15\n*\n***************************************\
****************\n*                  END          \
                      *\n*************************\
******************************\n","*TC_METHOD_FORM\
AT_01\n***************clustalw_method.tc_method***\
******\nEXECUTABLE	clustalw\nALN_MODE		pairwise\nI\
N_FLAG		-INFILE=\nOUT_FLAG		-OUTFILE=\nOUT_MODE		a\
ln\nPARAM		-gapopen=-10\nSEQ_TYPE		S\n************\
*************************************\n","$Version\
Tag =                                             \
                                                  \
                                    2.43;\nuse Env\
;\nuse FileHandle;\nuse Cwd;\nuse File::Path;\nuse\
 Sys::Hostname;\n\nour $PIDCHILD;\nour $ERROR_DONE\
;\nour @TMPFILE_LIST;\nour $EXIT_FAILURE=1;\nour $\
EXIT_SUCCESS=0;\n\nour $REFDIR=getcwd;\nour $EXIT_\
SUCCESS=0;\nour $EXIT_FAILURE=1;\n\nour $PROGRAM=\\
"extract_from_pdb\";\nour $CL=$PROGRAM;\n\nour $CL\
EAN_EXIT_STARTED;\nour $debug_lock=$ENV{\"DEBUG_LO\
CK\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\\
nif (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR\
=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$EN\
V{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (\
isshellpid(getppid())){lock4tc(getppid(), \"LLOCK\\
", \"LSET\", \"$$\\n\");}\n      \nour $SILENT=\" \
>/dev/null 2>/dev/null\";\nour $INTERNET=-1;\n\n\n\
\n\n\n\n\nour $BLAST_MAX_NRUNS=2;\nour $EXIT_SUCCE\
SS=0;\nour $EXIT_FAILURE=1;\nour $CONFIGURATION=-1\
;\nour $REF_EMAIL=\"\";\nour $PROGRAM=\"extract_fr\
om_pdb\";\n\n\nmy %onelett_prot=&fill_onelett_prot\
();\nmy %threelett_prot=&fill_threelett_prot();\nm\
y %onelett_RNA=&fill_onelett_RNA();\nmy %threelett\
_RNA=&fill_threelett_RNA();\nmy %onelett_DNA=&fill\
_onelett_DNA();\nmy %threelett_DNA=&fill_threelett\
_DNA();\n\n\n\n\n\nmy %onelett = (\n'P' => \\%onel\
ett_prot,\n'D' => \\%onelett_DNA,\n'R' => \\%onele\
tt_RNA\n);\n\n\nmy %threelett = (\n'P' => \\%three\
lett_prot,\n'D' => \\%threelett_DNA,\n'R' => \\%th\
reelett_RNA\n);\n\n\n\n\n\n\n\nif($ARGV[0]=~/help/\
 ||$ARGV[0]=~/man/ || $ARGV[0]=~/HELP/ || $ARGV[0]\
=~/Man/ || $ARGV[0] eq \"-h\"  || $ARGV[0] eq \"-H\
\"  )\n{die \"SYNTAX: extract_from_pdb Version $Ve\
rsionTag	\n	Minimum:             [extract_from_pdb\
 file] \n			   OR \n			     [... | extract_from_pd\
b]\n 	Flags (Default setting on the first line)\n	\
   -version...................[Returns the Version\
 Number]\n           -force.....................[F\
orces the file to be treated like a PDB file]\n   \
                                   [Regenerates th\
e header and SEQRES fields]\n           -force_nam\
e................[Forces the file to be named afte\
r name]]\n           -infile.....file...........[F\
lag can be omited]\n			              [File must be\
 pdb or fro pgm]\n                                \
      [File can also be compressed Z or gz]\n     \
                                 [In the case of a\
 compressed file, you can omit the gz|Z extension]\
\n           -netfile...................[File will\
 be fetch from the net using wget]\n              \
                        [wget or curl must be inst\
alled]\n                                      [ftp\
://ftp.gnu.org/pub/gnu/wget/]\n                   \
                   [http://curl.haxx.se/]\n       \
                               [Must also be used \
to retrieve the file from a local pdb copy (cf net\
address)]\n           -netaddress................[\
Address used for the retrieving the netfile]\n    \
                                  [http://www.rcsb\
.org/pdb/cgi/export.cgi/%%.pdb.gz?format=PDB&pdbId\
=%%&compression=gz]\n                             \
         [http://www.expasy.ch/cgi-bin/get-pdb-ent\
ry.pl?%%]\n                                      [\
local -> will get the file from pdb_dir (see pdb_d\
ir)]\n           -netcompression............[Exten\
sion if the netfile comes compressed]\n           \
                           [gz]\n           -pdb_d\
ir...................[address of the repertory whe\
re the pdb is installed]\n                        \
              [Supports standard ftp style install\
ation OR every stru in DIR]\n                     \
                 [Give the ..../pdb/structure/ dir\
]\n                                      [If value\
 omitted, the pg gets it from the env variable PDB\
_DIR]\n           -netcompression_pg.........[gunz\
ip]\n           -is_pdb_name..........name.[Return\
s 1 if the name is a PDB ID, 0 otherwise]\n       \
    -model_type...........name.[Returns the model \
type if valid PDB name]\n           -is_released_p\
db_name name.[Returns 1 if the name corresponds to\
 a released PDB file]\n           -get_pdb_chains.\
....name...[Returns the list of chains correspondi\
ng to the entry]\n           -get_pdb_id.........n\
ame...[Returns the PDB id within the provided pdb \
file]\n           -get_fugue_name.....name...[Turn\
s a name into a name valid for fugue]\n           \
                           [Uses the netaddress to\
 do so]\n	   -chain......FIRST..........[Extract t\
he first chain only]\n		       A B C..........[Ext\
ract Several chains if needed]\n		       ALL......\
......[Extract all the chains]	\n           -ligan\
d.....ALL............[Extract the ligands in the c\
hain (HETATM)]\n                       <name1>,<na\
me2>[Extract All the named lignds]\n	   -ligand_on\
ly...............[Extract only the ligands]\n     \
      -ligand_list...............[Extract the list\
 of ligands]\n	   -coor.......<start>..<end>.[Coor\
dinates of the fragment to extract]\n			          \
    [Omit end to include the Cter]\n           -nu\
m........absolute.......[absolute: relative to the\
 seq] \n                       file...........[fil\
e: relative to file]\n           -num_out....new..\
..........[new: start 1->L]\n                     \
  old............[old: keep the file coordinates]\\
n           -delete.....<start>..<end>.[Delete fro\
m residue start to residue end]\n	   -atom.......C\
A.............[Atoms to include, ALL for all of th\
em]\n		       CA O N.........[Indicate several ato\
ms if needed]\n	   -code.......3..............[Use\
 the 1 letter code or the 3 letters code]\n	   -mo\
de.......raw............[Output original pdb file]\
\n                       pdb............[Output so\
mething that looks like pdb]\n		       fasta......\
....[Output the sequences in fasta format]\n		    \
   simple.........[Output a format easy to parse i\
n C ]\n            -seq_field..ATOM...........[Fie\
ld used to extract the sequence]\n		       SEQRES.\
........[Use the complete sequence]\n	   -seq.....\
..................[Equivalent to  -mode fasta]\n	 \
  -model......1..............[Chosen Model in an N\
MR file]\n           -nodiagnostic..............[S\
witches Error Messages off]\n           -debug....\
.................[Sets the DEBUG ON]\n           -\
no_remote_pdb_dir.........[Do not look for a remot\
e file]\n           -cache_pdb.................[Ca\
che Value, default is $HOME/.t_coffee/cache, other\
 values: NO<=> No cache]\n\n      Environement Var\
iables\n           These variables can be set from\
 the environement\n           Command line values \
with the corresponding flag superseed evironement \
value\n           NO_REMOTE_PDB_DIR..........[Prev\
ents the program from searching remote file: faste\
r]\n           PDB_DIR....................[Indicat\
es where PDB file must be fetched (localy)]\n\n	 P\
ROBLEMS: please contact cedric.notredame\\@europe.\
com\\n\";\n	 exit ($EXIT_SUCCESS);\n}\n\n$np=0;\n$\
n_para=$#ARGV;\n$model=1;\n$pdb_dir=$ENV{'PDB_DIR'\
};if ($pdb_dir){$pdb_dir.=\"/\";}\n$debug=$ENV{'DE\
BUG_EXTRACT_FROM_PDB'};\n\n$no_remote_pdb_dir=$ENV\
{NO_REMOTE_PDB_DIR};\n$HOME=$ENV{'HOME'};\nif ( $E\
NV{CACHE_4_TCOFFEE})\n{$cache=$ENV{CACHE_4_TCOFFEE\
};}\nelse\n{\n    $cache=\"$HOME/.t_coffee/cache/\\
";\n}\n\n   \n$netaddress=\"http://www.rcsb.org/pd\
b/files/%%.pdb.gz\";\n$netcompression_pg=\"gunzip\\
";\n$netcompression=\"gz\";\n\nforeach ($np=0; $np\
<=$n_para; $np++)\n  {        \n    $value=$ARGV[$\
np];\n   \n    if  ($np==0 && !($value=~/^-.*/))\n\
      { \n	$pdb_file= $ARGV[$np];\n      }\n    el\
sif ( !($value=~/^-.*/))\n      {\n	print \"@ARGV\\
";\n	die;\n      } \n    \n    elsif ($value eq \"\
-nodiagnostic\"){$nodiagnostic=1;}\n    elsif ($va\
lue eq \"-force\")\n      {\n	$force_pdb=1;\n     \
 }\n    elsif ($value eq \"-force_name\")\n      {\
\n	$force_name=$ARGV[++$np];\n	$force_pdb=1;\n    \
  }\n    \n    elsif ($value eq \"-is_pdb_name\")\\
n      {\n	$pdb_file= $ARGV[++$np];	\n	$is_pdb_nam\
e=1;	\n      } \n    elsif ($value eq \"-is_releas\
ed_pdb_name\")\n      {\n	$pdb_file= $ARGV[++$np];\
\n	\n	if (!$pdb_file){print \"0\";exit (EXIT_SUCCE\
SS);}\n	$is_released_pdb_name=1;\n      }\n    els\
if ($value eq \"-model_type\")\n      {\n	$pdb_fil\
e= $ARGV[++$np];	\n	$model_type=1;\n      }\n    e\
lsif ($value eq \"-debug\")\n{\n	$debug=1;\n}\n   \
 elsif ($value eq \"-get_pdb_chains\")\n{\n	$pdb_f\
ile= $ARGV[++$np];\n	$get_pdb_chains=1;\n}\n    el\
sif ($value eq \"-get_pdb_ligands\")\n{\n	$get_pdb\
_ligands=1;\n}\n    \n    elsif ($value eq \"-get_\
pdb_id\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_pdb\
_id=1;\n	\n}\n    \n    elsif ( $value eq \"-get_f\
ugue_name\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_\
fugue_name=1;\n}\n    elsif ( $value eq \"-infile\\
")\n{\n       $pdb_file= $ARGV[++$np];\n} \n    el\
sif ($value eq \"-netfile\")\n{\n	$netfile=1;\n	if\
 ( !($ARGV[$np+1]=~/^-.*/)){$pdb_file= $ARGV[++$np\
];}\n}\n    elsif (  $value eq \"-num\")\n{\n     \
  $numbering= $ARGV[++$np];\n}\n    elsif (  $valu\
e eq \"-num_out\")\n{\n       $numbering_out= $ARG\
V[++$np];\n}\n    elsif ( $value eq \"-netaddress\\
")\n{\n	$netadress=$ARGV[++$np];\n}\n     \n    el\
sif ( $value eq \"-netcompression\")\n{\n	 $netcom\
pression=$ARGV[++$np];\n}\n    elsif ( $value eq \\
"-pdb_dir\")\n{\n	 if ( !($ARGV[$np+1]=~/^-.*/)){$\
pdb_dir= \"$ARGV[++$np]/\";}\n}\n     elsif ( $val\
ue eq \"-no_remote_pdb_dir\")\n{\n	$no_remote_pdb_\
dir=1;\n	if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_dir= \\
"$ARGV[++$np]/\";}\n}\n    elsif ( $value eq \"-ca\
che\")\n{\n	$cache=$ARGV[++$np];\n}\n    \n    els\
if ($value eq \"-netcompression_pg\")\n{\n	  $netc\
ompression_pg=$ARGV[++$np];\n}\n     elsif ($value\
 eq \"-mode\")\n{\n       $MODE=$ARGV[++$np];\n}\n\
\n    elsif ( $value eq \"-model\")\n{\n       $mo\
del= $ARGV[++$np];\n}\n    elsif ($value eq \"-seq\
_field\" )\n{\n       $seq_field= $ARGV[++$np];\n}\
   \n    elsif ($value eq \"-coor\" )\n{\n       $\
start= $ARGV[++$np];\n  \n       if (($ARGV[$np+1]\
 eq \"\") ||($ARGV[$np+1]=~/^-.*/)){$end=\"*\";} \\
n       else {$end=   $ARGV[++$np];}     \n       \
$coor_set=1;\n}\n    elsif ($value eq \"-delete\" \
)\n{\n       $delete_start= $ARGV[++$np];\n       \
$delete_end= $ARGV[++$np];\n       $delete_set=1;\\
n}\n    elsif  ($value eq \"-code\")\n{\n       $c\
ode= $ARGV[++$np];\n}\n    elsif  ($value eq \"-no\
_hetatm\")\n{\n       $no_hetatm=1;\n}\n    elsif \
($value eq \"-chain\")\n{\n       while (!($ARGV[$\
np+1] eq \"\") &&!($ARGV[$np+1]=~/^-.*/))\n{\n	   \
   ++$np;\n	      @c_chain=(@chain,  $ARGV[$np]);\\
n	      $hc_chain{$ARGV[$np]}=$#c_chain+1;\n}     \
      \n}\n    elsif ($value eq \"-atom\")\n{\n\n \
      while (!($ARGV[$np+1] eq \"\") && !($ARGV[$n\
p+1]=~/^-.*/))\n{\n	      ++$np;\n	      $atom[$n_\
atom++]=  $ARGV[$np];\n	      $atom_list{$ARGV[$np\
]}=1;	      \n} \n       \n}\n    elsif ( $value e\
q \"-unfold\")\n{\n	$unfold=1;\n}\n    elsif ($val\
ue eq \"-seq\" ||$value eq \"-fasta\" )\n{\n      \
 $MODE=\"fasta\";\n}\n    elsif ( $value eq \"-ver\
sion\")\n{\n	print STDERR  \"\\nextract_from_pdb: \
Version $VersionTag\\n\";\n	&myexit ($EXIT_SUCCESS\
);\n}\n    elsif ( $value eq \"-ligand\")\n{\n	whi\
le (!($ARGV[$np+1] eq \"\") && !($ARGV[$np+1]=~/^-\
.*/))\n{\n	    ++$np;\n	    $ligand=1;\n	    $liga\
nd_list{$ARGV[$np]}=1;	      \n} \n	$hc_chain{'LIG\
AND'}=1;\n}\n    elsif ( $value eq \"-ligand_only\\
")\n{\n	$ligand_only=1;\n}\n}\nif ( $debug)\n{\n  \
  print STDERR \"\\n[DEBUG:extract_from_pdb] NO_RE\
MOTE_PDB_DIR: $no_remote_pdb_dir\\n\";\n    print \
STDERR \"\\n[DEBUG:extract_from_pdb] PDB_DIR: $pdb\
_dir\\n\";\n}\n\n\nif ( $is_pdb_name)\n  {\n    if\
 (&remote_is_pdb_name($pdb_file))\n      {\n	print\
 \"1\";\n      }\n    else\n      {\n	print \"0\";\
\n      }\n    exit ($EXIT_SUCCESS);\n  }\n\nif ( \
$is_released_pdb_name)\n  {\n    \n    if (&is_rel\
eased($pdb_file))\n      {\n	print \"1\";\n      }\
\n    else\n      {\n	print \"0\";\n      }\n    e\
xit ($EXIT_SUCCESS);\n  }\nif ($model_type)\n  {\n\
   \n    printf \"%s\", &pdb2model_type($pdb_file)\
;\n    exit ($EXIT_SUCCESS);\n    \n  }\n    \n\ni\
f (!$force_name)\n{\n    $pdb_file=~/([^\\/]*)$/;\\
n    $force_name=$1;\n}\n\n$local_pdb_file=$pdb_fi\
le;\n\nif ( $debug){print STDERR \"\\n[DEBUG: extr\
act_from_pdb] Scan For $local_pdb_file\\n\";}\n\n$\
mem=$no_remote_pdb_dir;\n$no_remote_pdb_dir=1;\n$t\
mp_pdb_file=get_pdb_file ($local_pdb_file);\n\nif \
( !-e $tmp_pdb_file || $tmp_pdb_file eq \"\")\n  {\
\n    $local_pdb_file=$pdb_file;\n    ($local_pdb_\
file, $suffix_chain)=&pdb_name2name_and_chain($loc\
al_pdb_file);\n\n    if ($local_pdb_file)\n      {\
\n	if ( $debug){print STDERR \"\\nSplit $pdb_file \
into $local_pdb_file and $suffix_chain \\n\";}\n	$\
tmp_pdb_file=get_pdb_file ($local_pdb_file);\n	if \
( $tmp_pdb_file ne \"\")\n	  {\n	    @c_chain=();\\
n	    @c_chain=($suffix_chain);\n	    %hc_chain=()\
;\n	    $hc_chain{$suffix_chain}=1;\n	  }\n      }\
\n  }\n\n$no_remote_pdb_dir=$mem;\nif ($no_remote_\
pdb_dir==0)\n  {\n    \n    if ( !-e $tmp_pdb_file\
 || $tmp_pdb_file eq \"\")\n      {\n	\n	$local_pd\
b_file=$pdb_file;\n	($local_pdb_file, $suffix_chai\
n)=&pdb_name2name_and_chain($local_pdb_file);\n	if\
 ($local_pdb_file)\n	  {\n	    \n	    if ( $debug)\
{print STDERR \"\\nSplit $pdb_file into $local_pdb\
_file and $suffix_chain \\n\";}\n	    \n	    $tmp_\
pdb_file=get_pdb_file ($local_pdb_file);    \n	   \
 \n	    if ( $tmp_pdb_file ne \"\")\n	      {\n		@\
c_chain=();\n		@c_chain=($suffix_chain);\n		%hc_ch\
ain=();\n		$hc_chain{$suffix_chain}=1;\n	      }\n\
	  }\n      }\n  }\n\nif ( $debug){print STDERR \"\
\\n$pdb_file copied into ##$tmp_pdb_file##\\n\";}\\
n\nif ( !-e $tmp_pdb_file || $tmp_pdb_file eq \"\"\
)\n{\n	if ($is_pdb_name)\n{\n	    print \"0\\n\"; \
exit ($EXIT_SUCCESS);\n}\n	else\n{\n  \n	    print\
 \"\\nEXTRACT_FROM_PDB: NO RESULT for $pdb_file\\n\
\";\n	    &myexit ($EXIT_SUCCESS);	\n}\n}\n\n\n\n\\
n%molecule_type=&pdbfile2chaintype($tmp_pdb_file);\
\nif ( $debug){print STDERR \"\\n\\tSequence Type \
determined\\n\";}\n\n$pdb_id=&get_pdb_id ($tmp_pdb\
_file);\nif ( $debug){print STDERR \"\\n\\tID: $pd\
b_id (for $tmp_pdb_file)\\n\";}\n\nif ( $pdb_id eq\
 \"\"){$pdb_id=$force_name;}\n\n@f_chain=&get_chai\
n_list ($tmp_pdb_file);\nif ( $debug){print STDERR\
 \"\\n\\tChain_list:@f_chain\\n\";}\n\nif ( $get_p\
db_chains)\n{\n    print \"@f_chain\\n\";\n    &my\
exit ($EXIT_SUCCESS);\n}\nif ( $get_pdb_ligands)\n\
{\n    %complete_ligand_list=&get_ligand_list ($tm\
p_pdb_file);\n    print $complete_ligand_list{\"re\
sult\"};\n    &myexit ($EXIT_SUCCESS);\n}\n\nelsif\
 ( $get_pdb_id ||$get_fugue_name )\n{\n    if    (\
@c_chain && $c_chain[0] eq \"FIRST\"){$pdb_id=$pdb\
_id.$f_chain[0];}\n    elsif (@c_chain && $c_chain\
[0] ne \" \"){$pdb_id=$pdb_id.$c_chain[0];}\n    \\
n    print \"$pdb_id\\n\";\n    &myexit ($EXIT_SUC\
CESS);\n    \n}\nelsif ( $is_pdb_name)\n{\n    pri\
ntf \"1\\n\";\n    &myexit ($EXIT_SUCCESS);\n}\n\n\
\n\n$structure_file=vtmpnam();\n\nif ( $debug){pri\
nt STDERR \"\\n\\tCheck_point #1: $tmp_pdb_file  $\
structure_file\\n\";}\n\n$INFILE=vfopen (\"$tmp_pd\
b_file\", \"r\"); \nmy $TMP=vfopen (\"$structure_f\
ile\", \"w\");\n\n$print_model=1;\n$in_model=0;\n\\
nif ( $debug){print STDERR \"\\n\\tCheck_point #2\\
\n\";}\nwhile ( <$INFILE>)\n{\n  my $first_model=0\
;\n  $line=$_;\n\n  if ( !$first_model && ($line =\
~/^MODEL\\s*(\\d*)/))\n    {\n      $first_model=$\
1;\n      if ($model==1){$model=$first_model;}\n  \
  }\n  \n  if (($line =~/^MODEL\\s*(\\d*)/))\n    \
{\n      if ($1==$model)\n	{\n	  $in_model=1;\n	  \
$print_model=1;\n	  $is_nmr=1;\n	}\n      elsif ( \
$in_model==0)\n	{\n	  $print_model=0;\n	}\n      e\
lsif ( $in_model==1)\n	{\n	  last;\n	}\n    }\n  i\
f ($print_model){print $TMP $line;}  \n}\nclose ($\
TMP);\nclose ($INFILE);\n\nif ( $debug){print STDE\
RR \"\\n\\tCheck_point #3\\n\";}	\n\n  if ($number\
ing eq \"\"){$numbering=\"absolute\";}\n  if ($num\
bering_out eq \"\"){$numbering_out=\"new\";}\n\n  \
if ( $delete_set && $coor_set) {die \"-delete and \
-coor are mutually exclusive, sorry\\n\";}\n  if (\
 $n_atom==0){$atom_list[$n_atom++]=\"ALL\";$atom_l\
ist{$atom_list[0]}=1;}\n  if ( $seq_field eq \"\")\
{$seq_field=\"ATOM\";}\n  \n  if ( $MODE eq \"\"){\
$MODE=\"pdb\";}\n  elsif ( $MODE eq \"simple\" && \
$code==0){$code=1;}\n\n  if ( $code==0){$code=3;}\\
n\n\nif ($f_chain[0] eq \" \"){$hc_chain{' '}=1;$c\
_chain[0]=\" \";}\nelsif (!@c_chain){$hc_chain{FIR\
ST}=1;$c_chain[0]=\"FIRST\";}#make sure the first \
chain is taken by default\n\nif    ($hc_chain{ALL}\
) \n{\n      @c_chain=@f_chain;\n      foreach $e \
(@c_chain){$hc_chain{$e}=1;}\n}\nelsif($hc_chain{F\
IRST})\n{\n	@c_chain=($f_chain[0]);\n	$hc_chain{$f\
_chain[0]}=1;\n}\n\n\n$MAIN_HOM_CODE=&get_main_hom\
_code ($structure_file);\n$INFILE=vfopen ($structu\
re_file, \"r\");\n\n\nif ( $MODE eq \"raw_pdb\" ||\
 $MODE eq \"raw\")\n{\n    while (<$INFILE>)\n{	pr\
int \"$_\";}\n    close ( $INFILE);\n    &myexit($\
EXIT_SUCCESS);\n}    \nif ( $MODE eq \"raw4fugue\"\
 )\n{\n    while (<$INFILE>)\n{	\n	$l=$_;\n	if ($l\
=~/^SEQRES/)\n{\n	    \n	    $c= substr($l,11,1);\\
n	    if ($hc_chain {$c}){print \"$l\";}\n}\n	elsi\
f ( $l=~/^ATOM/)\n{\n	    $c=substr($l,21,1);\n	  \
  if ($hc_chain {$c}){print \"$l\";}\n}\n}\n    cl\
ose ( $INFILE);\n    &myexit($EXIT_SUCCESS);\n}   \
 \n\nif ( $MODE eq \"pdb\")\n{\n\n    $read_header\
=0;\n    while (<$INFILE>) \n{\n	    $line=$_;\n	 \
   if    ($line =~ /^HEADER/){print \"$line\";$rea\
d_header=1;}\n}\n    close ($INFILE);\n\n    if (!\
$read_header)\n{\n	print \"HEADER    UNKNOWN      \
                           00-JAN-00   $force_name\
\\n\";\n}\n\n    $INFILE=vfopen ($structure_file, \
\"r\");\n    \n    print \"COMPND   1 CHAIN:\";\n \
   $last=pop(@c_chain);\n    foreach $c ( @c_chain\
){ print \" $c,\";}\n    if ( $last eq \" \"){prin\
t \" NULL;\\n\";}\n    else \n{\n      print \" $l\
ast;\\n\";\n}\n    @c_chain=(@c_chain, $last);\n  \
  \n    print \"REMARK Output of the program extra\
ct_from_pdb (Version $VersionTag)\\n\";\n    print\
 \"REMARK Legal PDB format not Guaranteed\\n\";\n \
   print \"REMARK This format is not meant to be u\
sed in place of the PDB format\\n\";\n    print \"\
REMARK The header refers to the original entry\\n\\
";\n    print \"REMARK The sequence from the origi\
nal file has been taken in the field: $seq_field\\\
n\";\n    print \"REMARK extract_from_pdb, 2001, 2\
002, 2003, 2004, 2005 2006 (c) CNRS and Cedric Not\
redame\\n\";   \n    if ( $coor_set)\n{\n       pr\
int \"REMARK Partial chain: Start $start End $end\\
\n\";\n}\n    if ( $is_nmr)\n{\n       print \"REM\
ARK NMR structure: MODEL $model\\n\";\n}\n   \n   \
 if ( $n_atom!=0)\n{\n       print \"REMARK Contai\
ns Coordinates of: \";\n       foreach $a (@atom){\
print \"$a \";}\n       print \"\\n\";\n}  \n}\n\n\
\n\n\nmy $residue_index = -999;\nmy $old_c = \"Tem\
poraryChain\";\n\nwhile (<$INFILE>) \n{\n	$line=$_\
;\n\n\n	if ($line =~ /^SEQRES/)\n{\n\n		@field=/(\\
\S*)\\s*/g;\n\n		$c= substr($_,11,1);\n\n		\n		$l=\
$#field;\n		for ($a=4; $a<$#field ;)\n{\n			if (!$\
onelett{$molecule_type{$c}}->{$field[$a]})\n{\n			\
	splice @field, $a, 1;\n}\n			else \n{\n				$a++;\\
n}\n}\n	\n		if ( $c ne $in_chain)\n{\n			$pdb_chai\
n_list[$n_pdb_chains]=$c;\n			$pdb_chain_len [$n_p\
db_chains]=$len;\n			$in_chain=$c;\n			$n_pdb_chai\
ns++;\n}\n	\n		for ( $a=4; $a<$#field;$a++)\n{\n		\
	$complete_seq{$c}[$complete_seq_len{$c}++]=$field\
[$a];\n}\n}\n    elsif ( $line=~/^ATOM/ || ($line=\
~/^HETATM/ && &is_aa(substr($line,17,3),substr($li\
ne,21,1)) && !$no_hetatm))\n{\n\n	 \n    $RAW_AT_I\
D=$AT_ID=substr($line,12,4);\n	$RES_ID=&is_aa(subs\
tr($line,17,3),substr($line,21,1));\n	$CHAIN=subst\
r($line,21,1);\n\n    $RES_NO=substr($line,22,4);\\
n	$HOM_CODE=substr ($line, 26, 1);\n	$TEMP=substr(\
$line,60,6);\n	\n	$TEMP=~s/\\s//g;\n        $AT_ID\
=~s/\\s//g;\n	$RES_ID=~s/\\s//g;\n        $RES_NO=\
~s/\\s//g;\n		\n	if ( $HOM_CODE ne $MAIN_HOM_CODE)\
{next;}\n	elsif ( $already_read2{$CHAIN}{$RES_ID}{\
$AT_ID}{$RES_NO}){next;}\n	else{$already_read2{$CH\
AIN}{$RES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	\n	if ($co\
or_set && $numbering eq \"file\" && $residue_index\
 ne $RES_NO)\n{\n	    \n	    if ( $RES_NO<=$start)\
{$real_start{$CHAIN}++;}\n	    if ( $RES_NO<=$end)\
{$real_end{$CHAIN}++;}\n}\n	elsif ($numbering eq \\
"absolute\")\n{\n	    $real_start{$CHAIN}=$start;\\
n	    $real_end{$CHAIN}=$end;\n}\n\n        $KEY=\\
"ALL\";\n        if ( $CHAIN ne $in_atom_chain)\n{\
\n	    \n	  $pdb_atom_chain_list[$n_pdb_atom_chain\
s]=$c;\n	  $pdb_atom_chain_len [$n_pdb_atom_chains\
]=$len;\n	  $in_atom_chain=$c;\n	  $n_pdb_atom_cha\
ins++;\n}\n	\n	if ( $residue_index ne $RES_NO)\n{\\
n	     $residue_index = $RES_NO;\n	     $atom_seq{\
$CHAIN}[$atom_seq_len{$CHAIN}++]=$RES_ID;;\n}\n}\n\
\n}\nclose ($INFILE);\n\n\n\n\n\n\n$INFILE=vfopen \
($structure_file, \"r\");\nforeach $c (@c_chain)\n\
{\n\n	if    ( $seq_field eq \"SEQRES\"){@pdb_seq=@\
{$complete_seq{$c}};}\n	elsif ( $seq_field eq \"AT\
OM\")  {@pdb_seq=@{$atom_seq{$c}};}\n	\n\n	$full_l\
ength=$l=$#pdb_seq+1;\n		\n	if ( $real_end{$c}==\"\
*\"){$real_end{$c}=$full_length;}\n	if ( $coor_set\
)\n{	   \n\n	   if ( $real_end{$c} < $l){splice @p\
db_seq, $real_end{$c}, $l;}\n	   if ( $real_start{\
$c} < $l){splice @pdb_seq, 0, $real_start{$c}-1;}	\
  	   \n	   $l=$#pdb_seq;\n}\n\n	elsif ( $delete_s\
et)\n{\n	   splice @pdb_seq, $delete_start, $delet\
e_end-$delete_start+1;\n	   $l=$#pdb_seq;\n}\n	\n	\
$new_fasta_name=\"$pdb_id$c\";\n	if ( $coor_set)\n\
{\n	   if ( $n_pdb_chains==0){$new_fasta_name=\"$n\
ew_fasta_name$c\";}\n	   $new_fasta_name= $new_fas\
ta_name.\"\\_$start\\_$end\";\n}\n	   \n	if ( $MOD\
E eq \"pdb\")\n{\n	   $nl=1;\n	   $n=0;\n	   \n	  \
 foreach $res ( @pdb_seq)\n		{\n		if ( !$n)\n		{\n\
		\n		 printf \"SEQRES %3d %1s %4d  \", $nl,$c, $l\
;\n		 $nl++;\n	}\n	     $res=~s/\\s//g;\n	     \n	\
     if ($code==1){ printf \"%3s \",$onelett{$mole\
cule_type{$c}}->{$res};}\n	     elsif  ($code==3){\
 printf \"%3s \",$res};\n	     \n	     $n++;		  \n\
	     if ( $n==13){$n=0;print \"\\n\";}\n}\n	  if \
( $n!=0){print \"\\n\"; $n=0;}\n}\n	elsif ( $MODE \
eq \"simple\")\n{\n	  print \"# SIMPLE_PDB_FORMAT\\
\n\";\n	  if ( $new_fasta_name eq \" \"){$new_fast\
a_name=\"dummy_name\";}\n	  print \">$new_fasta_na\
me\\n\";\n\n	  foreach $res ( @pdb_seq)\n{\n	     \
 print \"$onelett{$molecule_type{$c}}->{$res}\";\n\
}\n	  print \"\\n\";\n}\n	elsif ( $MODE eq \"fasta\
\")\n{\n	  $n=0;\n	  print \">$new_fasta_name\\n\"\
;\n	  \n	  foreach $res ( @pdb_seq)\n{\n\n	      p\
rint \"$onelett{$molecule_type{$c}}->{$res}\";\n  \
            $n++;\n	      if ( $n==60){print \"\\n\
\"; $n=0;}\n}\n	  print \"\\n\"; \n}\n}\n\nif ( $M\
ODE eq \"fasta\")\n{\n     &myexit($EXIT_SUCCESS);\
\n  \n}\n\n  \n  $charcount=0;\n  $inchain=\"BEGIN\
\";\n  $n=0;\n  while (<$INFILE>) \n{\n    $line=$\
_;\n     \n    if ($line =~/^ATOM/  ||  ($line=~/^\
HETATM/))\n{\n	$line_header=\"UNKNWN\";\n	$RES_ID=\
substr($line,17,3);\n	$chain = substr($line,21,1);\
\n\n	if ($line =~/^ATOM/)\n{\n	    $line_header=\"\
ATOM\";\n	    $RES_ID=(&is_aa($RES_ID,$chain))?&is\
_aa($RES_ID,$chain):$RES_ID;\n}\n	elsif ($line=~/^\
HETATM/ && ($ligand_list {$RES_ID} ||$ligand_list \
{'ALL'} || !&is_aa($RES_ID,$chain)))\n{\n	    $lin\
e_header=\"HETATM\";\n}\n	elsif ($line=~/^HETATM/ \
&& (&is_aa($RES_ID,$chain) && !$no_hetatm))\n{\n	 \
   $line_header=\"ATOM\";\n	    $RES_ID=&is_aa($RE\
S_ID,$chain);\n}\n	else\n{\n	    next;\n}\n\n	\n\n\
	$X=substr($line,30,8);     \n	$Y=substr($line,38,\
8);\n	$Z=substr($line,46,8);\n	$TEMP=substr($line,\
60,6);\n	\n	$RAW_AT_ID=$AT_ID=substr($line,12,4);\\
n	$CHAIN=substr($line,21,1);\n	$RES_NO=substr($lin\
e,22,4);\n	$HOM_CODE=substr ($line, 26, 1);\n	\n	$\
X=~s/\\s//g;\n	$Y=~s/\\s//g;\n	$Z=~s/\\s//g;\n	$TE\
MP=~s/\\s//g;\n	\n	$AT_ID=~s/\\s//g;\n	$RES_ID=~s/\
\\s//g;\n	$RES_NO=~s/\\s//g;\n\n	\n	if ( $HOM_CODE\
 ne $MAIN_HOM_CODE){next;}\n	elsif ( $already_read\
{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}){next;}\n	else{\
$already_read{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}=1;\
}\n	\n	$KEY=\"ALL\";\n\n      	if ( $RES_NO ==0){$\
start_at_zero=1;}\n\n	$RES_NO+=$start_at_zero;    \
\n	\n	if ( $current_chain ne $CHAIN)\n{\n	    $cur\
rent_chain=$CHAIN;\n	    $pos=$current_residue=0;\\
n	    $offset=($coor_set)?($real_start{$CHAIN}-1):\
0;\n	    if    ( $seq_field eq \"SEQRES\"){@ref_se\
q=@{$complete_seq{$CHAIN}};}\n	    elsif ( $seq_fi\
eld eq \"ATOM\")  {@ref_seq=@{$atom_seq{$CHAIN}};}\
\n}\n	\n	if ($current_residue != $RES_NO)\n{\n	   \
 $current_residue=$RES_NO;\n	    if    ( $seq_fiel\
d eq \"SEQRES\"){$pos=$current_residue;}\n	    els\
if ( $seq_field eq \"ATOM\"){$pos++;}\n}\n	\n	\n	i\
f ($n_atom==0 || $atom_list{$AT_ID}==1 || $atom_li\
st{$KEY}==1)\n{ 	\n	    \n	    $do_it=(!@c_chain |\
| $hc_chain{$CHAIN} ||$hc_chain{'LIGAND'} );\n	   \
 \n	    $do_it= ($do_it==1) && ($coor_set==0 ||($p\
os>=$real_start{$CHAIN} && $pos<=$real_end{$CHAIN}\
));\n	    $do_it= ($do_it==1) && ($delete_set==0 |\
| $pos<$delete_start ||$pos>$delete_end );\n	    i\
f ($ligand==0 && $line_header eq \"HETATM\" ){$do_\
it=0;}\n	    if ($ligand_only==1 && $line_header e\
q \"ATOM\" ){$do_it=0;}\n	    if ($ligand==1 && $l\
ine_header eq \"HETATM\" && $ligand_list{$RES_ID}=\
=0 && $ligand_list{\"ALL\"}==0){$do_it=0;} \n	    \
\n	    \n	    if ( $do_it)\n{\n		$n++;\n		$out_pos\
=$pos;\n		\n	       if ( $delete_set)\n{\n		  if (\
 $out_pos< $delete_start){;}\n		  else {$offset=$d\
elete_end-$delete_start;}\n}       \n	       \n	  \
     if ( $numbering_out eq \"new\"){$out_pos-=$of\
fset;}\n	       elsif ( $numbering_out eq \"old\")\
{$out_pos=$RES_NO;}\n	       \n       \n	       \n\
	       if ( $code==1){$RES_ID=$onelett{$molecule_\
type{$c}}->{$RES_ID};}\n	    \n	       if ($unfold\
)\n{\n		   $unfolded_x+=5;\n		   $X=$unfolded_x;\n\
		   $Y=0;\n		   $Z=0;\n		   $float=1;\n}\n	      \
 else\n{\n		   $float=3;\n}\n\n	       if ( $MODE \
eq \"pdb\")\n{\n		   printf \"%-6s%5d %-4s %3s %s%\
4d    %8.3f%8.3f%8.3f  1.00 %5.2f\\n\",$line_heade\
r, $n, $RAW_AT_ID,$RES_ID,$CHAIN,$out_pos, $X, $Y,\
 $Z,$TEMP;		  \n}\n	       elsif ( $MODE eq \"simp\
le\")\n{\n		    if ( $RES_ID eq \"\"){$RES_ID=\"X\\
";}\n		  printf \"%-6s %5s %s %2s %4d    %8.3f %8.\
3f %8.3f\\n\",$line_header, $AT_ID, $RES_ID,($CHAI\
N eq\"\" || $CHAIN eq \" \")?\"A\":$CHAIN,$out_pos\
, $X, $Y, $Z,$TEMP;\n}\n\n}\n}\n}\n}\nprint \"\\n\\
";\nclose($INFILE);\n\n\nif ( $error ne \"\") \n{$\
error=$error.\"\\nDiagnostic:    SEQRES and the re\
sidues in ATOM are probably Incompatible\\n\";\n  \
  $error=$error.  \"Recomendation: Rerun with '-fi\
x 1' in order to ignore the SEQRES sequences\\n\";\
\n}\nif (!$nodiagnostic){print STDERR $error;}\n&m\
yexit ( $EXIT_SUCCESS);\n\nsub get_pdb_entry_type_\
file\n  {\n    my $cache_file=\"$cache/pdb_entry_t\
ype.txt\";\n    my $env_file  = $ENV{\"PDB_ENTRY_T\
YPE_FILE\"};\n    my $pdb_file  =\"$ENV{'PDB_DIR'}\
/derived_data/pdb_entry_type.txt\";\n    \n    \n \
   if (-z $cache_file){unlink ($cache_file);}#will\
 get updated\n    if (-z $env_file){$env_file=\"\"\
;}    #cannot update\n    if (-z $pdb_file){$pdb_f\
ile=\"\";}    #cannot update\n    \n    if    (-e \
$env_file){return $env_file;} #env wins: user deci\
des\n    elsif (-e $pdb_file){return $pdb_file;} #\
local database wins: network file may be out of sy\
nc\n    elsif ($no_remote_pdb_dir==1)\n      {\n	i\
f (-e $cache_file){return $cache_file;}\n	else\n	 \
 {add_warning($$,$$,\"PDB_ENTRY_TYPE_FILE must be \
set to the location of <pdb>/derived_data/pdb_entr\
y_type.txt when using NO_REMOTE_PDB_DIR=1\");\n	  \
 return \"\";\n	 }\n      }\n    else #update can \
only take place if the file lives in cache\n      \
{\n	my $new_file;\n	if (!-e $cache_file || (-M $ca\
che_file)>1)\n	  {\n	    $new_file=vtmpnam();\n	  \
  &url2file(\"ftp://ftp.wwpdb.org/pub/pdb/derived_\
data/pdb_entry_type.txt\", $new_file);\n	    if ( \
!-z $new_file){system (\"mv $new_file $cache_file\\
"); unlink ($new_file); $new_file=$cache_file;}\n	\
    else {unlink($new_file);}\n	  }\n	else\n	  {\n\
	    $new_file=$cache_file;\n	  }\n	\n	if (!-e $ca\
che_file && !-e $new_file)\n	  {\n	    add_warning\
($$,$$,\"Could not download ftp://ftp.wwpdb.org/pu\
b/pdb/derived_data/pdb_entry_type.txt\");\n	    re\
turn \"\";\n	  }\n	elsif (-e $cache_file && !-e $n\
ew_file)\n	  {\n	    my $m=(-M $cache_file);\n	   \
 add_warning($$,$$,\"Could not update file ftp://f\
tp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.t\
xt. Older Version [$cache_file]($m Month(s) old) w\
ill be used instead\");\n	    return $cache_file;\\
n	  }\n	else\n	  {\n	    return $new_file;\n	  }\n\
      }\n  }\n\n\n\nsub get_unrealeased_file\n  {\\
n    my $cache_file=\"$cache/unrealeased.xml\";\n \
   my $env_file  = $ENV{\"PDB_UNREALEASED_FILE\"};\
\n    my $pdb_file  =\"$ENV{'PDB_DIR'}/derived_dat\
a/unrealeased.xml\";\n    \n    \n    if ($env_fil\
e eq \"NO\" || $env_file eq \"No\" || $env_file eq\
 \"no\" || $env_file eq \"0\"){return \"NO\";}\n\n\
    if (-z $cache_file){unlink ($cache_file);}#wil\
l get updated\n    if (-z $env_file){unlink($env_f\
ile);}     #will update\n    if (-z $pdb_file){$pd\
b_file=\"\";}          #cannot update\n    \n    i\
f    (-e $env_file){return $env_file;} #env wins: \
user decides\n    elsif (-e $pdb_file){return $pdb\
_file;} #local database wins: network file may be \
out of sync\n    elsif ($no_remote_pdb_dir==1)    \
    \n      {\n	if (-e $cache_file){return $cache_\
file;}\n	elsif ( $env_file && ! -e $env_file)\n	  \
{\n	    &url2file(\"http://www.rcsb.org/pdb/rest/g\
etUnreleased\",$env_file);\n	    if ( -e $env_file\
 && !-z $env_file){return $env_file;}\n	  }\n	else\
\n	  {\n	    add_warning($$,$$,\"UNREALEASED_FILE \
must be set to the location of your unrealeased.xm\
l file as downloaded from http://www.rcsb.org/pdb/\
rest/getUnreleased when using NO_REMOTE_PDB_DIR=1\\
");\n	    return \"\";\n	  }\n      }\n    else #u\
pdate can only take place if the file lives in cac\
he\n      {\n	my $new_file=vtmpnam ();\n	if (!-e $\
cache_file || (-M $cache_file)>1)\n	  {\n	    &url\
2file(\"http://www.rcsb.org/pdb/rest/getUnreleased\
\",$new_file);\n	    if ( !-z $new_file){system (\\
"mv $new_file $cache_file\"); unlink ($new_file); \
$new_file=$cache_file;}\n	    else {unlink($new_fi\
le);}\n	  }\n	else\n	  {\n	    $new_file=$cache_fi\
le;\n	  }\n	\n	if (!-e $cache_file && !-e $new_fil\
e)\n	  {\n	    add_warning($$,$$,\"Could not downl\
oad http://www.rcsb.org/pdb/rest/getUnreleased\");\
\n	    return \"\";\n	  }\n	elsif (-e $cache_file \
&& !-e $new_file)\n	  {\n	    my $m=(-M $cache_fil\
e);\n	    add_warning($$,$$,\"Could not update fil\
e http://www.rcsb.org/pdb/rest/getUnreleased. Olde\
r Version [$cache_file]($m Month(s) ) will be used\
\");\n	    return $cache_file;\n	  }\n	else\n	  {\\
n	    return $new_file;\n	  }\n      }\n  }\n\nsub\
 is_released \n  {\n    my ($r);\n    my $in=@_[0]\
;\n    my $name=&remote_is_pdb_name ($in);\n    my\
 $hold=&remote_is_on_hold($in);\n    \n    $r=($na\
me && !$hold)?1:0;\n    return $r;\n  }\n\nsub rem\
ote_is_pdb_name \n  {\n    my $in=@_[0];\n    my (\
$pdb);\n    my ($value,$value1,$value2);\n    my $\
max=2;\n    \n    \n    \n    my $ref_file=&get_pd\
b_entry_type_file();\n    \n    if ( $in=~/[^\\w\\\
d\\:\\_]/){return 0;}\n    elsif (!-e $ref_file)\n\
      {\n	add_warning ($$,$$,\"Cannot find pdb_ent\
ry_type.txt;  $in is assumed to be valid; add ftp:\
//ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_typ\
e.txt in $cache to automatically check name status\
\");\n	return 1;\n      }\n    else\n      {\n	$pd\
b=substr ($in,0, 4);\n	chomp(($value1=`grep -c $pd\
b $ref_file`));\n	$pdb=lc($pdb);\n	chomp(($value2=\
`grep -c $pdb $ref_file`));\n	$value=($value1 || $\
value2)?1:0;\n	$value=($value>0)?1:0;\n	\n	return \
$value;\n      }\n  }\n\n\n\nsub pdb2model_type\n{\
\n    my $in=@_[0];\n    my ($ref_file, $pdb);\n  \
  my ($value, $ret);\n\n    if ( $in=~/[^\\w\\d\\:\
\\_]/){return 0;}\n    $ref_file=&get_pdb_entry_ty\
pe_file();\n    if (!-e $ref_file)\n      {\n	add_\
warning ($$,$$,\"Cannot find pdb_entry_type.txt;  \
$in is assumed to be diffraction; add ftp://ftp.ww\
pdb.org/pub/pdb/derived_data/pdb_entry_type.txt in\
 $cache to check name status\");\n	return \"diffra\
ction\";\n      }\n    else\n      {\n	$pdb=substr\
 ($in,0, 4);\n	$pdb=lc($pdb);\n	\n	chomp(($value=`\
grep $pdb $ref_file`));\n	\n	$value=~/^\\S+\\s+\\S\
+\\s+(\\S+)/;\n	$ret=$1;\n	if ( $ret eq\"\"){retur\
n \"UNKNOWN\";}\n	\n	return $ret;\n      }\n  }\ns\
ub remote_is_on_hold\n  {\n    my $in=@_[0];\n    \
my ($ref_file, $pdb);\n    my ($value1, $value2,$v\
alue);\n    \n\n\n    \n    $ref_file=&get_unreale\
ased_file();\n    if ($ref_file eq \"NO\"){return \
0;}\n\n\n    if ($no_remote_pdb==1){return 0;}\n  \
  if ( $in=~/[^\\w\\d\\:\\_]/){return 0;}\n    \n \
   $ref_file=&get_unrealeased_file();\n    if (!-e\
 $ref_file)\n      {\n	add_warning ($$,$$,\"Cannot\
 find unrealeased.xml;  $in is assumed to be relea\
sed;\");\n	return 1;\n      }\n    \n    $pdb=subs\
tr ($in,0, 4);\n    chomp(($value1=`grep -c $pdb $\
ref_file`));\n    $pdb=lc($pdb);\n    chomp(($valu\
e2=`grep -c $pdb $ref_file`));\n    $value=($value\
1 || $value2)?1:0;\n    $value=($value>0)?1:0;\n  \
  return $value;\n  }\n\nsub is_pdb_file\n  {\n   \
 my @arg=@_;\n    \n    if ( !-e $arg[0]){return 0\
;}\n    \n    $F=vfopen ($arg[0], \"r\");\n    whi\
le ( <$F>)\n      {\n	if (/^HEADER/)\n	  {\n	    c\
lose $F;\n	    return 1;\n	  }\n	elsif ( /^SEQRES/\
)\n	  {\n	    close $F;\n	    return 1;\n	  }\n	el\
sif ( /^ATOM/)\n	  {\n	    close $F;\n	    return \
1;\n	  }\n      }\n    return 0;\n  }\nsub get_pdb\
_id\n{\n    my $header_file=@_[0];\n    my $id;\n \
   my $F= new FileHandle;\n    \n    \n    $F=vfop\
en (\"$header_file\", \"r\");\n\n    while ( <$F>)\
\n      {\n	if ( /HEADER/)\n	  {\n	    if ($debug)\
{print \"$_\";}\n	    $id=substr($_,62,4 );\n	    \
return $id;\n	  }\n      }\n    close ($F);\n    \\
n    return \"\";\n}\n\nsub get_ligand_list\n{\n  \
  my $pdb_file=@_[0];\n    my $chain;\n    my $lig\
and;\n    my %complete_ligand_list;\n    \n\n    $\
F=vfopen ($pdb_file, \"r\");\n    while ( <$F>)\n{\
\n	if ( /^HETATM/)\n{\n	    $line=$_;\n	    $chain\
=substr($line,21,1);\n	    $ligand=substr($line,17\
,3);\n	    \n	    if (!$complete_ligand_list{$chai\
n}{$ligand})\n{\n		\n		$complete_ligand_list{\"res\
ult\"}.=\"CHAIN $chain LIGAND $ligand\\n\";\n		$co\
mplete_ligand_list{$chain}{$ligand}=1;\n}\n}\n}\n \
   close ($F);\n    return %complete_ligand_list;\\
n}\n\nsub get_chain_list \n{\n    my $header_file;\
\n    my @chain_list;\n    my @list;\n    my $n_ch\
ains;\n    my %chain_hasch;\n    my $pdb_file=@_[0\
];\n    my $c;\n    my %hasch;\n    my $chain;\n  \
\n    \n    $F=vfopen ($pdb_file, \"r\");\n    whi\
le ( <$F>)\n{\n\n\n	if (/SEQRES\\s+\\d+\\s+(\\S+)/\
)\n	  {\n	    $chain = substr($_,11,1);$chain=~s/\\
\s//g;if ( $chain eq \"\"){$chain=\" \";}\n	    if\
 (!$hasch{$chain}){$hasch{$chain}=1;push @chain_li\
st, $chain;}\n	  }\n	if (/^ATOM/ || /^HETATM/)\n	 \
 {\n	    $chain = substr($_,21,1); $chain=~s/\\s//\
g;if ( $chain eq \"\"){$chain=\" \";}\n	    if (!$\
hasch{$chain}){$hasch{$chain}=1;push @chain_list, \
$chain;}\n	  }\n      }\n\n\nclose ($F);\nif (!@ch\
ain_list)\n  {\n    @chain_list=(\"A\");\n  }\n\n\\
nreturn @chain_list;\n}\n\nsub token_is_in_list\n{\
\n\n    my @list=@_;\n    my $a;\n    \n    for ($\
a=1; $a<=$#list; $a++)\n{\n	if ( $list[$a] eq $lis\
t[0]){return $a;}\n}\n}\n\nsub pdb_name2name_and_c\
hain \n{\n    my $pdb_file=@_[0];\n    my $pdb_fil\
e_in;\n    my @array;\n    my $chain;\n    my $c;\\
n\n    $pdb_file_in=$pdb_file;\n\n    $pdb_file=~/\
^(.{4})/;$pdb_id=$1;\n    @array=($pdb_file=~/([\\\
w])/g);\n  \n  \n    $chain=uc ($array[4]);\n    $\
chain=($chain eq \"\")?\"FIRST\":$chain;\n    \n  \
  return ( $pdb_id, $chain);\n\n    if ( $#array==\
3){return ($pdb_id, \"FIRST\");}\n    elsif ( $#ar\
ray<4){ return ($pdb_id, \"\");}\n    else {return\
 ( $pdb_id, $chain);}\n      \n    \n    \n}\nsub \
get_main_hom_code \n{\n    my $pdb_file=@_[0];\n  \
  my %hom, $n, $best, $best_h;\n    open (F, $pdb_\
file);\n    while (<F>)\n{\n	if ( $_=~/^ATOM/)\n{\\
n	    $h=substr ($_,26, 1);\n	    $n=++$hom{$h};\n\
	    if ($n>$best)\n{\n		$best=$n;\n		$best_h=$h;\\
n}\n}\n}\n    close (F);\n    return $best_h;\n}\n\
\n\nsub get_pdb_file \n{\n    my ($pdb_file_in)=(@\
_);\n    my $result;\n    my @letter;\n    my @cha\
in;\n    my $v;\n    my $pdb_file=$pdb_file_in;\n\\
n    $pdb_file=($pdb_file_in=~/\\S+_S_(\\S+)/)?$1:\
$pdb_file_in;\n    \n    if ($no_remote_pdb_dir==0\
)\n      {\n	$no_remote_pdb_dir=1;\n	$result=get_p\
db_file3 ($pdb_file);\n	$no_remote_pdb_dir=0;\n	if\
 ( $result){return $result;}\n	else\n	  {\n	    \n\
	    lc ($pdb_file);\n	    $result=get_pdb_file3($\
pdb_file);\n	    return  $result;\n	  }\n      }\n\
    else\n      {\n	return get_pdb_file3 ($pdb_fil\
e);\n      }\n    \n  }\n\nsub get_pdb_file3 \n{\n\
    my $pdb_file_in=@_[0];\n    my $result;\n    m\
y @letter;\n    my @chain;\n    my $lcfile;\n    m\
y $ucfile;\n    my $pdb_file=$pdb_file_in;\n    \n\
    $lcfile=lc $pdb_file;\n    $ucfile=uc $pdb_fil\
e;\n\n    if ( ($result=get_pdb_file2 ($pdb_file))\
){return $result;}\n    \n\n    if ($lcfile ne $pd\
b_file && ($result=get_pdb_file2 ($lcfile))){retur\
n $result;}\n    if ($ucfile ne $pdb_file && ($res\
ult=get_pdb_file2 ($ucfile))){return $result;}\n  \
  \n   \n    \n    return \"\";\n}\nsub get_pdb_fi\
le2\n{\n    my $pdb_file=@_[0];\n    my $return_va\
lue;\n    \n    $return_value=\"\";\n    \n    if \
( ($result=get_pdb_file1 ($pdb_file))){$return_val\
ue=$result;}\n    elsif ( !($pdb_file=~/\\.pdb/) &\
& !($pdb_file=~/\\.PDB/))\n{\n	if ( ($result=get_p\
db_file1 (\"$pdb_file.pdb\"))){$return_value=$resu\
lt;}\n	elsif ( ($result=get_pdb_file1 (\"$pdb_file\
.PDB\"))){$return_value=$result;}\n\n	elsif ( ($re\
sult=get_pdb_file1 (\"pdb$pdb_file.pdb\"))){$retur\
n_value=$result;}	\n	elsif ( ($result=get_pdb_file\
1 (\"pdb$pdb_file.PDB\"))){$return_value=$result;}\
\n	elsif ( ($result=get_pdb_file1 (\"PDB$pdb_file.\
PDB\"))){$return_value=$result;}\n	elsif ( ($resul\
t=get_pdb_file1 (\"PDB$pdb_file.pdb\"))){$return_v\
alue=$result;}\n	\n	\n	elsif ( ($result=get_pdb_fi\
le1 (\"$pdb_file.ent\"))){$return_value=$result;}\\
n	elsif ( ($result=get_pdb_file1 (\"pdb$pdb_file.e\
nt\"))){$return_value=$result;}\n	elsif ( ($result\
=get_pdb_file1 (\"PDB$pdb_file.ent\"))){$return_va\
lue=$result;}\n\n	elsif ( ($result=get_pdb_file1 (\
\"$pdb_file.ENT\"))){$return_value=$result;}\n	els\
if ( ($result=get_pdb_file1 (\"pdb$pdb_file.ENT\")\
)){$return_value=$result;}\n	elsif ( ($result=get_\
pdb_file1 (\"PDB$pdb_file.ENT\"))){$return_value=$\
result;}\n	\n	\n	\n}\n    return $return_value;\n}\
\n    \nsub get_pdb_file1\n{\n    my ($pdb_file)=(\
@_);\n    my $return_value;\n    \n\n    $return_v\
alue=\"\";\n    if ( ($result=get_pdb_file0 ($pdb_\
file))){$return_value=$result;}\n    elsif ( ($res\
ult=get_pdb_file0 (\"$pdb_file.Z\"))){$return_valu\
e=$result;}\n    elsif ( ($result=get_pdb_file0 (\\
"$pdb_file.gz\"))){$return_value=$result;}\n    el\
sif ( ($result=get_pdb_file0 (\"$pdb_file.GZ\"))){\
$return_value=$result;}\n    return $return_value;\
\n}\nsub get_pdb_file0 \n{ \n    my ($pdb_file, $a\
ttempt)=(@_);\n    my $pdb_file=@_[0];\n    my $tm\
p_pdb_file;    \n    my $return_value;\n\n    if (\
 !$attempt){$attempt=1;}\n    \n    $local_pdb_fil\
e=\"$pdb_file\";\n    if ( $local_pdb_file eq \"\"\
)\n{\n	$tmp_pdb_file=vtmpnam();\n	open F, \">$tmp_\
pdb_file\";\n	\n	while (<STDIN>){print F \"$_\";}\\
n	close (F);\n	\n	if (-e $tmp_pdb_file && &is_pdb_\
file ( $local_pdb_file))\n{return $tmp_pdb_file;}\\
n}\n\n    $local_pdb_file=\"$pdb_file\";\n    &deb\
ug_print (\"\\nTry access local file: $local_pdb_f\
ile\");\n    \n    $local_pdb_file=&check_pdb_file\
4compression ($local_pdb_file);\n    if ( -e $loca\
l_pdb_file && (&is_pdb_file ($local_pdb_file) || $\
force_pdb))\n{\n	&debug_print ( \"\\n\\tIs in Curr\
ent Dir\");\n	$tmp_pdb_file=vtmpnam();\n	`cp $loca\
l_pdb_file $tmp_pdb_file`;\n	return $tmp_pdb_file;\
\n}\n    else\n{\n	&debug_print (\"\\n\\tFile Not \
in Current Dir\");\n}\n\n    if ($pdb_file=~/^pdb/\
||$pdb_file=~/^PDB/){$pdb_div=substr ($pdb_file, 4\
, 2);}\n    else\n{\n	  $pdb_div=substr ($pdb_file\
, 1, 2);\n}\n    $local_pdb_file=\"$pdb_dir/$pdb_d\
iv/$pdb_file\";\n    $local_pdb_file=&check_pdb_fi\
le4compression ( $local_pdb_file);\n    &debug_pri\
nt (\"\\nTry access file From PDB_DIR: $local_pdb_\
file\");\n    if ($pdb_dir && -e $local_pdb_file &\
& &is_pdb_file ($local_pdb_file))\n{\n	&debug_prin\
t ( \"\\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_fil\
e=vtmpnam();\n	`cp $local_pdb_file $tmp_pdb_file`;\
\n	return $tmp_pdb_file;\n}\n\n    $local_pdb_file\
=\"$pdb_dir/$pdb_file\";\n    $local_pdb_file=&che\
ck_pdb_file4compression ( $local_pdb_file);\n    &\
debug_print (\"\\nTry access file From PDB_DIR: lo\
cal_pdb_file\");\n    if ($pdb_dir && -e $local_pd\
b_file && &is_pdb_file ($local_pdb_file))\n{\n	&de\
bug_print ( \"\\n\\tIs in Local PDB DIR\");\n	$tmp\
_pdb_file=vtmpnam();\n	`cp $local_pdb_file $tmp_pd\
b_file`;\n	return $tmp_pdb_file;\n}\n\n    $local_\
pdb_file=\"$pdb_dir$pdb_file\";\n    $local_pdb_fi\
le=&check_pdb_file4compression ( $local_pdb_file);\
\n    &debug_print (\"\\nTry access file From PDB_\
DIR: $local_pdb_file\");\n    if ($pdb_dir && -e $\
local_pdb_file && &is_pdb_file ($local_pdb_file))\\
n{\n	&debug_print ( \"\\n\\tIs in Local PDB DIR\")\
;\n	$tmp_pdb_file=vtmpnam();\n	`cp $local_pdb_file\
 $tmp_pdb_file`;\n	return $tmp_pdb_file;\n}\n    e\
lse\n{&debug_print ( \"\\n\\tNot In Local Pdb Dir\\
");}\n\n    if ($cache ne \"NO\" && $cache ne \"no\
\")\n{\n\n	$local_pdb_file=\"$cache/$pdb_file\";\n\
	$local_pdb_file=&check_pdb_file4compression ( $lo\
cal_pdb_file);\n	&debug_print(\"\\nTry access file\
 From Cache: $local_pdb_file\");\n	if (-e $local_p\
db_file && &is_pdb_file ($local_pdb_file))\n{\n	  \
  &debug_print ( \"\\n\\tIs in T-Coffee Cache\");\\
n	    $tmp_pdb_file=vtmpnam();\n	    `cp $local_pd\
b_file $tmp_pdb_file`;\n	    return $tmp_pdb_file;\
\n}\n	else{&debug_print ( \"\\n\\tNot in Cache Dir\
\");}\n}\n\nif (!$no_remote_pdb_dir) \n  {\n    my\
 $value=&is_released ($pdb_file);\n    my $return_\
value=\"\";\n    if ($value==1)\n      {\n	\n	&deb\
ug_print (\"\\n***********************************\
******************\\nTry Remote Access for $pdb_fi\
le\");\n	$tmp_pdb_file=vtmpnam();\n	$netcommand=$n\
etaddress;\n	$netcommand=~s/%%/$pdb_file/g;\n	&url\
2file(\"$netcommand\", \"$tmp_pdb_file.$netcompres\
sion\");\n	&debug_print(\"\\nREMOTE: $netcommand\\\
n\");\n	\n	$compressed_tmp_file_name=\"$tmp_pdb_fi\
le.$netcompression\";\n	\n	if ($netcompression && \
-B $compressed_tmp_file_name && $attempt<5)\n	  {\\
n	    my $r;\n	    &debug_print (\"\\n\\tFile Foun\
d Remotely\");\n	    if (($r=safe_system ( \"$netc\
ompression_pg $compressed_tmp_file_name\")!=$EXIT_\
SUCCESS) && $attempts<5)\n	      {\n		&debug_print\
 (\"\\n\\tProper Download Failed Try again\");\n		\
unlink $compressed_tmp_file_name;\n		print \"\\nFa\
iled to Download $compressed_tmp_file_name. New At\
tempt $attempt/5\\n\";\n		return &get_pdb_file0($p\
db_file, $attempt+1);\n	      }\n	    elsif ($r== \
$EXIT_SUCCESS)\n	      {\n		&debug_print (\"\\n\\t\
Proper Download Succeeded \");\n		$return_value=$t\
mp_pdb_file;\n	      }\n	    else\n	      {\n		&de\
bug_print (\"\\n\\tProper Download Failed \");\n		\
&debug_print (\"\\nFile Not Found Remotely\");\n		\
unlink $compressed_tmp_file_name;\n	      }\n	  }\\
n	else\n	  {\n\n	    &debug_print (\"\\nFile Not F\
ound Remotely\");\n	    unlink $compressed_tmp_fil\
e_name;\n	  }\n	#Update cache if required\n	if ($c\
ache ne \"no\" && $cache ne \"update\" && -e $retu\
rn_value)\n	  {\n	    `cp $return_value $cache/$pd\
b_file.pdb`;\n	    #`t_coffee -other_pg clean_cach\
e.pl -file $pdb_file.pdb -dir $cache`;\n	  }\n    \
  }\n    &debug_print (\"\\nRemote Download Finish\
ed\");\n    return $return_value;\n  }\nreturn \"\\
";\n}\n\nsub check_pdb_file4compression \n{\n    m\
y $file=@_[0];\n    my $tmp;\n    my $r;\n    \n  \
  $tmp=&vtmpnam();\n    if (-e $tmp){unlink $tmp;}\
\n    \n    $file=~s/\\/\\//\\//g;\n    if    (-B \
$file && ($file=~/\\.Z/)) {`cp $file $tmp.Z`;`rm $\
tmp`;`gunzip $tmp.Z $SILENT`;$r=$tmp;}\n    elsif \
(-B $file && ($file=~/\\.gz/)){`cp $file $tmp.gz`;\
`gunzip $tmp.gz $SILENT`;return $r=$tmp;}\n    els\
if (-B $file ){`cp $file $tmp.gz`;`gunzip $tmp.gz \
$SILENT`;$r=$tmp;}\n    elsif ( -e $file ) {$r= $f\
ile;}\n    elsif ( -e \"$file.gz\" ){ `cp $file.gz\
 $tmp.gz`;`gunzip     $tmp.gz $SILENT`;$r=$tmp;}  \
  \n    elsif ( -e \"$file.Z\") {`cp $file.Z  $tmp\
.Z`; `gunzip $tmp.Z $SILENT`;$r=$tmp;}\n    else  \
{$r= $file;}\n\n    if ( -e \"$tmp.Z\"){unlink \"$\
tmp.Z\";}\n    if ( -e \"$tmp.gz\"){unlink \"$tmp.\
gz\";}\n    \n    return $r;\n    \n}\n\n\n\n\n\n \
   \n\n\n\n\n\n\n\nsub vfopen \n{\n    my $file=@_\
[0];\n    my $mode=@_[1];\n    my $tmp;\n    my $F\
 = new FileHandle;\n    \n    \n    $tmp=$file;\n	\
\n    \n    if ( $mode eq \"r\" && !-e $file){ mye\
xit(flush_error (\"Cannot open file $file\"));}\n \
   elsif ($mode eq \"w\"){$tmp=\">$file\";}\n    e\
lsif ($mode eq \"a\"){$tmp=\">>$file\";}\n    \n  \
  \n    open ($F,$tmp);\n    return $F;\n}\nsub de\
bug_print\n{\n    my $message =@_[0];\n    if ($de\
bug){print STDERR \"NO_REMOTE_PDB_DIR: $no_remote_\
pdb_dir - $message [DEBUG:extract_from_pdb]\";}\n \
   return;\n}\nsub is_aa \n{\n    my ($aa, $chain)\
 =@_;\n\n    my $one;\n    my $trhee;\n    \n    i\
f ( $onelett{$molecule_type{$chain}}->{$aa} eq 'X'\
 || !$onelett{$molecule_type{$chain}}->{$aa} ){ret\
urn '';}\n    else\n      {\n	$one=$onelett{$molec\
ule_type{$chain}}->{$aa};\n\n	$three=$threelett{$m\
olecule_type{$chain}}->{$one};\n	\n\n	return $thre\
e;\n      }\n  }\n\n\n\n\n\nsub url2file\n{\n    m\
y ($address, $out, $wget_arg, $curl_arg)=(@_);\n  \
  my ($pg, $flag, $r, $arg, $count);\n    \n    if\
 (!$CONFIGURATION){&check_configuration (\"wget\",\
 \"INTERNET\", \"gzip\");$CONFIGURATION=1;}\n    \\
n    if (&pg_is_installed (\"wget\"))   {$pg=\"wge\
t\"; $flag=\"-O\";$arg=$wget_arg;}\n    elsif (&pg\
_is_installed (\"curl\")){$pg=\"curl\"; $flag=\"-o\
\";$arg=$curl_arg;}\n    return safe_system (\"$pg\
 $flag$out $address >/dev/null 2>/dev/null\");\n\n\
}\n\n\n\n\nsub pdbfile2chaintype\n  {\n    my $fil\
e=@_[0];\n    my %ct;\n    my $F;\n    \n    $F=vf\
open ($file, \"r\");\n    while (<$F>)\n      {\n	\
my $line=$_;\n	if ($line =~/^ATOM/)\n	  {\n	    my\
 $C=substr($line,21,1);\n	    if (!$ct{$C})\n	    \
  {	\n		my $r=substr($line,17,3);\n		$r=~s/\\s+//;\
\n		if (length ($r)==1){$ct{$C}=\"R\";}\n		elsif (\
length ($r)==2){$ct{$C}=\"D\";}\n		elsif (length (\
$r)==3){$ct{$C}=\"P\";}\n		else \n		  {\n		    mye\
xit(flush_error(\"ERROR: Could not read RES_ID fie\
ld in file $file\"));\n		  }\n	      }\n	  }\n    \
  }\n    close ($F);\n    return %ct;\n  }\n   \n \
   \n\n\n\nsub fill_threelett_RNA\n{\n\n	my %three\
lett=(\n	'A', '  A',\n	'T', '  T',\n	'U', '  U',\n\
	'C', '  C',\n	'G', '  G',\n	'I', '  I', #Inosine\\
n	);\n	\n	return %threelett;\n\n}\n\n\nsub fill_on\
elett_RNA\n{\n	my   %onelett=(\n	'  A' => 'A',\n	'\
  T' => 'T',\n	'  U' => 'U',\n	'  C' => 'C',\n	'  \
G' => 'G',\n	'CSL' => 'X',\n	'UMS' => 'X',\n	'  I'\
 => 'I',\n	'A' => 'A',\n	'T' => 'T',\n	'U' => 'U',\
\n	'C' => 'C',\n	'G' => 'G',\n	'I' => 'I',\n	);\n\\
n	return %onelett;\n\n}\n\n\nsub fill_onelett_DNA\\
n{\n	my   %onelett=(\n	' DA', 'A',\n	' DT', 'T',\n\
	' DC', 'C',\n	' DG', 'G',\n	'DA', 'A',\n	'DT', 'T\
',\n	'DC', 'C',\n	'DG', 'G',\n	);\n\n	return %onel\
ett;\n\n}\n\nsub fill_threelett_DNA\n{\n\n	my %thr\
eelett=(\n	'A', ' DA',\n	'T', ' DT',\n	'C', ' DC',\
\n	'G', ' DG',\n	);\n\n	return %threelett;\n\n}\n\\
n\n\n\nsub fill_threelett_prot\n{  \n  my %threele\
tt;\n\n  %threelett=(\n'A', 'ALA',\n'C', 'CYS',\n'\
D', 'ASP',\n'E', 'GLU',\n'F', 'PHE',\n'G', 'GLY',\\
n'H', 'HIS',\n'I', 'ILE',\n'K', 'LYS',\n'L', 'LEU'\
,\n'N', 'ASN',\n'M', 'MET',\n'P', 'PRO',\n'Q', 'GL\
N',\n'R', 'ARG',\n'S', 'SER',\n'T', 'THR',\n'V', '\
VAL',\n'W', 'TRP',\n'Y', 'TYR',\n);\n\nreturn %thr\
eelett;\n\n\n}\n\nsub fill_onelett_prot\n{\n    my\
 %onelett;\n    \n    %onelett=(\n\n'10A', 'X',\n'\
11O', 'X',\n'12A', 'X',\n'13P', 'X',\n'13R', 'X',\\
n'13S', 'X',\n'14W', 'X',\n'15P', 'X',\n'16A', 'X'\
,\n'16G', 'X',\n'1AN', 'X',\n'1AP', 'X',\n'1AR', '\
X',\n'1BH', 'X',\n'1BO', 'X',\n'1C5', 'X',\n'1CU',\
 'X',\n'1DA', 'X',\n'1GL', 'X',\n'1GN', 'X',\n'1IN\
', 'X',\n'1LU', 'L',\n'1MA', 'X',\n'1MC', 'X',\n'1\
MG', 'X',\n'1MZ', 'X',\n'1NA', 'X',\n'1NB', 'X',\n\
'1NI', 'X',\n'1PA', 'A',\n'1PC', 'X',\n'1PE', 'X',\
\n'1PG', 'X',\n'1PI', 'A',\n'1PM', 'X',\n'1PN', 'X\
',\n'1PU', 'X',\n'1PY', 'X',\n'1UN', 'X',\n'24T', \
'X',\n'25T', 'X',\n'26P', 'X',\n'2AB', 'X',\n'2AM'\
, 'X',\n'2AN', 'X',\n'2AP', 'X',\n'2AR', 'X',\n'2A\
S', 'D',\n'2BL', 'X',\n'2BM', 'X',\n'2CP', 'X',\n'\
2DA', 'X',\n'2DG', 'X',\n'2DP', 'X',\n'2DT', 'X',\\
n'2EP', 'X',\n'2EZ', 'X',\n'2FG', 'X',\n'2FL', 'X'\
,\n'2FP', 'X',\n'2FU', 'X',\n'2GL', 'X',\n'2GP', '\
X',\n'2HP', 'X',\n'2IB', 'X',\n'2IP', 'X',\n'2LU',\
 'L',\n'2MA', 'X',\n'2MD', 'X',\n'2ME', 'X',\n'2MG\
', 'X',\n'2ML', 'L',\n'2MO', 'X',\n'2MR', 'R',\n'2\
MU', 'X',\n'2MZ', 'X',\n'2NO', 'X',\n'2NP', 'X',\n\
'2OG', 'X',\n'2PA', 'X',\n'2PC', 'X',\n'2PE', 'X',\
\n'2PG', 'X',\n'2PH', 'X',\n'2PI', 'X',\n'2PL', 'X\
',\n'2PP', 'X',\n'2PU', 'X',\n'2SI', 'X',\n'2TB', \
'X',\n'34C', 'X',\n'35G', 'X',\n'3AA', 'X',\n'3AD'\
, 'X',\n'3AH', 'H',\n'3AN', 'X',\n'3AP', 'X',\n'3A\
T', 'X',\n'3BT', 'X',\n'3CH', 'X',\n'3CN', 'X',\n'\
3CO', 'X',\n'3CP', 'X',\n'3DR', 'X',\n'3EP', 'X',\\
n'3FM', 'X',\n'3GA', 'X',\n'3GP', 'X',\n'3HB', 'X'\
,\n'3HC', 'X',\n'3HP', 'X',\n'3IB', 'X',\n'3ID', '\
X',\n'3IN', 'X',\n'3MA', 'X',\n'3MB', 'X',\n'3MC',\
 'X',\n'3MD', 'D',\n'3MF', 'X',\n'3MP', 'X',\n'3MT\
', 'X',\n'3OL', 'X',\n'3PA', 'X',\n'3PG', 'X',\n'3\
PO', 'X',\n'3PP', 'X',\n'3PY', 'X',\n'49A', 'X',\n\
'4AB', 'X',\n'4AM', 'X',\n'4AN', 'X',\n'4AP', 'X',\
\n'4BA', 'X',\n'4BT', 'X',\n'4CA', 'X',\n'4CO', 'X\
',\n'4HP', 'X',\n'4IP', 'X',\n'4MO', 'X',\n'4MV', \
'X',\n'4MZ', 'X',\n'4NC', 'X',\n'4NP', 'X',\n'4OX'\
, 'X',\n'4PB', 'X',\n'4PN', 'X',\n'4PP', 'X',\n'4S\
C', 'X',\n'4SU', 'X',\n'4TB', 'X',\n'55C', 'X',\n'\
5AD', 'X',\n'5AN', 'X',\n'5AT', 'X',\n'5CM', 'X',\\
n'5GP', 'X',\n'5HP', 'E',\n'5HT', 'X',\n'5IT', 'X'\
,\n'5IU', 'X',\n'5MB', 'X',\n'5MC', 'X',\n'5MD', '\
X',\n'5MP', 'X',\n'5MU', 'X',\n'5NC', 'X',\n'5OB',\
 'X',\n'5PA', 'X',\n'5PV', 'X',\n'6AB', 'X',\n'6CT\
', 'X',\n'6HA', 'X',\n'6HC', 'X',\n'6HG', 'X',\n'6\
HT', 'X',\n'6IN', 'X',\n'6MO', 'X',\n'6MP', 'X',\n\
'6PG', 'X',\n'6WO', 'X',\n'70U', 'X',\n'7DG', 'X',\
\n'7HP', 'X',\n'7I2', 'X',\n'7MG', 'X',\n'7MQ', 'X\
',\n'7NI', 'X',\n'87Y', 'X',\n'8AD', 'X',\n'8BR', \
'X',\n'8IG', 'X',\n'8IN', 'X',\n'8OG', 'X',\n'95A'\
, 'X',\n'9AD', 'X',\n'9AM', 'X',\n'9AP', 'X',\n'9D\
G', 'X',\n'9DI', 'X',\n'9HX', 'X',\n'9OH', 'X',\n'\
9TA', 'X',\n'A12', 'X',\n'A15', 'X',\n'A23', 'X',\\
n'A24', 'X',\n'A26', 'X',\n'A2G', 'X',\n'A2P', 'X'\
,\n'A32', 'X',\n'A3P', 'X',\n'A4P', 'X',\n'A5P', '\
X',\n'A70', 'X',\n'A76', 'X',\n'A77', 'X',\n'A78',\
 'X',\n'A79', 'X',\n'A80', 'X',\n'A85', 'X',\n'A88\
', 'X',\n'A9A', 'X',\n'AA3', 'X',\n'AA4', 'X',\n'A\
A6', 'X',\n'AAA', 'X',\n'AAB', 'X',\n'AAC', 'X',\n\
'AAE', 'X',\n'AAG', 'R',\n'AAH', 'X',\n'AAM', 'X',\
\n'AAN', 'X',\n'AAP', 'X',\n'AAR', 'R',\n'AAS', 'X\
',\n'AAT', 'X',\n'ABA', 'X',\n'ABC', 'X',\n'ABD', \
'X',\n'ABE', 'X',\n'ABH', 'X',\n'ABI', 'X',\n'ABK'\
, 'X',\n'ABM', 'X',\n'ABN', 'X',\n'ABP', 'X',\n'AB\
R', 'X',\n'ABS', 'X',\n'ABU', 'X',\n'AC1', 'X',\n'\
AC2', 'X',\n'ACA', 'X',\n'ACB', 'D',\n'ACC', 'C',\\
n'ACD', 'X',\n'ACE', 'X',\n'ACH', 'X',\n'ACI', 'X'\
,\n'ACL', 'R',\n'ACM', 'X',\n'ACN', 'X',\n'ACO', '\
X',\n'ACP', 'X',\n'ACQ', 'X',\n'ACR', 'X',\n'ACS',\
 'X',\n'ACT', 'X',\n'ACV', 'V',\n'ACX', 'X',\n'ACY\
', 'X',\n'AD2', 'X',\n'AD3', 'X',\n'ADC', 'X',\n'A\
DD', 'X',\n'ADE', 'X',\n'ADH', 'X',\n'ADI', 'X',\n\
'ADM', 'X',\n'ADN', 'X',\n'ADP', 'X',\n'ADQ', 'X',\
\n'ADR', 'X',\n'ADS', 'X',\n'ADT', 'X',\n'ADU', 'X\
',\n'ADW', 'X',\n'ADX', 'X',\n'AE2', 'X',\n'AEA', \
'X',\n'AEB', 'X',\n'AEI', 'D',\n'AEN', 'X',\n'AET'\
, 'T',\n'AF1', 'X',\n'AF3', 'X',\n'AFA', 'D',\n'AF\
P', 'X',\n'AG7', 'X',\n'AGB', 'X',\n'AGF', 'X',\n'\
AGL', 'X',\n'AGM', 'R',\n'AGN', 'X',\n'AGP', 'X',\\
n'AGS', 'X',\n'AGU', 'X',\n'AH0', 'X',\n'AH1', 'X'\
,\n'AHA', 'X',\n'AHB', 'D',\n'AHC', 'X',\n'AHF', '\
X',\n'AHG', 'X',\n'AHH', 'X',\n'AHM', 'X',\n'AHO',\
 'X',\n'AHP', 'X',\n'AHS', 'X',\n'AHT', 'Y',\n'AHU\
', 'X',\n'AHX', 'X',\n'AI1', 'X',\n'AI2', 'X',\n'A\
IB', 'X',\n'AIC', 'X',\n'AIM', 'X',\n'AIP', 'X',\n\
'AIQ', 'X',\n'AIR', 'X',\n'AJ3', 'X',\n'AKB', 'X',\
\n'AKG', 'X',\n'AKR', 'X',\n'AL1', 'X',\n'AL2', 'X\
',\n'AL3', 'X',\n'AL4', 'X',\n'AL5', 'X',\n'AL6', \
'X',\n'AL7', 'X',\n'AL8', 'X',\n'AL9', 'X',\n'ALA'\
, 'A',\n'ALB', 'X',\n'ALC', 'X',\n'ALD', 'L',\n'AL\
E', 'X',\n'ALF', 'X',\n'ALG', 'X',\n'ALL', 'X',\n'\
ALM', 'A',\n'ALN', 'A',\n'ALO', 'T',\n'ALP', 'X',\\
n'ALQ', 'X',\n'ALR', 'X',\n'ALS', 'X',\n'ALT', 'A'\
,\n'ALY', 'K',\n'ALZ', 'X',\n'AMA', 'X',\n'AMB', '\
X',\n'AMC', 'X',\n'AMD', 'X',\n'AMG', 'X',\n'AMH',\
 'X',\n'AMI', 'X',\n'AML', 'X',\n'AMN', 'X',\n'AMO\
', 'X',\n'AMP', 'X',\n'AMQ', 'X',\n'AMR', 'X',\n'A\
MS', 'X',\n'AMT', 'X',\n'AMU', 'X',\n'AMW', 'X',\n\
'AMX', 'X',\n'AMY', 'X',\n'ANA', 'X',\n'ANB', 'X',\
\n'ANC', 'X',\n'AND', 'X',\n'ANE', 'X',\n'ANI', 'X\
',\n'ANL', 'X',\n'ANO', 'X',\n'ANP', 'X',\n'ANS', \
'X',\n'ANT', 'X',\n'AOE', 'X',\n'AOP', 'X',\n'AP1'\
, 'X',\n'AP2', 'X',\n'AP3', 'X',\n'AP4', 'X',\n'AP\
5', 'X',\n'AP6', 'X',\n'APA', 'X',\n'APB', 'X',\n'\
APC', 'X',\n'APE', 'F',\n'APF', 'X',\n'APG', 'X',\\
n'APH', 'A',\n'API', 'X',\n'APL', 'X',\n'APM', 'X'\
,\n'APN', 'G',\n'APP', 'X',\n'APQ', 'X',\n'APR', '\
X',\n'APS', 'X',\n'APT', 'X',\n'APU', 'X',\n'APX',\
 'X',\n'APY', 'X',\n'APZ', 'X',\n'AQS', 'X',\n'AR1\
', 'X',\n'AR2', 'X',\n'ARA', 'X',\n'ARB', 'X',\n'A\
RC', 'X',\n'ARD', 'X',\n'ARG', 'R',\n'ARH', 'X',\n\
'ARI', 'X',\n'ARM', 'R',\n'ARN', 'X',\n'ARO', 'R',\
\n'ARP', 'X',\n'ARQ', 'X',\n'ARS', 'X',\n'AS1', 'R\
',\n'AS2', 'X',\n'ASA', 'D',\n'ASB', 'D',\n'ASC', \
'X',\n'ASD', 'X',\n'ASE', 'X',\n'ASF', 'X',\n'ASI'\
, 'X',\n'ASK', 'D',\n'ASL', 'X',\n'ASM', 'N',\n'AS\
O', 'X',\n'ASP', 'D',\n'ASQ', 'X',\n'ASU', 'X',\n'\
ATA', 'X',\n'ATC', 'X',\n'ATD', 'X',\n'ATF', 'X',\\
n'ATG', 'X',\n'ATH', 'X',\n'ATM', 'X',\n'ATO', 'X'\
,\n'ATP', 'X',\n'ATQ', 'X',\n'ATR', 'X',\n'ATT', '\
X',\n'ATY', 'X',\n'ATZ', 'X',\n'AUC', 'X',\n'AUR',\
 'X',\n'AVG', 'X',\n'AXP', 'X',\n'AYA', 'A',\n'AZ2\
', 'X',\n'AZA', 'X',\n'AZC', 'X',\n'AZD', 'X',\n'A\
ZE', 'X',\n'AZI', 'X',\n'AZL', 'X',\n'AZM', 'X',\n\
'AZR', 'X',\n'AZT', 'X',\n'B12', 'X',\n'B1F', 'F',\
\n'B2A', 'A',\n'B2F', 'F',\n'B2I', 'I',\n'B2V', 'V\
',\n'B3I', 'X',\n'B3P', 'X',\n'B7G', 'X',\n'B96', \
'X',\n'B9A', 'X',\n'BA1', 'X',\n'BAA', 'X',\n'BAB'\
, 'X',\n'BAC', 'X',\n'BAF', 'X',\n'BAH', 'X',\n'BA\
I', 'X',\n'BAK', 'X',\n'BAL', 'A',\n'BAM', 'X',\n'\
BAO', 'X',\n'BAP', 'X',\n'BAR', 'X',\n'BAS', 'X',\\
n'BAT', 'F',\n'BAY', 'X',\n'BAZ', 'X',\n'BB1', 'X'\
,\n'BB2', 'X',\n'BBA', 'X',\n'BBH', 'X',\n'BBS', '\
X',\n'BBT', 'X',\n'BBZ', 'X',\n'BCA', 'X',\n'BCB',\
 'X',\n'BCC', 'X',\n'BCD', 'X',\n'BCL', 'X',\n'BCN\
', 'X',\n'BCR', 'X',\n'BCS', 'C',\n'BCT', 'X',\n'B\
CY', 'X',\n'BCZ', 'X',\n'BDA', 'X',\n'BDG', 'X',\n\
'BDK', 'X',\n'BDM', 'X',\n'BDN', 'X',\n'BDS', 'X',\
\n'BE1', 'X',\n'BE2', 'X',\n'BEA', 'X',\n'BEF', 'X\
',\n'BEN', 'X',\n'BEO', 'X',\n'BEP', 'X',\n'BER', \
'X',\n'BES', 'X',\n'BET', 'X',\n'BEZ', 'X',\n'BF2'\
, 'X',\n'BFA', 'X',\n'BFD', 'X',\n'BFP', 'X',\n'BF\
S', 'X',\n'BFU', 'X',\n'BG6', 'X',\n'BGF', 'X',\n'\
BGG', 'X',\n'BGL', 'X',\n'BGN', 'X',\n'BGP', 'X',\\
n'BGX', 'X',\n'BH4', 'X',\n'BHA', 'X',\n'BHC', 'X'\
,\n'BHD', 'D',\n'BHO', 'X',\n'BHS', 'X',\n'BIC', '\
X',\n'BIN', 'X',\n'BIO', 'X',\n'BIP', 'X',\n'BIS',\
 'X',\n'BIZ', 'X',\n'BJH', 'X',\n'BJI', 'X',\n'BJP\
', 'X',\n'BLA', 'X',\n'BLB', 'X',\n'BLE', 'L',\n'B\
LG', 'P',\n'BLI', 'X',\n'BLM', 'X',\n'BLV', 'X',\n\
'BLY', 'K',\n'BM1', 'X',\n'BM2', 'X',\n'BM5', 'X',\
\n'BM9', 'X',\n'BMA', 'X',\n'BMD', 'X',\n'BME', 'X\
',\n'BMP', 'X',\n'BMQ', 'X',\n'BMS', 'X',\n'BMT', \
'T',\n'BMU', 'X',\n'BMY', 'X',\n'BMZ', 'X',\n'BNA'\
, 'X',\n'BNG', 'X',\n'BNI', 'X',\n'BNN', 'F',\n'BN\
O', 'L',\n'BNS', 'X',\n'BNZ', 'X',\n'BO3', 'X',\n'\
BO4', 'X',\n'BOC', 'X',\n'BOG', 'X',\n'BOM', 'X',\\
n'BOT', 'X',\n'BOX', 'X',\n'BOZ', 'X',\n'BPA', 'X'\
,\n'BPB', 'X',\n'BPD', 'X',\n'BPG', 'X',\n'BPH', '\
X',\n'BPI', 'X',\n'BPJ', 'X',\n'BPM', 'X',\n'BPN',\
 'X',\n'BPO', 'X',\n'BPP', 'X',\n'BPT', 'X',\n'BPY\
', 'X',\n'BRB', 'X',\n'BRC', 'X',\n'BRE', 'X',\n'B\
RI', 'X',\n'BRL', 'X',\n'BRM', 'X',\n'BRN', 'X',\n\
'BRO', 'X',\n'BRS', 'X',\n'BRU', 'X',\n'BRZ', 'X',\
\n'BSB', 'X',\n'BSI', 'X',\n'BSP', 'X',\n'BT1', 'X\
',\n'BT2', 'X',\n'BT3', 'X',\n'BTA', 'L',\n'BTB', \
'X',\n'BTC', 'C',\n'BTD', 'X',\n'BTN', 'X',\n'BTP'\
, 'X',\n'BTR', 'W',\n'BU1', 'X',\n'BUA', 'X',\n'BU\
B', 'X',\n'BUC', 'X',\n'BUG', 'X',\n'BUL', 'X',\n'\
BUM', 'X',\n'BUQ', 'X',\n'BUT', 'X',\n'BVD', 'X',\\
n'BX3', 'X',\n'BYS', 'X',\n'BZ1', 'X',\n'BZA', 'X'\
,\n'BZB', 'X',\n'BZC', 'X',\n'BZD', 'X',\n'BZF', '\
X',\n'BZI', 'X',\n'BZM', 'X',\n'BZO', 'X',\n'BZP',\
 'X',\n'BZQ', 'X',\n'BZS', 'X',\n'BZT', 'X',\n'C02\
', 'X',\n'C11', 'X',\n'C1O', 'X',\n'C20', 'X',\n'C\
24', 'X',\n'C2F', 'X',\n'C2O', 'X',\n'C2P', 'X',\n\
'C3M', 'X',\n'C3P', 'X',\n'C3X', 'X',\n'C48', 'X',\
\n'C4M', 'X',\n'C4X', 'X',\n'C5C', 'X',\n'C5M', 'X\
',\n'C5P', 'X',\n'C5X', 'X',\n'C60', 'X',\n'C6C', \
'X',\n'C6M', 'X',\n'C78', 'X',\n'C8E', 'X',\n'CA3'\
, 'X',\n'CA5', 'X',\n'CAA', 'X',\n'CAB', 'X',\n'CA\
C', 'X',\n'CAD', 'X',\n'CAF', 'C',\n'CAG', 'X',\n'\
CAH', 'X',\n'CAL', 'X',\n'CAM', 'X',\n'CAN', 'X',\\
n'CAO', 'X',\n'CAP', 'X',\n'CAQ', 'X',\n'CAR', 'X'\
,\n'CAS', 'C',\n'CAT', 'X',\n'CAV', 'X',\n'CAY', '\
C',\n'CAZ', 'X',\n'CB3', 'X',\n'CB4', 'X',\n'CBA',\
 'X',\n'CBD', 'X',\n'CBG', 'X',\n'CBI', 'X',\n'CBL\
', 'X',\n'CBM', 'X',\n'CBN', 'X',\n'CBO', 'X',\n'C\
BP', 'X',\n'CBS', 'X',\n'CBX', 'X',\n'CBZ', 'X',\n\
'CC0', 'X',\n'CC1', 'X',\n'CCC', 'X',\n'CCH', 'X',\
\n'CCI', 'X',\n'CCM', 'X',\n'CCN', 'X',\n'CCO', 'X\
',\n'CCP', 'X',\n'CCR', 'X',\n'CCS', 'C',\n'CCV', \
'X',\n'CCY', 'X',\n'CD1', 'X',\n'CDC', 'X',\n'CDE'\
, 'X',\n'CDF', 'X',\n'CDI', 'X',\n'CDL', 'X',\n'CD\
M', 'X',\n'CDP', 'X',\n'CDR', 'X',\n'CDU', 'X',\n'\
CE1', 'X',\n'CEA', 'C',\n'CEB', 'X',\n'CEC', 'X',\\
n'CED', 'X',\n'CEF', 'X',\n'CEH', 'X',\n'CEM', 'X'\
,\n'CEO', 'X',\n'CEP', 'X',\n'CEQ', 'X',\n'CER', '\
X',\n'CES', 'G',\n'CET', 'X',\n'CFC', 'X',\n'CFF',\
 'X',\n'CFM', 'X',\n'CFO', 'X',\n'CFP', 'X',\n'CFS\
', 'X',\n'CFX', 'X',\n'CGN', 'X',\n'CGP', 'X',\n'C\
GS', 'X',\n'CGU', 'E',\n'CH2', 'X',\n'CH3', 'X',\n\
'CHA', 'X',\n'CHB', 'X',\n'CHD', 'X',\n'CHF', 'X',\
\n'CHG', 'G',\n'CHI', 'X',\n'CHN', 'X',\n'CHO', 'X\
',\n'CHP', 'G',\n'CHR', 'X',\n'CHS', 'F',\n'CHT', \
'X',\n'CHX', 'X',\n'CIC', 'X',\n'CIN', 'X',\n'CIP'\
, 'X',\n'CIR', 'X',\n'CIT', 'X',\n'CIU', 'X',\n'CK\
I', 'X',\n'CL1', 'X',\n'CL2', 'X',\n'CLA', 'X',\n'\
CLB', 'A',\n'CLC', 'S',\n'CLD', 'A',\n'CLE', 'L',\\
n'CLF', 'X',\n'CLK', 'S',\n'CLL', 'X',\n'CLM', 'X'\
,\n'CLN', 'X',\n'CLO', 'X',\n'CLP', 'X',\n'CLQ', '\
X',\n'CLR', 'X',\n'CLS', 'X',\n'CLT', 'X',\n'CLX',\
 'X',\n'CLY', 'X',\n'CMA', 'R',\n'CMC', 'X',\n'CMD\
', 'X',\n'CME', 'C',\n'CMG', 'X',\n'CMK', 'X',\n'C\
MN', 'X',\n'CMO', 'X',\n'CMP', 'X',\n'CMR', 'X',\n\
'CMS', 'X',\n'CMT', 'C',\n'CMX', 'X',\n'CNA', 'X',\
\n'CNC', 'X',\n'CND', 'X',\n'CNH', 'X',\n'CNM', 'X\
',\n'CNN', 'X',\n'CNO', 'X',\n'CNP', 'X',\n'CO2', \
'X',\n'CO3', 'X',\n'CO5', 'X',\n'CO8', 'X',\n'COA'\
, 'X',\n'COB', 'X',\n'COC', 'X',\n'COD', 'X',\n'CO\
E', 'X',\n'COF', 'X',\n'COH', 'X',\n'COI', 'X',\n'\
COJ', 'X',\n'COL', 'X',\n'COM', 'X',\n'CON', 'X',\\
n'COP', 'X',\n'COR', 'X',\n'COS', 'X',\n'COT', 'X'\
,\n'COY', 'X',\n'CP1', 'G',\n'CP2', 'X',\n'CP4', '\
X',\n'CPA', 'X',\n'CPB', 'X',\n'CPC', 'X',\n'CPD',\
 'X',\n'CPG', 'X',\n'CPH', 'X',\n'CPI', 'X',\n'CPM\
', 'X',\n'CPN', 'G',\n'CPO', 'X',\n'CPP', 'X',\n'C\
PQ', 'X',\n'CPR', 'X',\n'CPS', 'X',\n'CPT', 'X',\n\
'CPU', 'X',\n'CPV', 'X',\n'CPY', 'X',\n'CR1', 'X',\
\n'CR6', 'X',\n'CRA', 'X',\n'CRB', 'X',\n'CRC', 'X\
',\n'CRG', 'X',\n'CRH', 'X',\n'CRO', 'T',\n'CRP', \
'X',\n'CRQ', 'X',\n'CRS', 'X',\n'CRT', 'X',\n'CRY'\
, 'X',\n'CSA', 'C',\n'CSB', 'X',\n'CSD', 'C',\n'CS\
E', 'C',\n'CSH', 'X',\n'CSI', 'X',\n'CSN', 'X',\n'\
CSO', 'C',\n'CSP', 'C',\n'CSR', 'C',\n'CSS', 'C',\\
n'CST', 'X',\n'CSW', 'C',\n'CSX', 'C',\n'CSY', 'X'\
,\n'CSZ', 'C',\n'CT3', 'X',\n'CTA', 'X',\n'CTB', '\
X',\n'CTC', 'X',\n'CTD', 'X',\n'CTH', 'T',\n'CTO',\
 'X',\n'CTP', 'X',\n'CTR', 'X',\n'CTS', 'X',\n'CTT\
', 'X',\n'CTY', 'X',\n'CTZ', 'X',\n'CU1', 'X',\n'C\
UA', 'X',\n'CUC', 'X',\n'CUL', 'X',\n'CUO', 'X',\n\
'CUZ', 'X',\n'CVI', 'X',\n'CXF', 'X',\n'CXL', 'X',\
\n'CXM', 'M',\n'CXN', 'X',\n'CXP', 'X',\n'CXS', 'X\
',\n'CY1', 'C',\n'CY3', 'X',\n'CYB', 'X',\n'CYC', \
'X',\n'CYF', 'C',\n'CYG', 'C',\n'CYH', 'X',\n'CYL'\
, 'X',\n'CYM', 'C',\n'CYN', 'X',\n'CYO', 'X',\n'CY\
P', 'X',\n'CYQ', 'C',\n'CYS', 'C',\n'CYU', 'X',\n'\
CYY', 'X',\n'CYZ', 'X',\n'CZH', 'X',\n'CZZ', 'C',\\
n'D12', 'X',\n'D13', 'X',\n'D16', 'X',\n'D18', 'X'\
,\n'D19', 'X',\n'D1P', 'X',\n'D24', 'X',\n'D34', '\
X',\n'D35', 'X',\n'D4D', 'X',\n'D4T', 'X',\n'D6G',\
 'X',\n'DA2', 'R',\n'DA3', 'X',\n'DA6', 'X',\n'DA7\
', 'X',\n'DAA', 'X',\n'DAB', 'X',\n'DAC', 'X',\n'D\
AD', 'X',\n'DAE', 'X',\n'DAF', 'X',\n'DAG', 'X',\n\
'DAH', 'A',\n'DAJ', 'X',\n'DAK', 'X',\n'DAL', 'A',\
\n'DAM', 'A',\n'DAN', 'X',\n'DAO', 'X',\n'DAP', 'X\
',\n'DAQ', 'X',\n'DAR', 'R',\n'DAS', 'D',\n'DAT', \
'X',\n'DAU', 'X',\n'DAV', 'X',\n'DBA', 'X',\n'DBD'\
, 'X',\n'DBF', 'X',\n'DBG', 'X',\n'DBI', 'X',\n'DB\
V', 'X',\n'DBY', 'Y',\n'DCA', 'X',\n'DCB', 'X',\n'\
DCE', 'X',\n'DCF', 'X',\n'DCG', 'X',\n'DCH', 'X',\\
n'DCI', 'I',\n'DCL', 'X',\n'DCM', 'X',\n'DCP', 'X'\
,\n'DCS', 'X',\n'DCT', 'X',\n'DCY', 'C',\n'DCZ', '\
X',\n'DDA', 'X',\n'DDB', 'X',\n'DDC', 'X',\n'DDF',\
 'X',\n'DDG', 'X',\n'DDH', 'X',\n'DDL', 'X',\n'DDM\
', 'X',\n'DDO', 'L',\n'DDP', 'X',\n'DDQ', 'X',\n'D\
DT', 'Y',\n'DDU', 'X',\n'DEA', 'X',\n'DEB', 'X',\n\
'DEC', 'X',\n'DEF', 'X',\n'DEL', 'X',\n'DEM', 'X',\
\n'DEN', 'X',\n'DEP', 'X',\n'DEQ', 'X',\n'DES', 'X\
',\n'DET', 'X',\n'DFC', 'X',\n'DFG', 'X',\n'DFI', \
'X',\n'DFL', 'X',\n'DFO', 'X',\n'DFP', 'X',\n'DFR'\
, 'X',\n'DFT', 'X',\n'DFV', 'X',\n'DFX', 'X',\n'DG\
2', 'X',\n'DG3', 'X',\n'DG6', 'X',\n'DGA', 'X',\n'\
DGD', 'X',\n'DGG', 'X',\n'DGL', 'E',\n'DGN', 'Q',\\
n'DGP', 'X',\n'DGT', 'X',\n'DGX', 'X',\n'DH2', 'X'\
,\n'DHA', 'A',\n'DHB', 'X',\n'DHC', 'X',\n'DHD', '\
X',\n'DHE', 'X',\n'DHF', 'X',\n'DHG', 'X',\n'DHI',\
 'H',\n'DHL', 'X',\n'DHM', 'X',\n'DHN', 'V',\n'DHP\
', 'X',\n'DHQ', 'X',\n'DHR', 'X',\n'DHS', 'X',\n'D\
HT', 'X',\n'DHU', 'X',\n'DHY', 'X',\n'DHZ', 'X',\n\
'DI2', 'X',\n'DI3', 'G',\n'DI4', 'X',\n'DI5', 'X',\
\n'DIA', 'X',\n'DIC', 'X',\n'DIF', 'X',\n'DIG', 'X\
',\n'DII', 'X',\n'DIL', 'I',\n'DIM', 'X',\n'DIO', \
'X',\n'DIP', 'X',\n'DIQ', 'X',\n'DIS', 'X',\n'DIT'\
, 'X',\n'DIV', 'V',\n'DIX', 'X',\n'DIY', 'X',\n'DK\
A', 'X',\n'DLA', 'X',\n'DLE', 'L',\n'DLF', 'X',\n'\
DLS', 'K',\n'DLY', 'K',\n'DM1', 'X',\n'DM2', 'X',\\
n'DM3', 'X',\n'DM4', 'X',\n'DM5', 'X',\n'DM6', 'X'\
,\n'DM7', 'X',\n'DM8', 'X',\n'DM9', 'X',\n'DMA', '\
X',\n'DMB', 'X',\n'DMC', 'X',\n'DMD', 'X',\n'DME',\
 'X',\n'DMF', 'X',\n'DMG', 'G',\n'DMH', 'N',\n'DMI\
', 'X',\n'DMJ', 'X',\n'DML', 'X',\n'DMM', 'X',\n'D\
MN', 'X',\n'DMO', 'X',\n'DMP', 'X',\n'DMQ', 'X',\n\
'DMR', 'X',\n'DMS', 'X',\n'DMT', 'X',\n'DMV', 'X',\
\n'DMY', 'X',\n'DNC', 'X',\n'DND', 'X',\n'DNH', 'X\
',\n'DNJ', 'X',\n'DNN', 'X',\n'DNP', 'X',\n'DNQ', \
'X',\n'DNR', 'X',\n'DO2', 'X',\n'DO3', 'X',\n'DOA'\
, 'X',\n'DOB', 'X',\n'DOC', 'X',\n'DOH', 'D',\n'DO\
M', 'X',\n'DOS', 'X',\n'DOX', 'X',\n'DP5', 'X',\n'\
DP7', 'X',\n'DPA', 'X',\n'DPC', 'X',\n'DPD', 'X',\\
n'DPE', 'X',\n'DPG', 'X',\n'DPH', 'F',\n'DPM', 'X'\
,\n'DPN', 'F',\n'DPO', 'X',\n'DPP', 'X',\n'DPR', '\
P',\n'DPS', 'X',\n'DPT', 'X',\n'DPX', 'X',\n'DPY',\
 'X',\n'DPZ', 'X',\n'DQH', 'X',\n'DQN', 'X',\n'DR1\
', 'X',\n'DRB', 'X',\n'DRC', 'X',\n'DRI', 'X',\n'D\
RP', 'X',\n'DRT', 'X',\n'DRU', 'X',\n'DSA', 'X',\n\
'DSB', 'X',\n'DSC', 'X',\n'DSD', 'X',\n'DSE', 'S',\
\n'DSI', 'X',\n'DSN', 'S',\n'DSP', 'D',\n'DSR', 'X\
',\n'DSS', 'X',\n'DSX', 'X',\n'DSY', 'X',\n'DTB', \
'X',\n'DTD', 'X',\n'DTH', 'T',\n'DTN', 'X',\n'DTO'\
, 'X',\n'DTP', 'X',\n'DTQ', 'X',\n'DTR', 'W',\n'DT\
T', 'X',\n'DTY', 'Y',\n'DUD', 'X',\n'DUO', 'X',\n'\
DUR', 'X',\n'DUT', 'X',\n'DVA', 'V',\n'DVR', 'X',\\
n'DX9', 'X',\n'DXA', 'X',\n'DXB', 'X',\n'DXC', 'X'\
,\n'DXG', 'X',\n'DXX', 'X',\n'DZF', 'X',\n'E09', '\
X',\n'E20', 'X',\n'E2P', 'X',\n'E3G', 'X',\n'E4N',\
 'X',\n'E4P', 'X',\n'E64', 'X',\n'E6C', 'X',\n'E96\
', 'X',\n'E97', 'X',\n'EA2', 'X',\n'EAA', 'X',\n'E\
AP', 'X',\n'EBP', 'X',\n'EBW', 'X',\n'ECO', 'X',\n\
'EDA', 'X',\n'EDC', 'X',\n'EDE', 'X',\n'EDO', 'X',\
\n'EDR', 'X',\n'EEB', 'X',\n'EEE', 'X',\n'EFC', 'X\
',\n'EFZ', 'X',\n'EG1', 'X',\n'EG2', 'X',\n'EG3', \
'X',\n'EGC', 'X',\n'EGL', 'X',\n'EHP', 'A',\n'EIC'\
, 'X',\n'EJT', 'X',\n'ELA', 'X',\n'EMB', 'X',\n'EM\
C', 'X',\n'EMD', 'X',\n'EMM', 'X',\n'EMO', 'X',\n'\
EMP', 'X',\n'EMR', 'X',\n'ENA', 'X',\n'ENC', 'X',\\
n'ENH', 'X',\n'ENO', 'X',\n'ENP', 'X',\n'EOA', 'X'\
,\n'EOH', 'X',\n'EOT', 'X',\n'EOX', 'X',\n'EPA', '\
X',\n'EPE', 'X',\n'EPH', 'X',\n'EPI', 'X',\n'EPN',\
 'X',\n'EPO', 'X',\n'EPT', 'X',\n'EPU', 'X',\n'EPX\
', 'X',\n'EPY', 'X',\n'EQI', 'X',\n'EQP', 'X',\n'E\
QU', 'X',\n'ERG', 'X',\n'ERI', 'X',\n'ERY', 'X',\n\
'ESC', 'X',\n'ESD', 'X',\n'ESI', 'X',\n'ESO', 'X',\
\n'ESP', 'X',\n'EST', 'X',\n'ESX', 'X',\n'ETA', 'X\
',\n'ETC', 'X',\n'ETD', 'X',\n'ETF', 'X',\n'ETH', \
'X',\n'ETI', 'X',\n'ETN', 'X',\n'ETO', 'X',\n'ETP'\
, 'X',\n'ETR', 'X',\n'ETS', 'X',\n'ETY', 'X',\n'EU\
3', 'X',\n'EUG', 'X',\n'EYS', 'C',\n'F09', 'X',\n'\
F2B', 'X',\n'F3S', 'X',\n'F42', 'X',\n'F43', 'X',\\
n'F4S', 'X',\n'F6B', 'X',\n'F6P', 'X',\n'F89', 'X'\
,\n'FA1', 'X',\n'FA5', 'F',\n'FAA', 'X',\n'FAB', '\
X',\n'FAC', 'X',\n'FAD', 'X',\n'FAF', 'X',\n'FAG',\
 'X',\n'FAM', 'X',\n'FAR', 'X',\n'FAS', 'X',\n'FAT\
', 'X',\n'FBA', 'X',\n'FBE', 'X',\n'FBI', 'X',\n'F\
BP', 'X',\n'FBQ', 'X',\n'FBS', 'X',\n'FBT', 'X',\n\
'FBU', 'X',\n'FCA', 'X',\n'FCB', 'X',\n'FCI', 'X',\
\n'FCN', 'X',\n'FCO', 'X',\n'FCR', 'X',\n'FCT', 'X\
',\n'FCX', 'X',\n'FCY', 'C',\n'FD1', 'F',\n'FD2', \
'F',\n'FD3', 'F',\n'FD4', 'F',\n'FDA', 'X',\n'FDC'\
, 'X',\n'FDI', 'X',\n'FDP', 'X',\n'FDS', 'X',\n'FE\
2', 'X',\n'FEA', 'X',\n'FEL', 'X',\n'FEM', 'X',\n'\
FEN', 'X',\n'FEO', 'X',\n'FEP', 'X',\n'FER', 'X',\\
n'FES', 'X',\n'FFB', 'X',\n'FFC', 'X',\n'FFF', 'X'\
,\n'FFO', 'X',\n'FGL', 'G',\n'FHB', 'X',\n'FHC', '\
X',\n'FHP', 'X',\n'FHU', 'X',\n'FID', 'X',\n'FII',\
 'X',\n'FIP', 'X',\n'FK5', 'X',\n'FKA', 'X',\n'FKI\
', 'X',\n'FKP', 'X',\n'FL2', 'X',\n'FL9', 'X',\n'F\
LA', 'A',\n'FLC', 'X',\n'FLD', 'X',\n'FLE', 'L',\n\
'FLF', 'X',\n'FLO', 'X',\n'FLP', 'X',\n'FLT', 'Y',\
\n'FLU', 'X',\n'FLX', 'X',\n'FM1', 'X',\n'FM2', 'X\
',\n'FMA', 'X',\n'FMB', 'X',\n'FMC', 'X',\n'FME', \
'M',\n'FMN', 'X',\n'FMP', 'X',\n'FMR', 'X',\n'FMS'\
, 'X',\n'FMT', 'X',\n'FNE', 'X',\n'FNP', 'X',\n'FN\
S', 'X',\n'FOC', 'X',\n'FOE', 'X',\n'FOG', 'F',\n'\
FOH', 'X',\n'FOK', 'X',\n'FOL', 'X',\n'FON', 'X',\\
n'FOP', 'X',\n'FOR', 'X',\n'FOS', 'X',\n'FPA', 'X'\
,\n'FPC', 'X',\n'FPI', 'X',\n'FPO', 'X',\n'FPP', '\
X',\n'FPT', 'X',\n'FQP', 'X',\n'FRA', 'X',\n'FRD',\
 'F',\n'FRU', 'X',\n'FS3', 'X',\n'FS4', 'X',\n'FSB\
', 'X',\n'FSO', 'X',\n'FSX', 'X',\n'FTC', 'X',\n'F\
TP', 'X',\n'FTR', 'W',\n'FTT', 'X',\n'FTY', 'Y',\n\
'FUA', 'X',\n'FUC', 'X',\n'FUM', 'X',\n'FUP', 'X',\
\n'FVF', 'X',\n'FXP', 'X',\n'FXV', 'X',\n'FYA', 'F\
',\n'G16', 'X',\n'G1P', 'X',\n'G20', 'X',\n'G21', \
'X',\n'G23', 'X',\n'G26', 'X',\n'G28', 'X',\n'G2F'\
, 'X',\n'G37', 'X',\n'G39', 'X',\n'G3H', 'X',\n'G3\
P', 'X',\n'G4D', 'X',\n'G6D', 'X',\n'G6P', 'X',\n'\
G6Q', 'X',\n'G7M', 'X',\n'GA2', 'X',\n'GAA', 'X',\\
n'GAB', 'X',\n'GAC', 'X',\n'GAI', 'X',\n'GAL', 'X'\
,\n'GAM', 'X',\n'GAN', 'X',\n'GAO', 'X',\n'GAP', '\
X',\n'GAR', 'G',\n'GAS', 'X',\n'GAT', 'X',\n'GBC',\
 'X',\n'GBI', 'X',\n'GBP', 'X',\n'GBS', 'X',\n'GBX\
', 'X',\n'GC4', 'X',\n'GCA', 'X',\n'GCD', 'X',\n'G\
CG', 'G',\n'GCH', 'G',\n'GCK', 'X',\n'GCL', 'X',\n\
'GCM', 'X',\n'GCN', 'X',\n'GCO', 'X',\n'GCP', 'X',\
\n'GCR', 'X',\n'GCS', 'X',\n'GCU', 'X',\n'GD3', 'X\
',\n'GDB', 'X',\n'GDM', 'X',\n'GDN', 'X',\n'GDP', \
'X',\n'GDS', 'X',\n'GDU', 'X',\n'GE1', 'X',\n'GE2'\
, 'X',\n'GE3', 'X',\n'GEA', 'X',\n'GEL', 'X',\n'GE\
M', 'X',\n'GEN', 'X',\n'GEP', 'X',\n'GER', 'X',\n'\
GFP', 'X',\n'GGB', 'X',\n'GGL', 'E',\n'GGP', 'X',\\
n'GHP', 'G',\n'GIP', 'X',\n'GIS', 'X',\n'GKR', 'X'\
,\n'GL2', 'X',\n'GL3', 'G',\n'GL4', 'X',\n'GL5', '\
X',\n'GL7', 'X',\n'GL9', 'X',\n'GLA', 'X',\n'GLB',\
 'X',\n'GLC', 'X',\n'GLD', 'X',\n'GLE', 'X',\n'GLF\
', 'X',\n'GLG', 'X',\n'GLH', 'Q',\n'GLI', 'X',\n'G\
LL', 'X',\n'GLM', 'G',\n'GLN', 'Q',\n'GLO', 'X',\n\
'GLP', 'X',\n'GLR', 'X',\n'GLS', 'X',\n'GLT', 'X',\
\n'GLU', 'E',\n'GLV', 'X',\n'GLW', 'X',\n'GLY', 'G\
',\n'GLZ', 'X',\n'GM1', 'X',\n'GMA', 'X',\n'GMC', \
'X',\n'GMH', 'X',\n'GMP', 'X',\n'GMY', 'X',\n'GN7'\
, 'X',\n'GNA', 'X',\n'GNB', 'X',\n'GNH', 'X',\n'GN\
P', 'X',\n'GNT', 'X',\n'GOA', 'X',\n'GOL', 'X',\n'\
GOX', 'X',\n'GP1', 'X',\n'GP3', 'X',\n'GP4', 'X',\\
n'GP6', 'X',\n'GP8', 'X',\n'GPB', 'E',\n'GPC', 'X'\
,\n'GPE', 'X',\n'GPG', 'X',\n'GPI', 'X',\n'GPJ', '\
X',\n'GPL', 'K',\n'GPM', 'X',\n'GPN', 'G',\n'GPP',\
 'X',\n'GPR', 'X',\n'GPS', 'X',\n'GPX', 'X',\n'GR1\
', 'X',\n'GR3', 'X',\n'GR4', 'X',\n'GSA', 'X',\n'G\
SB', 'X',\n'GSC', 'G',\n'GSE', 'S',\n'GSH', 'X',\n\
'GSP', 'X',\n'GSR', 'X',\n'GSS', 'X',\n'GT9', 'C',\
\n'GTA', 'X',\n'GTB', 'X',\n'GTD', 'X',\n'GTE', 'X\
',\n'GTH', 'T',\n'GTN', 'X',\n'GTO', 'X',\n'GTP', \
'X',\n'GTR', 'X',\n'GTS', 'X',\n'GTT', 'X',\n'GTX'\
, 'X',\n'GTZ', 'X',\n'GU7', 'X',\n'GUA', 'X',\n'GU\
D', 'X',\n'GUM', 'X',\n'GUN', 'X',\n'GUP', 'X',\n'\
GUR', 'X',\n'GW3', 'X',\n'GZZ', 'X',\n'H2B', 'X',\\
n'H2P', 'H',\n'H2S', 'X',\n'H2U', 'X',\n'H4B', 'X'\
,\n'H5M', 'P',\n'H5P', 'X',\n'HAA', 'X',\n'HAB', '\
X',\n'HAC', 'A',\n'HAD', 'X',\n'HAE', 'X',\n'HAG',\
 'X',\n'HAI', 'X',\n'HAM', 'X',\n'HAP', 'X',\n'HAQ\
', 'X',\n'HAR', 'R',\n'HAS', 'X',\n'HAV', 'V',\n'H\
AX', 'X',\n'HAZ', 'X',\n'HBA', 'X',\n'HBC', 'X',\n\
'HBD', 'X',\n'HBI', 'X',\n'HBO', 'X',\n'HBU', 'X',\
\n'HBY', 'X',\n'HC0', 'X',\n'HC1', 'X',\n'HC4', 'X\
',\n'HCA', 'X',\n'HCC', 'X',\n'HCI', 'X',\n'HCS', \
'X',\n'HDA', 'X',\n'HDD', 'X',\n'HDF', 'X',\n'HDN'\
, 'X',\n'HDS', 'X',\n'HDZ', 'X',\n'HE1', 'X',\n'HE\
6', 'X',\n'HEA', 'X',\n'HEB', 'X',\n'HEC', 'X',\n'\
HED', 'X',\n'HEE', 'X',\n'HEF', 'X',\n'HEG', 'X',\\
n'HEM', 'X',\n'HEN', 'X',\n'HEO', 'X',\n'HEP', 'X'\
,\n'HEU', 'X',\n'HEV', 'X',\n'HEX', 'X',\n'HEZ', '\
X',\n'HF1', 'X',\n'HFA', 'X',\n'HFP', 'X',\n'HGA',\
 'Q',\n'HGB', 'X',\n'HGC', 'X',\n'HGI', 'X',\n'HGU\
', 'X',\n'HHO', 'X',\n'HHP', 'X',\n'HIB', 'X',\n'H\
IC', 'H',\n'HII', 'X',\n'HIN', 'X',\n'HIO', 'X',\n\
'HIP', 'H',\n'HIS', 'H',\n'HLE', 'X',\n'HLT', 'X',\
\n'HMA', 'A',\n'HMB', 'X',\n'HMC', 'X',\n'HMD', 'X\
',\n'HMF', 'A',\n'HMG', 'X',\n'HMH', 'X',\n'HMI', \
'L',\n'HMM', 'X',\n'HMN', 'X',\n'HMO', 'X',\n'HMP'\
, 'X',\n'HMR', 'R',\n'HNI', 'X',\n'HNP', 'X',\n'HO\
A', 'X',\n'HOE', 'X',\n'HOH', 'X',\n'HOM', 'X',\n'\
HOP', 'X',\n'HOQ', 'X',\n'HP1', 'A',\n'HP2', 'A',\\
n'HP3', 'X',\n'HPA', 'X',\n'HPB', 'X',\n'HPC', 'X'\
,\n'HPD', 'X',\n'HPE', 'A',\n'HPG', 'X',\n'HPH', '\
F',\n'HPP', 'X',\n'HPQ', 'F',\n'HPR', 'X',\n'HPT',\
 'X',\n'HPY', 'X',\n'HQO', 'X',\n'HQQ', 'X',\n'HQU\
', 'X',\n'HRG', 'R',\n'HRI', 'X',\n'HSA', 'X',\n'H\
SE', 'S',\n'HSF', 'X',\n'HSM', 'X',\n'HSO', 'H',\n\
'HSP', 'X',\n'HT1', 'X',\n'HT2', 'X',\n'HTA', 'X',\
\n'HTL', 'X',\n'HTO', 'X',\n'HTP', 'X',\n'HTR', 'W\
',\n'HUP', 'X',\n'HUX', 'X',\n'HV5', 'A',\n'HV7', \
'X',\n'HV8', 'X',\n'HXA', 'X',\n'HXC', 'X',\n'HXP'\
, 'X',\n'HY1', 'X',\n'HYA', 'X',\n'HYB', 'X',\n'HY\
D', 'X',\n'HYG', 'X',\n'HYP', 'P',\n'I06', 'X',\n'\
I10', 'X',\n'I11', 'X',\n'I17', 'X',\n'I2P', 'X',\\
n'I3N', 'X',\n'I3P', 'X',\n'I40', 'X',\n'I48', 'X'\
,\n'I4B', 'X',\n'I52', 'X',\n'I5P', 'X',\n'I84', '\
G',\n'IAG', 'G',\n'IAS', 'X',\n'IB2', 'X',\n'IBB',\
 'X',\n'IBP', 'X',\n'IBR', 'X',\n'IBS', 'X',\n'IBZ\
', 'X',\n'IC1', 'X',\n'ICA', 'X',\n'ICI', 'X',\n'I\
CL', 'X',\n'ICP', 'X',\n'ICT', 'X',\n'ICU', 'X',\n\
'ID2', 'X',\n'IDC', 'X',\n'IDG', 'X',\n'IDH', 'X',\
\n'IDM', 'X',\n'IDO', 'X',\n'IDP', 'X',\n'IDR', 'X\
',\n'IDS', 'X',\n'IDT', 'X',\n'IDU', 'X',\n'IFG', \
'X',\n'IFP', 'X',\n'IGL', 'X',\n'IGN', 'X',\n'IGP'\
, 'X',\n'IGU', 'X',\n'IH1', 'X',\n'IH2', 'X',\n'IH\
3', 'X',\n'IHB', 'X',\n'IHN', 'X',\n'IHP', 'X',\n'\
IIC', 'X',\n'IIL', 'I',\n'IIP', 'X',\n'IK2', 'X',\\
n'IKT', 'X',\n'ILA', 'I',\n'ILE', 'I',\n'ILG', 'X'\
,\n'ILO', 'X',\n'ILX', 'I',\n'IM1', 'X',\n'IM2', '\
X',\n'IMC', 'X',\n'IMD', 'X',\n'IME', 'X',\n'IMF',\
 'X',\n'IMG', 'X',\n'IMH', 'X',\n'IMI', 'X',\n'IML\
', 'I',\n'IMM', 'X',\n'IMN', 'X',\n'IMO', 'X',\n'I\
MP', 'X',\n'IMR', 'X',\n'IMU', 'X',\n'IN0', 'D',\n\
'IN1', 'R',\n'IN2', 'K',\n'IN3', 'L',\n'IN4', 'X',\
\n'IN5', 'A',\n'IN6', 'L',\n'IN7', 'X',\n'IN8', 'X\
',\n'IN9', 'X',\n'INA', 'L',\n'INB', 'X',\n'INC', \
'X',\n'IND', 'X',\n'INE', 'X',\n'INF', 'F',\n'ING'\
, 'F',\n'INH', 'R',\n'INI', 'X',\n'INJ', 'X',\n'IN\
K', 'X',\n'INL', 'X',\n'INM', 'X',\n'INN', 'A',\n'\
INO', 'X',\n'INP', 'X',\n'INQ', 'X',\n'INR', 'X',\\
n'INS', 'X',\n'INT', 'V',\n'INU', 'X',\n'INV', 'X'\
,\n'INW', 'X',\n'INX', 'X',\n'INY', 'X',\n'INZ', '\
X',\n'IOA', 'X',\n'IOB', 'X',\n'IOC', 'X',\n'IOD',\
 'X',\n'IOE', 'X',\n'IOF', 'X',\n'IOH', 'X',\n'IOL\
', 'X',\n'IOP', 'X',\n'IP1', 'X',\n'IP2', 'X',\n'I\
P3', 'X',\n'IP4', 'X',\n'IPA', 'X',\n'IPB', 'X',\n\
'IPD', 'X',\n'IPG', 'G',\n'IPH', 'X',\n'IPL', 'X',\
\n'IPM', 'X',\n'IPN', 'X',\n'IPO', 'F',\n'IPP', 'X\
',\n'IPS', 'X',\n'IPT', 'X',\n'IPU', 'X',\n'IPY', \
'A',\n'IQB', 'X',\n'IQP', 'X',\n'IQS', 'X',\n'IR3'\
, 'X',\n'IRI', 'X',\n'IRP', 'X',\n'ISA', 'X',\n'IS\
F', 'X',\n'ISO', 'X',\n'ISP', 'X',\n'ISQ', 'X',\n'\
ISU', 'X',\n'ITM', 'X',\n'ITP', 'X',\n'ITR', 'W',\\
n'ITS', 'X',\n'ITU', 'X',\n'IU5', 'X',\n'IUM', 'X'\
,\n'IUR', 'X',\n'IVA', 'X',\n'IYG', 'G',\n'IYR', '\
Y',\n'J77', 'X',\n'J78', 'X',\n'J80', 'X',\n'JE2',\
 'X',\n'JEN', 'X',\n'JST', 'X',\n'K21', 'X',\n'KAH\
', 'X',\n'KAI', 'X',\n'KAM', 'X',\n'KAN', 'X',\n'K\
AP', 'X',\n'KCP', 'X',\n'KCX', 'K',\n'KDO', 'X',\n\
'KEF', 'X',\n'KET', 'X',\n'KGR', 'X',\n'KH1', 'X',\
\n'KIF', 'X',\n'KIV', 'V',\n'KNI', 'X',\n'KPH', 'K\
',\n'KTH', 'X',\n'KTN', 'X',\n'KTP', 'X',\n'KWT', \
'X',\n'L04', 'X',\n'L1P', 'X',\n'L24', 'E',\n'L2P'\
, 'X',\n'L34', 'E',\n'L37', 'E',\n'L3P', 'X',\n'L4\
P', 'X',\n'L75', 'X',\n'LAC', 'X',\n'LAD', 'X',\n'\
LAK', 'X',\n'LAM', 'X',\n'LAR', 'X',\n'LAT', 'X',\\
n'LAX', 'X',\n'LCO', 'X',\n'LCP', 'X',\n'LCS', 'X'\
,\n'LDA', 'X',\n'LDO', 'L',\n'LDP', 'X',\n'LEA', '\
X',\n'LEO', 'X',\n'LEU', 'L',\n'LG2', 'X',\n'LG6',\
 'X',\n'LGC', 'X',\n'LGP', 'X',\n'LHG', 'X',\n'LHY\
', 'F',\n'LI1', 'X',\n'LIG', 'X',\n'LIL', 'X',\n'L\
IM', 'X',\n'LIN', 'X',\n'LIO', 'X',\n'LIP', 'X',\n\
'LLA', 'X',\n'LLP', 'K',\n'LLY', 'K',\n'LMG', 'X',\
\n'LML', 'X',\n'LMT', 'X',\n'LMU', 'X',\n'LMZ', 'X\
',\n'LNK', 'X',\n'LNL', 'X',\n'LNO', 'X',\n'LOF', \
'X',\n'LOL', 'L',\n'LOM', 'X',\n'LOR', 'X',\n'LOS'\
, 'X',\n'LOV', 'L',\n'LOX', 'X',\n'LP1', 'X',\n'LP\
2', 'R',\n'LPA', 'X',\n'LPC', 'X',\n'LPF', 'X',\n'\
LPL', 'X',\n'LPM', 'X',\n'LPP', 'X',\n'LRB', 'X',\\
n'LRU', 'X',\n'LS1', 'X',\n'LS2', 'X',\n'LS3', 'X'\
,\n'LS4', 'X',\n'LS5', 'X',\n'LTA', 'X',\n'LTL', '\
X',\n'LTR', 'W',\n'LUM', 'X',\n'LVS', 'L',\n'LXC',\
 'X',\n'LY2', 'X',\n'LY3', 'X',\n'LYA', 'X',\n'LYB\
', 'X',\n'LYC', 'X',\n'LYD', 'X',\n'LYM', 'K',\n'L\
YN', 'X',\n'LYS', 'K',\n'LYT', 'X',\n'LYW', 'X',\n\
'LYZ', 'K',\n'M1A', 'X',\n'M1G', 'X',\n'M2G', 'X',\
\n'M3L', 'K',\n'M6P', 'X',\n'M6T', 'X',\n'M7G', 'X\
',\n'MA1', 'X',\n'MA2', 'X',\n'MA3', 'X',\n'MA4', \
'X',\n'MA6', 'X',\n'MAA', 'A',\n'MAB', 'X',\n'MAC'\
, 'X',\n'MAE', 'X',\n'MAG', 'X',\n'MAH', 'X',\n'MA\
I', 'R',\n'MAK', 'X',\n'MAL', 'X',\n'MAM', 'X',\n'\
MAN', 'X',\n'MAO', 'X',\n'MAP', 'X',\n'MAR', 'X',\\
n'MAS', 'X',\n'MAT', 'X',\n'MAU', 'X',\n'MAZ', 'X'\
,\n'MBA', 'X',\n'MBD', 'X',\n'MBG', 'X',\n'MBH', '\
X',\n'MBN', 'X',\n'MBO', 'X',\n'MBR', 'X',\n'MBS',\
 'X',\n'MBV', 'X',\n'MBZ', 'X',\n'MCA', 'X',\n'MCD\
', 'X',\n'MCE', 'X',\n'MCG', 'G',\n'MCI', 'X',\n'M\
CN', 'X',\n'MCP', 'X',\n'MCT', 'X',\n'MCY', 'X',\n\
'MD2', 'X',\n'MDA', 'X',\n'MDC', 'X',\n'MDG', 'X',\
\n'MDH', 'X',\n'MDL', 'X',\n'MDM', 'X',\n'MDN', 'X\
',\n'MDP', 'X',\n'ME6', 'X',\n'MEB', 'X',\n'MEC', \
'X',\n'MEL', 'X',\n'MEN', 'N',\n'MEP', 'X',\n'MER'\
, 'X',\n'MES', 'X',\n'MET', 'M',\n'MEV', 'X',\n'MF\
2', 'X',\n'MF3', 'M',\n'MFB', 'X',\n'MFD', 'X',\n'\
MFU', 'X',\n'MG7', 'X',\n'MGA', 'X',\n'MGB', 'X',\\
n'MGD', 'X',\n'MGG', 'R',\n'MGL', 'X',\n'MGN', 'Q'\
,\n'MGO', 'X',\n'MGP', 'X',\n'MGR', 'X',\n'MGS', '\
X',\n'MGT', 'X',\n'MGU', 'X',\n'MGY', 'G',\n'MHB',\
 'X',\n'MHF', 'X',\n'MHL', 'L',\n'MHM', 'X',\n'MHO\
', 'M',\n'MHS', 'H',\n'MHZ', 'X',\n'MIA', 'X',\n'M\
IC', 'X',\n'MID', 'X',\n'MIL', 'X',\n'MIM', 'X',\n\
'MIN', 'G',\n'MIP', 'X',\n'MIS', 'S',\n'MIT', 'X',\
\n'MJI', 'X',\n'MK1', 'X',\n'MKC', 'X',\n'MLA', 'X\
',\n'MLC', 'X',\n'MLE', 'L',\n'MLN', 'X',\n'MLT', \
'X',\n'MLY', 'K',\n'MLZ', 'K',\n'MM3', 'X',\n'MM4'\
, 'X',\n'MMA', 'X',\n'MMC', 'X',\n'MME', 'M',\n'MM\
O', 'R',\n'MMP', 'X',\n'MMQ', 'X',\n'MMT', 'X',\n'\
MN1', 'X',\n'MN2', 'X',\n'MN3', 'X',\n'MN5', 'X',\\
n'MN7', 'X',\n'MN8', 'X',\n'MNA', 'X',\n'MNB', 'X'\
,\n'MNC', 'X',\n'MNG', 'X',\n'MNL', 'L',\n'MNO', '\
X',\n'MNP', 'X',\n'MNQ', 'X',\n'MNS', 'X',\n'MNT',\
 'X',\n'MNV', 'V',\n'MO1', 'X',\n'MO2', 'X',\n'MO3\
', 'X',\n'MO4', 'X',\n'MO5', 'X',\n'MO6', 'X',\n'M\
OA', 'X',\n'MOB', 'X',\n'MOC', 'X',\n'MOE', 'X',\n\
'MOG', 'X',\n'MOH', 'X',\n'MOL', 'X',\n'MOO', 'X',\
\n'MOP', 'X',\n'MOR', 'X',\n'MOS', 'X',\n'MOT', 'X\
',\n'MOX', 'X',\n'MP1', 'X',\n'MP3', 'X',\n'MPA', \
'X',\n'MPB', 'X',\n'MPC', 'X',\n'MPD', 'X',\n'MPG'\
, 'X',\n'MPH', 'M',\n'MPI', 'X',\n'MPJ', 'M',\n'MP\
L', 'X',\n'MPN', 'X',\n'MPO', 'X',\n'MPP', 'X',\n'\
MPQ', 'G',\n'MPR', 'X',\n'MPS', 'X',\n'MQ0', 'X',\\
n'MQ7', 'X',\n'MQ8', 'X',\n'MQ9', 'X',\n'MQI', 'X'\
,\n'MR2', 'X',\n'MRC', 'X',\n'MRM', 'X',\n'MRP', '\
X',\n'MS2', 'X',\n'MSA', 'X',\n'MSB', 'X',\n'MSD',\
 'X',\n'MSE', 'M',\n'MSF', 'X',\n'MSI', 'X',\n'MSO\
', 'M',\n'MSQ', 'X',\n'MST', 'X',\n'MSU', 'X',\n'M\
TA', 'X',\n'MTB', 'X',\n'MTC', 'X',\n'MTD', 'X',\n\
'MTE', 'X',\n'MTF', 'X',\n'MTG', 'X',\n'MTO', 'X',\
\n'MTS', 'X',\n'MTT', 'X',\n'MTX', 'X',\n'MTY', 'Y\
',\n'MUG', 'X',\n'MUP', 'X',\n'MUR', 'X',\n'MVA', \
'V',\n'MW1', 'X',\n'MW2', 'X',\n'MXA', 'X',\n'MXY'\
, 'X',\n'MYA', 'X',\n'MYC', 'X',\n'MYG', 'X',\n'MY\
R', 'X',\n'MYS', 'X',\n'MYT', 'X',\n'MZM', 'X',\n'\
N1T', 'X',\n'N25', 'X',\n'N2B', 'X',\n'N3T', 'X',\\
n'N4B', 'X',\n'NA2', 'X',\n'NA5', 'X',\n'NA6', 'X'\
,\n'NAA', 'X',\n'NAB', 'X',\n'NAC', 'X',\n'NAD', '\
X',\n'NAE', 'X',\n'NAF', 'X',\n'NAG', 'X',\n'NAH',\
 'X',\n'NAI', 'X',\n'NAL', 'A',\n'NAM', 'A',\n'NAN\
', 'X',\n'NAO', 'X',\n'NAP', 'X',\n'NAQ', 'X',\n'N\
AR', 'X',\n'NAS', 'X',\n'NAU', 'X',\n'NAV', 'X',\n\
'NAW', 'X',\n'NAX', 'X',\n'NAY', 'X',\n'NBA', 'X',\
\n'NBD', 'X',\n'NBE', 'X',\n'NBG', 'X',\n'NBN', 'X\
',\n'NBP', 'X',\n'NBS', 'X',\n'NBU', 'X',\n'NCA', \
'X',\n'NCB', 'A',\n'NCD', 'X',\n'NCH', 'X',\n'NCM'\
, 'X',\n'NCN', 'X',\n'NCO', 'X',\n'NCR', 'X',\n'NC\
S', 'X',\n'ND4', 'X',\n'NDA', 'X',\n'NDC', 'X',\n'\
NDD', 'X',\n'NDO', 'X',\n'NDP', 'X',\n'NDT', 'X',\\
n'NEA', 'X',\n'NEB', 'X',\n'NED', 'X',\n'NEM', 'H'\
,\n'NEN', 'X',\n'NEO', 'X',\n'NEP', 'H',\n'NEQ', '\
X',\n'NES', 'X',\n'NET', 'X',\n'NEV', 'X',\n'NFA',\
 'F',\n'NFE', 'X',\n'NFG', 'X',\n'NFP', 'X',\n'NFS\
', 'X',\n'NG6', 'X',\n'NGA', 'X',\n'NGL', 'X',\n'N\
GM', 'X',\n'NGO', 'X',\n'NGP', 'X',\n'NGT', 'X',\n\
'NGU', 'X',\n'NH2', 'X',\n'NH3', 'X',\n'NH4', 'X',\
\n'NHD', 'X',\n'NHE', 'X',\n'NHM', 'X',\n'NHP', 'X\
',\n'NHR', 'X',\n'NHS', 'X',\n'NI1', 'X',\n'NI2', \
'X',\n'NIC', 'X',\n'NID', 'X',\n'NIK', 'X',\n'NIO'\
, 'X',\n'NIP', 'X',\n'NIT', 'X',\n'NIU', 'X',\n'NI\
Y', 'Y',\n'NLA', 'X',\n'NLE', 'L',\n'NLG', 'X',\n'\
NLN', 'L',\n'NLP', 'L',\n'NM1', 'X',\n'NMA', 'A',\\
n'NMB', 'X',\n'NMC', 'G',\n'NMD', 'X',\n'NME', 'X'\
,\n'NMN', 'X',\n'NMO', 'X',\n'NMQ', 'X',\n'NMX', '\
X',\n'NMY', 'X',\n'NNH', 'R',\n'NNO', 'X',\n'NO2',\
 'X',\n'NO3', 'X',\n'NOA', 'X',\n'NOD', 'X',\n'NOJ\
', 'X',\n'NON', 'X',\n'NOP', 'X',\n'NOR', 'X',\n'N\
OS', 'X',\n'NOV', 'X',\n'NOX', 'X',\n'NP3', 'X',\n\
'NPA', 'X',\n'NPC', 'X',\n'NPD', 'X',\n'NPE', 'X',\
\n'NPF', 'X',\n'NPH', 'C',\n'NPI', 'X',\n'NPL', 'X\
',\n'NPN', 'X',\n'NPO', 'X',\n'NPP', 'X',\n'NPT', \
'X',\n'NPY', 'X',\n'NRG', 'R',\n'NRI', 'X',\n'NS1'\
, 'X',\n'NS5', 'X',\n'NSP', 'X',\n'NTA', 'X',\n'NT\
B', 'X',\n'NTC', 'X',\n'NTH', 'X',\n'NTM', 'X',\n'\
NTP', 'X',\n'NTS', 'X',\n'NTU', 'X',\n'NTZ', 'X',\\
n'NU1', 'X',\n'NVA', 'V',\n'NVI', 'X',\n'NVP', 'X'\
,\n'NW1', 'X',\n'NYP', 'X',\n'O4M', 'X',\n'OAA', '\
X',\n'OAI', 'X',\n'OAP', 'X',\n'OAR', 'X',\n'OAS',\
 'S',\n'OBA', 'X',\n'OBN', 'X',\n'OC1', 'X',\n'OC2\
', 'X',\n'OC3', 'X',\n'OC4', 'X',\n'OC5', 'X',\n'O\
C6', 'X',\n'OC7', 'X',\n'OCL', 'X',\n'OCM', 'X',\n\
'OCN', 'X',\n'OCO', 'X',\n'OCP', 'X',\n'OCS', 'C',\
\n'OCT', 'X',\n'OCV', 'K',\n'OCY', 'C',\n'ODA', 'X\
',\n'ODS', 'X',\n'OES', 'X',\n'OET', 'X',\n'OF1', \
'X',\n'OF2', 'X',\n'OF3', 'X',\n'OFL', 'X',\n'OFO'\
, 'X',\n'OHE', 'X',\n'OHO', 'X',\n'OHT', 'X',\n'OI\
C', 'X',\n'OIP', 'X',\n'OKA', 'X',\n'OLA', 'X',\n'\
OLE', 'X',\n'OLI', 'X',\n'OLO', 'X',\n'OMB', 'X',\\
n'OMC', 'X',\n'OMD', 'X',\n'OME', 'X',\n'OMG', 'X'\
,\n'OMP', 'X',\n'OMT', 'M',\n'OMU', 'X',\n'ONE', '\
X',\n'ONL', 'L',\n'ONP', 'X',\n'OPA', 'X',\n'OPD',\
 'X',\n'OPE', 'X',\n'OPG', 'X',\n'OPH', 'X',\n'OPN\
', 'X',\n'OPP', 'X',\n'OPR', 'R',\n'ORN', 'X',\n'O\
RO', 'X',\n'ORP', 'X',\n'OSB', 'X',\n'OSS', 'X',\n\
'OTA', 'X',\n'OTB', 'X',\n'OTE', 'X',\n'OTG', 'X',\
\n'OUT', 'X',\n'OVA', 'X',\n'OWQ', 'X',\n'OXA', 'X\
',\n'OXE', 'X',\n'OXI', 'X',\n'OXL', 'X',\n'OXM', \
'X',\n'OXN', 'X',\n'OXO', 'X',\n'OXP', 'X',\n'OXS'\
, 'X',\n'OXY', 'X',\n'P11', 'A',\n'P24', 'X',\n'P2\
8', 'X',\n'P2P', 'X',\n'P2U', 'X',\n'P3M', 'X',\n'\
P4C', 'X',\n'P4P', 'X',\n'P5P', 'X',\n'P6G', 'X',\\
n'PA1', 'X',\n'PA2', 'X',\n'PA3', 'X',\n'PA4', 'X'\
,\n'PA5', 'X',\n'PAA', 'X',\n'PAB', 'X',\n'PAC', '\
X',\n'PAD', 'X',\n'PAE', 'X',\n'PAG', 'X',\n'PAH',\
 'X',\n'PAI', 'X',\n'PAL', 'D',\n'PAM', 'X',\n'PAN\
', 'X',\n'PAO', 'X',\n'PAP', 'A',\n'PAQ', 'F',\n'P\
AR', 'X',\n'PAS', 'X',\n'PAT', 'W',\n'PBA', 'X',\n\
'PBB', 'X',\n'PBC', 'X',\n'PBF', 'F',\n'PBG', 'X',\
\n'PBI', 'X',\n'PBM', 'X',\n'PBN', 'X',\n'PBP', 'X\
',\n'PBR', 'X',\n'PBZ', 'X',\n'PC2', 'X',\n'PCA', \
'E',\n'PCB', 'X',\n'PCD', 'X',\n'PCE', 'X',\n'PCG'\
, 'X',\n'PCH', 'X',\n'PCL', 'X',\n'PCM', 'X',\n'PC\
P', 'X',\n'PCR', 'X',\n'PCS', 'X',\n'PCU', 'X',\n'\
PCV', 'X',\n'PCY', 'X',\n'PD1', 'X',\n'PDA', 'X',\\
n'PDC', 'X',\n'PDD', 'A',\n'PDE', 'A',\n'PDI', 'X'\
,\n'PDL', 'A',\n'PDN', 'X',\n'PDO', 'X',\n'PDP', '\
X',\n'PDT', 'X',\n'PDU', 'X',\n'PE2', 'X',\n'PE6',\
 'X',\n'PEA', 'X',\n'PEB', 'X',\n'PEC', 'X',\n'PED\
', 'X',\n'PEE', 'X',\n'PEF', 'X',\n'PEG', 'X',\n'P\
EL', 'X',\n'PEO', 'X',\n'PEP', 'X',\n'PEQ', 'X',\n\
'PER', 'X',\n'PET', 'X',\n'PFB', 'X',\n'PFC', 'X',\
\n'PFG', 'X',\n'PFL', 'X',\n'PFM', 'X',\n'PFZ', 'X\
',\n'PG4', 'X',\n'PG5', 'X',\n'PG6', 'X',\n'PGA', \
'X',\n'PGC', 'X',\n'PGD', 'X',\n'PGE', 'X',\n'PGG'\
, 'G',\n'PGH', 'X',\n'PGL', 'X',\n'PGO', 'X',\n'PG\
P', 'X',\n'PGQ', 'X',\n'PGR', 'X',\n'PGS', 'X',\n'\
PGU', 'X',\n'PGX', 'X',\n'PGY', 'G',\n'PH1', 'X',\\
n'PH2', 'X',\n'PH3', 'X',\n'PHA', 'F',\n'PHB', 'X'\
,\n'PHC', 'X',\n'PHD', 'X',\n'PHE', 'F',\n'PHG', '\
X',\n'PHH', 'X',\n'PHI', 'F',\n'PHL', 'F',\n'PHM',\
 'X',\n'PHN', 'X',\n'PHO', 'X',\n'PHP', 'X',\n'PHQ\
', 'X',\n'PHS', 'H',\n'PHT', 'X',\n'PHW', 'P',\n'P\
HY', 'X',\n'PI1', 'X',\n'PI2', 'X',\n'PI3', 'X',\n\
'PI4', 'X',\n'PI5', 'X',\n'PI6', 'X',\n'PI7', 'X',\
\n'PI8', 'X',\n'PI9', 'X',\n'PIA', 'X',\n'PIB', 'X\
',\n'PIC', 'X',\n'PID', 'X',\n'PIG', 'X',\n'PIH', \
'X',\n'PIM', 'X',\n'PIN', 'X',\n'PIO', 'X',\n'PIP'\
, 'X',\n'PIQ', 'X',\n'PIR', 'X',\n'PIV', 'X',\n'PK\
F', 'X',\n'PL1', 'X',\n'PL9', 'X',\n'PLA', 'D',\n'\
PLC', 'X',\n'PLE', 'L',\n'PLG', 'G',\n'PLH', 'X',\\
n'PLM', 'X',\n'PLP', 'X',\n'PLS', 'S',\n'PLT', 'W'\
,\n'PLU', 'L',\n'PLY', 'X',\n'PMA', 'X',\n'PMB', '\
X',\n'PMC', 'X',\n'PME', 'F',\n'PML', 'X',\n'PMM',\
 'X',\n'PMO', 'X',\n'PMP', 'X',\n'PMS', 'X',\n'PMY\
', 'X',\n'PN2', 'X',\n'PNA', 'X',\n'PNB', 'X',\n'P\
NC', 'G',\n'PND', 'X',\n'PNE', 'A',\n'PNF', 'X',\n\
'PNG', 'X',\n'PNI', 'X',\n'PNL', 'X',\n'PNM', 'X',\
\n'PNN', 'X',\n'PNO', 'X',\n'PNP', 'X',\n'PNQ', 'X\
',\n'PNS', 'X',\n'PNT', 'X',\n'PNU', 'X',\n'PO2', \
'X',\n'PO4', 'X',\n'POB', 'X',\n'POC', 'X',\n'POL'\
, 'X',\n'POM', 'P',\n'PON', 'X',\n'POP', 'X',\n'PO\
R', 'X',\n'POS', 'X',\n'PP1', 'X',\n'PP2', 'X',\n'\
PP3', 'A',\n'PP4', 'X',\n'PP5', 'X',\n'PP6', 'X',\\
n'PP7', 'X',\n'PP8', 'N',\n'PP9', 'X',\n'PPB', 'X'\
,\n'PPC', 'X',\n'PPD', 'X',\n'PPE', 'E',\n'PPG', '\
X',\n'PPH', 'F',\n'PPI', 'X',\n'PPJ', 'V',\n'PPL',\
 'X',\n'PPM', 'X',\n'PPN', 'A',\n'PPO', 'X',\n'PPP\
', 'X',\n'PPQ', 'X',\n'PPR', 'X',\n'PPS', 'X',\n'P\
PT', 'X',\n'PPU', 'X',\n'PPX', 'F',\n'PPY', 'X',\n\
'PPZ', 'X',\n'PQ0', 'X',\n'PQN', 'X',\n'PQQ', 'X',\
\n'PR1', 'X',\n'PR2', 'X',\n'PR3', 'X',\n'PRA', 'X\
',\n'PRB', 'X',\n'PRC', 'X',\n'PRD', 'X',\n'PRE', \
'X',\n'PRF', 'X',\n'PRH', 'X',\n'PRI', 'P',\n'PRL'\
, 'X',\n'PRN', 'X',\n'PRO', 'P',\n'PRP', 'X',\n'PR\
R', 'A',\n'PRS', 'P',\n'PRZ', 'X',\n'PS0', 'X',\n'\
PSA', 'X',\n'PSD', 'X',\n'PSE', 'X',\n'PSF', 'S',\\
n'PSG', 'X',\n'PSI', 'X',\n'PSO', 'X',\n'PSQ', 'X'\
,\n'PSS', 'X',\n'PST', 'X',\n'PSU', 'X',\n'PT1', '\
X',\n'PT3', 'X',\n'PTA', 'X',\n'PTC', 'X',\n'PTD',\
 'X',\n'PTE', 'X',\n'PTH', 'Y',\n'PTL', 'X',\n'PTM\
', 'Y',\n'PTN', 'X',\n'PTO', 'X',\n'PTP', 'X',\n'P\
TR', 'Y',\n'PTS', 'X',\n'PTT', 'X',\n'PTU', 'X',\n\
'PTY', 'X',\n'PUA', 'X',\n'PUB', 'X',\n'PUR', 'X',\
\n'PUT', 'X',\n'PVA', 'X',\n'PVB', 'X',\n'PVH', 'H\
',\n'PVL', 'X',\n'PXA', 'X',\n'PXF', 'X',\n'PXG', \
'X',\n'PXP', 'X',\n'PXY', 'X',\n'PXZ', 'X',\n'PY2'\
, 'X',\n'PY4', 'X',\n'PY5', 'X',\n'PY6', 'X',\n'PY\
A', 'A',\n'PYC', 'X',\n'PYD', 'X',\n'PYE', 'X',\n'\
PYL', 'X',\n'PYM', 'X',\n'PYO', 'X',\n'PYP', 'X',\\
n'PYQ', 'X',\n'PYR', 'X',\n'PYS', 'X',\n'PYT', 'X'\
,\n'PYX', 'X',\n'PYY', 'X',\n'PYZ', 'X',\n'PZQ', '\
X',\n'Q82', 'X',\n'QNC', 'X',\n'QND', 'X',\n'QSI',\
 'Q',\n'QTR', 'X',\n'QUA', 'X',\n'QUE', 'X',\n'QUI\
', 'X',\n'QUO', 'X',\n'R11', 'X',\n'R12', 'X',\n'R\
13', 'X',\n'R18', 'X',\n'R1P', 'X',\n'R56', 'X',\n\
'R5P', 'X',\n'RA2', 'X',\n'RAD', 'X',\n'RAI', 'X',\
\n'RAL', 'X',\n'RAM', 'X',\n'RAN', 'X',\n'RAP', 'X\
',\n'RBF', 'X',\n'RBU', 'X',\n'RCA', 'X',\n'RCL', \
'X',\n'RCO', 'X',\n'RDC', 'X',\n'RDF', 'W',\n'RE9'\
, 'X',\n'REA', 'X',\n'RED', 'K',\n'REO', 'X',\n'RE\
P', 'X',\n'RET', 'X',\n'RFA', 'X',\n'RFB', 'X',\n'\
RFL', 'X',\n'RFP', 'X',\n'RG1', 'X',\n'RGS', 'X',\\
n'RH1', 'X',\n'RHA', 'X',\n'RHC', 'X',\n'RHD', 'X'\
,\n'RHM', 'X',\n'RHO', 'X',\n'RHQ', 'X',\n'RHS', '\
X',\n'RIA', 'X',\n'RIB', 'X',\n'RIC', 'X',\n'RIF',\
 'X',\n'RIN', 'X',\n'RIP', 'X',\n'RIT', 'X',\n'RMB\
', 'X',\n'RMN', 'X',\n'RMP', 'X',\n'RNG', 'X',\n'R\
NS', 'X',\n'RNT', 'X',\n'RO2', 'X',\n'RO4', 'X',\n\
'ROC', 'N',\n'ROI', 'X',\n'ROM', 'X',\n'RON', 'V',\
\n'ROP', 'X',\n'ROS', 'X',\n'ROX', 'X',\n'RPA', 'X\
',\n'RPD', 'X',\n'RPH', 'X',\n'RPL', 'X',\n'RPP', \
'X',\n'RPR', 'X',\n'RPX', 'X',\n'RQ3', 'X',\n'RR1'\
, 'X',\n'RR6', 'X',\n'RRS', 'X',\n'RS1', 'X',\n'RS\
2', 'X',\n'RS7', 'X',\n'RSS', 'X',\n'RTA', 'X',\n'\
RTB', 'X',\n'RTC', 'X',\n'RTL', 'X',\n'RUB', 'X',\\
n'RUN', 'X',\n'RWJ', 'X',\n'RXP', 'X',\n'S02', 'X'\
,\n'S11', 'X',\n'S1H', 'S',\n'S27', 'X',\n'S2C', '\
C',\n'S3P', 'X',\n'S4U', 'X',\n'S57', 'X',\n'S58',\
 'X',\n'S5H', 'X',\n'S6G', 'X',\n'S80', 'X',\n'SAA\
', 'X',\n'SAB', 'X',\n'SAC', 'S',\n'SAD', 'X',\n'S\
AE', 'X',\n'SAF', 'X',\n'SAH', 'C',\n'SAI', 'C',\n\
'SAL', 'X',\n'SAM', 'M',\n'SAN', 'X',\n'SAP', 'X',\
\n'SAR', 'X',\n'SAS', 'X',\n'SB1', 'X',\n'SB2', 'X\
',\n'SB3', 'X',\n'SB4', 'X',\n'SB5', 'X',\n'SB6', \
'X',\n'SBA', 'L',\n'SBB', 'X',\n'SBD', 'A',\n'SBI'\
, 'X',\n'SBL', 'A',\n'SBN', 'X',\n'SBO', 'X',\n'SB\
R', 'X',\n'SBS', 'X',\n'SBT', 'X',\n'SBU', 'X',\n'\
SBX', 'X',\n'SC4', 'X',\n'SCA', 'X',\n'SCC', 'X',\\
n'SCD', 'X',\n'SCH', 'C',\n'SCI', 'X',\n'SCL', 'X'\
,\n'SCM', 'X',\n'SCN', 'X',\n'SCO', 'X',\n'SCP', '\
S',\n'SCR', 'X',\n'SCS', 'X',\n'SCV', 'C',\n'SCY',\
 'C',\n'SD8', 'X',\n'SDK', 'X',\n'SDZ', 'X',\n'SE4\
', 'X',\n'SEA', 'X',\n'SEB', 'S',\n'SEC', 'X',\n'S\
EG', 'A',\n'SEI', 'X',\n'SEL', 'S',\n'SEM', 'X',\n\
'SEO', 'X',\n'SEP', 'S',\n'SER', 'S',\n'SES', 'X',\
\n'SET', 'S',\n'SEU', 'X',\n'SF4', 'X',\n'SFG', 'X\
',\n'SFN', 'X',\n'SFO', 'X',\n'SGA', 'X',\n'SGC', \
'X',\n'SGL', 'X',\n'SGM', 'X',\n'SGN', 'X',\n'SGP'\
, 'X',\n'SHA', 'X',\n'SHC', 'X',\n'SHF', 'X',\n'SH\
H', 'X',\n'SHP', 'G',\n'SHR', 'E',\n'SHT', 'T',\n'\
SHU', 'X',\n'SI2', 'X',\n'SIA', 'X',\n'SIF', 'X',\\
n'SIG', 'X',\n'SIH', 'X',\n'SIM', 'X',\n'SIN', 'X'\
,\n'SKD', 'X',\n'SKF', 'X',\n'SLB', 'X',\n'SLE', '\
X',\n'SLZ', 'K',\n'SMA', 'X',\n'SMC', 'C',\n'SME',\
 'M',\n'SML', 'X',\n'SMM', 'M',\n'SMN', 'X',\n'SMP\
', 'X',\n'SMS', 'X',\n'SN1', 'X',\n'SN6', 'X',\n'S\
N7', 'X',\n'SNC', 'C',\n'SNN', 'X',\n'SNP', 'X',\n\
'SO1', 'X',\n'SO2', 'X',\n'SO3', 'X',\n'SO4', 'X',\
\n'SOA', 'X',\n'SOC', 'C',\n'SOM', 'X',\n'SOR', 'X\
',\n'SOT', 'X',\n'SOX', 'X',\n'SPA', 'X',\n'SPB', \
'X',\n'SPC', 'X',\n'SPD', 'X',\n'SPE', 'X',\n'SPG'\
, 'X',\n'SPH', 'X',\n'SPI', 'X',\n'SPK', 'X',\n'SP\
M', 'X',\n'SPN', 'X',\n'SPO', 'X',\n'SPP', 'X',\n'\
SPS', 'X',\n'SPY', 'X',\n'SQU', 'X',\n'SRA', 'X',\\
n'SRB', 'X',\n'SRD', 'X',\n'SRL', 'X',\n'SRM', 'X'\
,\n'SRS', 'X',\n'SRY', 'X',\n'SSA', 'X',\n'SSB', '\
X',\n'SSG', 'X',\n'SSP', 'X',\n'ST1', 'X',\n'ST2',\
 'X',\n'ST3', 'X',\n'ST4', 'X',\n'ST5', 'X',\n'ST6\
', 'X',\n'STA', 'X',\n'STB', 'X',\n'STE', 'X',\n'S\
TG', 'X',\n'STI', 'X',\n'STL', 'X',\n'STN', 'X',\n\
'STO', 'X',\n'STP', 'X',\n'STR', 'X',\n'STU', 'X',\
\n'STY', 'Y',\n'SU1', 'X',\n'SU2', 'X',\n'SUC', 'X\
',\n'SUI', 'X',\n'SUL', 'X',\n'SUR', 'X',\n'SVA', \
'S',\n'SWA', 'X',\n'T16', 'X',\n'T19', 'X',\n'T23'\
, 'X',\n'T29', 'X',\n'T33', 'X',\n'T3P', 'X',\n'T4\
2', 'A',\n'T44', 'X',\n'T5A', 'X',\n'T6A', 'T',\n'\
T6P', 'X',\n'T80', 'X',\n'T87', 'X',\n'TA1', 'X',\\
n'TAA', 'X',\n'TAB', 'X',\n'TAC', 'X',\n'TAD', 'X'\
,\n'TAF', 'X',\n'TAM', 'X',\n'TAP', 'X',\n'TAR', '\
X',\n'TAS', 'X',\n'TAU', 'X',\n'TAX', 'X',\n'TAZ',\
 'X',\n'TB9', 'X',\n'TBA', 'X',\n'TBD', 'X',\n'TBG\
', 'G',\n'TBH', 'X',\n'TBM', 'T',\n'TBO', 'X',\n'T\
BP', 'X',\n'TBR', 'X',\n'TBS', 'X',\n'TBT', 'X',\n\
'TBU', 'X',\n'TBZ', 'X',\n'TC4', 'X',\n'TCA', 'X',\
\n'TCB', 'X',\n'TCH', 'X',\n'TCK', 'X',\n'TCL', 'X\
',\n'TCM', 'X',\n'TCN', 'X',\n'TCP', 'X',\n'TCR', \
'W',\n'TCS', 'X',\n'TCZ', 'X',\n'TDA', 'X',\n'TDB'\
, 'X',\n'TDG', 'X',\n'TDP', 'X',\n'TDR', 'X',\n'TD\
X', 'X',\n'TEA', 'X',\n'TEM', 'X',\n'TEN', 'X',\n'\
TEO', 'X',\n'TEP', 'X',\n'TER', 'X',\n'TES', 'X',\\
n'TET', 'X',\n'TFA', 'X',\n'TFB', 'X',\n'TFH', 'X'\
,\n'TFI', 'X',\n'TFK', 'X',\n'TFP', 'X',\n'THA', '\
X',\n'THB', 'X',\n'THC', 'T',\n'THD', 'X',\n'THE',\
 'X',\n'THF', 'X',\n'THJ', 'X',\n'THK', 'X',\n'THM\
', 'X',\n'THN', 'X',\n'THO', 'T',\n'THP', 'X',\n'T\
HQ', 'X',\n'THR', 'T',\n'THS', 'X',\n'THT', 'X',\n\
'THU', 'X',\n'THX', 'X',\n'THZ', 'X',\n'TI1', 'X',\
\n'TI2', 'X',\n'TI3', 'P',\n'TIA', 'X',\n'TIH', 'A\
',\n'TK4', 'X',\n'TLA', 'X',\n'TLC', 'X',\n'TLM', \
'X',\n'TLN', 'X',\n'TLX', 'X',\n'TM5', 'X',\n'TM6'\
, 'X',\n'TMA', 'X',\n'TMB', 'T',\n'TMC', 'X',\n'TM\
D', 'T',\n'TME', 'X',\n'TMF', 'X',\n'TML', 'K',\n'\
TMM', 'X',\n'TMN', 'X',\n'TMP', 'X',\n'TMQ', 'X',\\
n'TMR', 'X',\n'TMT', 'X',\n'TMZ', 'X',\n'TNB', 'C'\
,\n'TND', 'X',\n'TNK', 'X',\n'TNP', 'X',\n'TNT', '\
X',\n'TOA', 'X',\n'TOB', 'X',\n'TOC', 'X',\n'TOL',\
 'X',\n'TOP', 'X',\n'TOS', 'X',\n'TOT', 'X',\n'TP1\
', 'G',\n'TP2', 'P',\n'TP3', 'E',\n'TP4', 'E',\n'T\
P7', 'T',\n'TPA', 'X',\n'TPE', 'X',\n'TPF', 'X',\n\
'TPI', 'X',\n'TPL', 'W',\n'TPM', 'X',\n'TPN', 'G',\
\n'TPO', 'T',\n'TPP', 'X',\n'TPQ', 'A',\n'TPR', 'P\
',\n'TPS', 'X',\n'TPT', 'X',\n'TPV', 'X',\n'TPX', \
'X',\n'TPY', 'X',\n'TQ3', 'X',\n'TQ4', 'X',\n'TQ5'\
, 'X',\n'TQ6', 'X',\n'TR1', 'X',\n'TRA', 'X',\n'TR\
B', 'X',\n'TRC', 'X',\n'TRD', 'X',\n'TRE', 'X',\n'\
TRF', 'W',\n'TRG', 'K',\n'TRH', 'X',\n'TRI', 'X',\\
n'TRJ', 'X',\n'TRM', 'X',\n'TRN', 'W',\n'TRO', 'W'\
,\n'TRP', 'W',\n'TRQ', 'X',\n'TRS', 'X',\n'TRX', '\
W',\n'TRZ', 'X',\n'TS2', 'X',\n'TS3', 'X',\n'TS4',\
 'X',\n'TS5', 'X',\n'TSA', 'X',\n'TSB', 'X',\n'TSI\
', 'X',\n'TSM', 'X',\n'TSN', 'X',\n'TSP', 'X',\n'T\
SU', 'X',\n'TTA', 'X',\n'TTE', 'X',\n'TTN', 'X',\n\
'TTO', 'X',\n'TTP', 'X',\n'TTX', 'X',\n'TXL', 'X',\
\n'TYA', 'Y',\n'TYB', 'Y',\n'TYD', 'X',\n'TYI', 'Y\
',\n'TYL', 'X',\n'TYM', 'W',\n'TYN', 'Y',\n'TYQ', \
'Y',\n'TYR', 'Y',\n'TYS', 'Y',\n'TYV', 'X',\n'TYY'\
, 'A',\n'TZB', 'X',\n'TZC', 'X',\n'TZE', 'X',\n'TZ\
L', 'X',\n'TZO', 'X',\n'TZP', 'X',\n'U01', 'X',\n'\
U02', 'X',\n'U03', 'X',\n'U04', 'X',\n'U05', 'X',\\
n'U0E', 'X',\n'U10', 'X',\n'U18', 'X',\n'U2G', 'X'\
,\n'U3P', 'X',\n'U49', 'X',\n'U55', 'X',\n'U5P', '\
X',\n'U66', 'X',\n'U89', 'X',\n'U8U', 'X',\n'UAA',\
 'X',\n'UAG', 'A',\n'UAP', 'X',\n'UAR', 'X',\n'UC1\
', 'X',\n'UC2', 'X',\n'UC3', 'X',\n'UC4', 'X',\n'U\
D1', 'X',\n'UD2', 'X',\n'UDP', 'X',\n'UDX', 'X',\n\
'UFG', 'X',\n'UFM', 'X',\n'UFP', 'X',\n'UGA', 'X',\
\n'UIN', 'X',\n'UKP', 'A',\n'UM3', 'X',\n'UMA', 'A\
',\n'UMG', 'X',\n'UMP', 'X',\n'UNA', 'X',\n'UND', \
'X',\n'UNI', 'X',\n'UNK', 'X',\n'UNN', 'X',\n'UNX'\
, 'X',\n'UP5', 'X',\n'UP6', 'X',\n'UPA', 'X',\n'UP\
F', 'X',\n'UPG', 'X',\n'UPP', 'X',\n'UQ1', 'X',\n'\
UQ2', 'X',\n'UQ6', 'X',\n'UR2', 'X',\n'URA', 'X',\\
n'URE', 'X',\n'URF', 'X',\n'URI', 'X',\n'URS', 'X'\
,\n'UTP', 'X',\n'UVC', 'X',\n'UVW', 'X',\n'V35', '\
X',\n'V36', 'X',\n'V4O', 'X',\n'V7O', 'X',\n'VAA',\
 'V',\n'VAC', 'X',\n'VAD', 'V',\n'VAF', 'V',\n'VAG\
', 'X',\n'VAL', 'V',\n'VAN', 'X',\n'VAS', 'X',\n'V\
AX', 'X',\n'VDX', 'X',\n'VDY', 'X',\n'VG1', 'X',\n\
'VIB', 'X',\n'VIR', 'X',\n'VIT', 'X',\n'VK3', 'X',\
\n'VO3', 'X',\n'VO4', 'X',\n'VS1', 'F',\n'VS2', 'F\
',\n'VS3', 'F',\n'VS4', 'F',\n'VXA', 'X',\n'W01', \
'X',\n'W02', 'X',\n'W03', 'X',\n'W11', 'X',\n'W33'\
, 'X',\n'W35', 'X',\n'W42', 'X',\n'W43', 'X',\n'W5\
4', 'X',\n'W56', 'X',\n'W59', 'X',\n'W71', 'X',\n'\
W84', 'X',\n'W8R', 'X',\n'W91', 'X',\n'WAY', 'X',\\
n'WCC', 'X',\n'WO2', 'X',\n'WO4', 'X',\n'WRB', 'X'\
,\n'WRR', 'X',\n'WRS', 'X',\n'WW7', 'X',\n'X2F', '\
X',\n'X7O', 'X',\n'XAA', 'X',\n'XAN', 'X',\n'XAO',\
 'X',\n'XBB', 'X',\n'XBP', 'X',\n'XDN', 'X',\n'XDP\
', 'X',\n'XIF', 'X',\n'XIM', 'X',\n'XK2', 'X',\n'X\
L1', 'X',\n'XLS', 'X',\n'XMP', 'X',\n'XN1', 'X',\n\
'XN2', 'X',\n'XN3', 'X',\n'XUL', 'X',\n'XV6', 'X',\
\n'XYD', 'X',\n'XYH', 'X',\n'XYL', 'X',\n'XYP', 'X\
',\n'XYS', 'X',\n'YOF', 'Y',\n'YRR', 'X',\n'YT3', \
'X',\n'YZ9', 'X',\n'Z34', 'G',\n'Z5A', 'X',\n'ZAF'\
, 'X',\n'ZAP', 'X',\n'ZEB', 'X',\n'ZEN', 'X',\n'ZE\
S', 'X',\n'ZID', 'X',\n'ZMR', 'X',\n'ZN3', 'X',\n'\
ZNH', 'X',\n'ZNO', 'X',\n'ZO3', 'X',\n'ZPR', 'P',\\
n'ZRA', 'A',\n'ZST', 'X',\n'ZYA', 'A',\n\n\n'ASN',\
'N');\n} \n\n\nsub file2head\n      {\n	my $file =\
 shift;\n	my $size = shift;\n	my $f= new FileHandl\
e;\n	my $line;\n	open ($f,$file);\n	read ($f,$line\
, $size);\n	close ($f);\n	return $line;\n      }\n\
sub file2tail\n      {\n	my $file = shift;\n	my $s\
ize = shift;\n	my $f= new FileHandle;\n	my $line;\\
n	\n	open ($f,$file);\n	seek ($f,$size*-1, 2);\n	r\
ead ($f,$line, $size);\n	close ($f);\n	return $lin\
e;\n      }\n\n\nsub vtmpnam\n      {\n	my $r=rand\
(100000);\n	my $f=\"file.$r.$$\";\n	while (-e $f)\\
n	  {\n	    $f=vtmpnam();\n	  }\n	push (@TMPFILE_L\
IST, $f);\n	return $f;\n      }\n\nsub myexit\n  {\
\n    my $code=@_[0];\n    if ($CLEAN_EXIT_STARTED\
==1){return;}\n    else {$CLEAN_EXIT_STARTED=1;}\n\
    ### ONLY BARE EXIT\n    exit ($code);\n  }\nsu\
b set_error_lock\n    {\n      my $name = shift;\n\
      my $pid=$$;\n\n      \n      &lock4tc ($$,\"\
LERROR\", \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\
\n\");\n      return;\n    }\nsub set_lock\n  {\n \
   my $pid=shift;\n    my $msg= shift;\n    my $p=\
getppid();\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\\
",\"$p$msg\\n\");\n  }\nsub unset_lock\n   {\n    \
 \n    my $pid=shift;\n    &lock4tc ($pid,\"LLOCK\\
",\"LRELEASE\",\"\");\n  }\nsub shift_lock\n  {\n \
   my $from=shift;\n    my $to=shift;\n    my $fro\
m_type=shift;\n    my $to_type=shift;\n    my $act\
ion=shift;\n    my $msg;\n    \n    if (!&lock4tc(\
$from, $from_type, \"LCHECK\", \"\")){return 0;}\n\
    $msg=&lock4tc ($from, $from_type, \"LREAD\", \\
"\");\n    &lock4tc ($from, $from_type,\"LRELEASE\\
", $msg);\n    &lock4tc ($to, $to_type, $action, $\
msg);\n    return;\n  }\nsub isshellpid\n  {\n    \
my $p=shift;\n    if (!lock4tc ($p, \"LLOCK\", \"L\
CHECK\")){return 0;}\n    else\n      {\n	my $c=lo\
ck4tc($p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-SHEL\
L-/){return 1;}\n      }\n    return 0;\n  }\nsub \
isrootpid\n  {\n    if(lock4tc (getppid(), \"LLOCK\
\", \"LCHECK\")){return 0;}\n    else {return 1;}\\
n  }\nsub lock4tc\n	{\n	  my ($pid,$type,$action,$\
value)=@_;\n	  my $fname;\n	  my $host=hostname;\n\
	  \n	  if ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/\
.$pid.$host.lock4tcoffee\";}\n	  elsif ( $type eq \
\"LERROR\"){ $fname=\"$LOCKDIR/.$pid.$host.error4t\
coffee\";}\n	  elsif ( $type eq \"LWARNING\"){ $fn\
ame=\"$LOCKDIR/.$pid.$host.warning4tcoffee\";}\n	 \
 \n	  if ($debug_lock)\n	    {\n	      print STDER\
R \"\\n\\t---lock4tc(tcg): $action => $fname =>$va\
lue (RD: $LOCKDIR)\\n\";\n	    }\n\n	  if    ($act\
ion eq \"LCHECK\") {return -e $fname;}\n	  elsif (\
$action eq \"LREAD\"){return file2string($fname);}\
\n	  elsif ($action eq \"LSET\") {return string2fi\
le ($value, $fname, \">>\");}\n	  elsif ($action e\
q \"LRESET\") {return string2file ($value, $fname,\
 \">\");}\n	  elsif ($action eq \"LRELEASE\") \n	 \
   {\n	      if ( $debug_lock)\n		{\n		  my $g=new\
 FileHandle;\n		  open ($g, \">>$fname\");\n		  pr\
int $g \"\\nDestroyed by $$\\n\";\n		  close ($g);\
\n		  safe_system (\"mv $fname $fname.old\");\n		}\
\n	      else\n		{\n		  unlink ($fname);\n		}\n	  \
  }\n	  return \"\";\n	}\n	\nsub file2string\n	{\n\
	  my $file=@_[0];\n	  my $f=new FileHandle;\n	  m\
y $r;\n	  open ($f, \"$file\");\n	  while (<$f>){$\
r.=$_;}\n	  close ($f);\n	  return $r;\n	}\nsub st\
ring2file \n    {\n    my ($s,$file,$mode)=@_;\n  \
  my $f=new FileHandle;\n    \n    open ($f, \"$mo\
de$file\");\n    print $f  \"$s\";\n    close ($f)\
;\n  }\n\nBEGIN\n    {\n      srand;\n    \n      \
$SIG{'SIGUP'}='signal_cleanup';\n      $SIG{'SIGIN\
T'}='signal_cleanup';\n      $SIG{'SIGQUIT'}='sign\
al_cleanup';\n      $SIG{'SIGILL'}='signal_cleanup\
';\n      $SIG{'SIGTRAP'}='signal_cleanup';\n     \
 $SIG{'SIGABRT'}='signal_cleanup';\n      $SIG{'SI\
GEMT'}='signal_cleanup';\n      $SIG{'SIGFPE'}='si\
gnal_cleanup';\n      \n      $SIG{'SIGKILL'}='sig\
nal_cleanup';\n      $SIG{'SIGPIPE'}='signal_clean\
up';\n      $SIG{'SIGSTOP'}='signal_cleanup';\n   \
   $SIG{'SIGTTIN'}='signal_cleanup';\n      $SIG{'\
SIGXFSZ'}='signal_cleanup';\n      $SIG{'SIGINFO'}\
='signal_cleanup';\n      \n      $SIG{'SIGBUS'}='\
signal_cleanup';\n      $SIG{'SIGALRM'}='signal_cl\
eanup';\n      $SIG{'SIGTSTP'}='signal_cleanup';\n\
      $SIG{'SIGTTOU'}='signal_cleanup';\n      $SI\
G{'SIGVTALRM'}='signal_cleanup';\n      $SIG{'SIGU\
SR1'}='signal_cleanup';\n\n\n      $SIG{'SIGSEGV'}\
='signal_cleanup';\n      $SIG{'SIGTERM'}='signal_\
cleanup';\n      $SIG{'SIGCONT'}='signal_cleanup';\
\n      $SIG{'SIGIO'}='signal_cleanup';\n      $SI\
G{'SIGPROF'}='signal_cleanup';\n      $SIG{'SIGUSR\
2'}='signal_cleanup';\n\n      $SIG{'SIGSYS'}='sig\
nal_cleanup';\n      $SIG{'SIGURG'}='signal_cleanu\
p';\n      $SIG{'SIGCHLD'}='signal_cleanup';\n    \
  $SIG{'SIGXCPU'}='signal_cleanup';\n      $SIG{'S\
IGWINCH'}='signal_cleanup';\n      \n      $SIG{'I\
NT'}='signal_cleanup';\n      $SIG{'TERM'}='signal\
_cleanup';\n      $SIG{'KILL'}='signal_cleanup';\n\
      $SIG{'QUIT'}='signal_cleanup';\n      \n    \
  our $debug_lock=$ENV{\"DEBUG_LOCK\"};\n      \n \
     \n      \n      \n      foreach my $a (@ARGV)\
{$CL.=\" $a\";}\n      if ( $debug_lock ){print ST\
DERR \"\\n\\n\\n********** START PG: $PROGRAM ****\
*********\\n\";}\n      if ( $debug_lock ){print S\
TDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR\
 $$ *************\\n\";}\n      if ( $debug_lock )\
{print STDERR \"\\n --- $$ -- $CL\\n\";}\n      \n\
	     \n      \n      \n    }\nsub flush_error\n  \
{\n    my $msg=shift;\n    return add_error ($EXIT\
_FAILURE,$$, $$,getppid(), $msg, $CL);\n  }\nsub a\
dd_error \n  {\n    my $code=shift;\n    my $rpid=\
shift;\n    my $pid=shift;\n    my $ppid=shift;\n \
   my $type=shift;\n    my $com=shift;\n    \n    \
$ERROR_DONE=1;\n    lock4tc ($rpid, \"LERROR\",\"L\
SET\",\"$pid -- ERROR: $type\\n\");\n    lock4tc (\
$$, \"LERROR\",\"LSET\", \"$pid -- COM: $com\\n\")\
;\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pid --\
 STACK: $ppid -> $pid\\n\");\n   \n    return $cod\
e;\n  }\nsub add_warning \n  {\n    my $rpid=shift\
;\n    my $pid =shift;\n    my $command=shift;\n  \
  my $msg=\"$$ -- WARNING: $command\\n\";\n    pri\
nt STDERR \"$msg\";\n    lock4tc ($$, \"LWARNING\"\
, \"LSET\", $msg);\n  }\n\nsub signal_cleanup\n  {\
\n    print dtderr \"\\n**** $$ (tcg) was killed\\\
n\";\n    &cleanup;\n    exit ($EXIT_FAILURE);\n  \
}\nsub clean_dir\n  {\n    my $dir=@_[0];\n    if \
( !-d $dir){return ;}\n    elsif (!($dir=~/tmp/)){\
return ;}#safety check 1\n    elsif (($dir=~/\\*/)\
){return ;}#safety check 2\n    else\n      {\n	`r\
m -rf $dir`;\n      }\n    return;\n  }\nsub clean\
up\n  {\n    #print stderr \"\\n----tc: $$ Kills $\
PIDCHILD\\n\";\n    #kill (SIGTERM,$PIDCHILD);\n  \
  my $p=getppid();\n    $CLEAN_EXIT_STARTED=1;\n  \
  \n    \n    \n    if (&lock4tc($$,\"LERROR\", \"\
LCHECK\", \"\"))\n      {\n	my $ppid=getppid();\n	\
if (!$ERROR_DONE) \n	  {\n	    &lock4tc($$,\"LERRO\
R\", \"LSET\", \"$$ -- STACK: $p -> $$\\n\");\n	  \
  &lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- COM: $\
CL\\n\");\n	  }\n      }\n    my $warning=&lock4tc\
($$, \"LWARNING\", \"LREAD\", \"\");\n    my $erro\
r=&lock4tc($$,  \"LERROR\", \"LREAD\", \"\");\n   \
 #release error and warning lock if root\n    \n  \
  if (isrootpid() && ($warning || $error) )\n     \
 {\n	\n	print STDERR \"**************** Summary **\
***********\\n$error\\n$warning\\n\";\n\n	&lock4tc\
($$,\"LERROR\",\"RELEASE\",\"\");\n	&lock4tc($$,\"\
LWARNING\",\"RELEASE\",\"\");\n      } \n    \n   \
 \n    foreach my $f (@TMPFILE_LIST)\n      {\n	if\
 (-e $f){unlink ($f);} \n      }\n    foreach my $\
d (@TMPDIR_LIST)\n      {\n	clean_dir ($d);\n     \
 }\n    #No More Lock Release\n    #&lock4tc($$,\"\
LLOCK\",\"LRELEASE\",\"\"); #release lock \n\n    \
if ( $debug_lock ){print STDERR \"\\n\\n\\n*******\
*** END PG: $PROGRAM ($$) *************\\n\";}\n  \
  if ( $debug_lock ){print STDERR \"\\n\\n\\n*****\
*****(tcg) LOCKDIR: $LOCKDIR $$ *************\\n\"\
;}\n  }\nEND \n  {\n    \n    &cleanup();\n  }\n  \
 \n\nsub safe_system \n{\n  my $com=shift;\n  my $\
ntry=shift;\n  my $ctry=shift;\n  my $pid;\n  my $\
status;\n  my $ppid=getppid();\n  if ($com eq \"\"\
){return 1;}\n  \n  \n\n  if (($pid = fork ()) < 0\
){return (-1);}\n  if ($pid == 0)\n    {\n      se\
t_lock($$, \" -SHELL- $com (tcg)\");\n      exec (\
$com);\n    }\n  else\n    {\n      lock4tc ($$, \\
"LLOCK\", \"LSET\", \"$pid\\n\");#update parent\n \
     $PIDCHILD=$pid;\n    }\n  if ($debug_lock){pr\
intf STDERR \"\\n\\t .... safe_system (fasta_seq2h\
mm)  p: $$ c: $pid COM: $com\\n\";}\n\n  waitpid (\
$pid,WTERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNI\
NG\",\"LWARNING\", \"LSET\");\n\n  if ($? == $EXIT\
_FAILURE || lock4tc($pid, \"LERROR\", \"LCHECK\", \
\"\"))\n    {\n      if ($ntry && $ctry <$ntry)\n	\
{\n	  add_warning ($$,$$,\"$com failed [retry: $ct\
ry]\");\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\\
", \"\");\n	  return safe_system ($com, $ntry, ++$\
ctry);\n	}\n      elsif ($ntry == -1)\n	{\n	  if (\
!shift_lock ($pid, $$, \"LERROR\", \"LWARNING\", \\
"LSET\"))\n	    {\n	      add_warning ($$,$$,\"$co\
m failed\");\n	    }\n	  else\n	    {\n	      lock\
4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	    \
}\n	  return $?;}\n      else\n	{\n	  if (!shift_l\
ock ($pid,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	\
    {\n	      myexit(add_error ($EXIT_FAILURE,$$,$\
pid,getppid(), \"UNSPECIFIED system\", $com));\n	 \
   }\n	}\n    }\n  return $?;\n}\n\nsub check_conf\
iguration \n    {\n      my @l=@_;\n      my $v;\n\
      foreach my $p (@l)\n	{\n	  \n	  if   ( $p eq\
 \"EMAIL\")\n	    { \n	      if ( !($EMAIL=~/@/))\\
n		{\n		add_warning($$,$$,\"Could Not Use EMAIL\")\
;\n		myexit(add_error ($EXIT_FAILURE,$$,$$,getppid\
(),\"EMAIL\",\"$CL\"));\n	      }\n	    }\n	  elsi\
f( $p eq \"INTERNET\")\n	    {\n	      if ( !&chec\
k_internet_connection())\n		{\n		  myexit(add_erro\
r ($EXIT_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$C\
L\"));\n		}\n	    }\n	  elsif( $p eq \"wget\")\n	 \
   {\n	      if (!&pg_is_installed (\"wget\") && !\
&pg_is_installed (\"curl\"))\n		{\n		  myexit(add_\
error ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INST\
ALLED:wget\",\"$CL\"));\n		}\n	    }\n	  elsif( !(\
&pg_is_installed ($p)))\n	    {\n	      myexit(add\
_error ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INS\
TALLED:$p\",\"$CL\"));\n	    }\n	}\n      return 1\
;\n    }\nsub pg_is_installed\n  {\n    my @ml=@_;\
\n    my $r, $p, $m;\n    my $supported=0;\n    \n\
    my $p=shift (@ml);\n    if ($p=~/::/)\n      {\
\n	if (safe_system (\"perl -M$p -e 1\")==$EXIT_SUC\
CESS){return 1;}\n	else {return 0;}\n      }\n    \
else\n      {\n	$r=`which $p 2>/dev/null`;\n	if ($\
r eq \"\"){return 0;}\n	else {return 1;}\n      }\\
n  }\n\n\n\nsub check_internet_connection\n  {\n  \
  my $internet;\n    my $tmp;\n    &check_configur\
ation ( \"wget\"); \n    \n    $tmp=&vtmpnam ();\n\
    \n    if     (&pg_is_installed    (\"wget\")){\
`wget www.google.com -O$tmp >/dev/null 2>/dev/null\
`;}\n    elsif  (&pg_is_installed    (\"curl\")){`\
curl www.google.com -o$tmp >/dev/null 2>/dev/null`\
;}\n    \n    if ( !-e $tmp || -s $tmp < 10){$inte\
rnet=0;}\n    else {$internet=1;}\n    if (-e $tmp\
){unlink $tmp;}\n\n    return $internet;\n  }\nsub\
 check_pg_is_installed\n  {\n    my @ml=@_;\n    m\
y $r=&pg_is_installed (@ml);\n    if (!$r && $p=~/\
::/)\n      {\n	print STDERR \"\\nYou Must Install\
 the perl package $p on your system.\\nRUN:\\n\\ts\
udo perl -MCPAN -e 'install $pg'\\n\";\n      }\n \
   elsif (!$r)\n      {\n	myexit(flush_error(\"\\n\
Program $p Supported but Not Installed on your sys\
tem\"));\n      }\n    else\n      {\n	return 1;\n\
      }\n  }\n\n\nsub remote_is_pdb_name_deprecate\
d\n{\n    my $in=@_[0];\n    my ($ref_file, $pdb);\
\n    my ($value,$value1,$value2);\n    my $max=2;\
\n    \n    \n    \n    $ref_file=\"$cache/pdb_ent\
ry_type.txt\";\n    \n    if ( $in=~/[^\\w\\d\\:\\\
_]/){return 0;}\n    elsif ($no_remote_pdb_dir==1)\
\n      {\n	my $pdbdir=$ENV{'PDB_DIR'};\n	\n	my $r\
1=\"$pdbdir/derived_data/pdb_entry_type.txt\";\n	m\
y $r2=$ref_file;\n	if    (-e $r1){$ref_file=$r1;}\\
n	elsif (-e $r2){$ref_file=$r2;}\n	else\n	  {\n	  \
  my $p=substr ($in,0, 4);\n	    add_warning ($$, \
$$, \"Cannot find pdb_entry_type.txt;  $p is assum\
ed to be valid; add ftp://ftp.wwpdb.org/pub/pdb/de\
rived_data/pdb_entry_type.txt in $cache to check n\
ame status\");\n	  }\n      }\n    elsif ( !-e $re\
f_file || (-M $ref_file)>$max || -z $ref_file)\n  \
    {\n	&url2file(\"ftp://ftp.wwpdb.org/pub/pdb/de\
rived_data/pdb_entry_type.txt\", $ref_file);\n    \
  }\n    $pdb=substr ($in,0, 4);\n    chomp(($valu\
e1=`grep -c $pdb $ref_file`));\n    $pdb=lc($pdb);\
\n    chomp(($value2=`grep -c $pdb $ref_file`));\n\
    $value=($value1 || $value2)?1:0;\n    $value=(\
$value>0)?1:0;\n    \n    return $value;\n  }\n","\
use Cwd;\nuse Env;\nuse File::Path;\nuse FileHandl\
e;\nuse strict;\n\n\nour (%MODE, %PG, %ENV_SET, %S\
UPPORTED_OS);\n\nour $VERSION=\"0\"; #_#UPDATE_VER\
SION\n\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\
\nour $INTERNET=0;\n\nour $CP=\"cp \"; #was causin\
g a crash on MacOSX\nour $SILENT=\">/dev/null 2>/d\
ev/null\";\nour $WEB_BASE=\"http://www.tcoffee.org\
\";\nour $TCLINKDB_ADDRESS=\"$WEB_BASE/Resources/t\
clinkdb.txt\";\nour $OS=get_os();\nour $ROOT=&get_\
root();\nour $CD=cwd();\nour $CDIR=$CD;\nour $HOME\
=$ENV{'HOME'};\n\nour $OSNAME=$ENV{'OSNAME'};\nour\
 $OSARCH=$ENV{'OSARCH'};\nour $REPO_ROOT=\"\";\n\n\
our $TCDIR;\nour $TCCACHE;\nour $TCTMP;\nour $TCM;\
\nour $TCMETHODS;\nour $TCPLUGINS;\nour $PLUGINS_D\
IR=\"\";\nour $INSTALL_DIR=\"\";\nour $email;\nour\
 $recompile;\n\nour $CXX=\"g++\";\nour $CXXFLAGS=\\
"\";\n\nour $CPP=\"g++\";\nour $CPPFLAGS=\"\";\n\n\
our $CC=\"gcc\";\nour $CFLAGS=$ENV{'CFLAGS'};\n\no\
ur $FC=\"f77\";\nour $FFLAGS=\"\";\n\nmy $install=\
\"all\";\nmy $default_update_action=\"no_update\";\
\nmy @required_applications=(\"wget_OR_curl\");\nm\
y @smode=(\"all\", \"clean\", \"install\");\n\n&in\
itialize_PG();\nmy $cl=join( \" \", @ARGV);\nif ($\
#ARGV==-1 || ($cl=~/-h/) ||($cl=~/-H/) )\n  {\n   \
  print \"\\n!!!!!!! ./install  t_coffee          \
   --> installs t_coffee only\";\n     print \"\\n\
!!!!!!! ./install  all                  --> instal\
ls all the modes [mcoffee, expresso, psicoffee,rco\
ffee..]\";\n     print \"\\n!!!!!!! ./install  [mc\
offee|rcoffee|..] --> installs the specified mode\\
";\n     print \"\\n!!!!!!! ./install  -h         \
          --> print usage\\n\\n\";\n     if ( $#AR\
GV==-1){exit ($EXIT_FAILURE);}\n   }\n     \nif ((\
$cl=~/-h/) ||($cl=~/-H/) )\n  {\n    my $m;\n    p\
rint \"\\n\\n!!!!!!! advanced mode\\n\";\n    fore\
ach $m ((keys (%MODE)),@smode)\n      {\n	print \"\
!!!!!!!       ./install $m\\n\";\n      }\n    \n \
   print \"!!!!!!! ./install [target:package|mode|\
] [-update|-force|-exec=dir|-dis=dir|-root|-tclink\
db=file|-] [CC=|FCC=|CXX=|CFLAGS=|CXXFLAGS=]\\n\";\
\n    print \"!!!!!!! ./install clean    [removes \
all executables]\\n\";\n    print \"!!!!!!! ./inst\
all [optional:target] -update               [updat\
es package already installed]\\n\";\n    print \"!\
!!!!!! ./install [optional:target] -recompile     \
       [forces the recompilation of T-Coffee]\\n\"\
;\n\n    print \"!!!!!!! ./install [optional:targe\
t] -force                [Forces recompilation ove\
r everything]\\n\";\n    \n    print \"!!!!!!! ./i\
nstall [optional:target] -root                 [Yo\
u are running as root]\\n\";\n    print \"!!!!!!! \
./install [optional:target] -exec=/foo/bar/       \
[address for the T-Coffee executable]\\n\";\n    p\
rint \"!!!!!!! ./install [optional:target] -dis=/f\
oo/bar/        [Address where distributions should\
 be stored]\\n\";\n    print \"!!!!!!! ./install [\
optional:target] -tclinkdb=foo|update  [file conta\
ining all the packages to be installed]\\n\";\n   \
 print \"!!!!!!! ./install [optional:target] -clea\
n                [clean everything]\\n\";\n    pri\
nt \"!!!!!!! ./install [optional:target] -plugins \
             [plugins directory]\\n\";\n    print \
\"!!!!!!! ./install [optional:target] -tcdir=/foor\
/bar      [base path where T-Coffee will be instal\
led - default ~/.t_coffee]\\n\";\n    print \"!!!!\
!!! ./install [optional:target] -repo=/path/to/rep\
o   [binaries repository root directory]\\n\";\n  \
  print \"!!!!!!! ./install [optional:target] -ema\
il=<your email>   [needed for remote BLAST]\\n\";\\
n    print \"!!!!!!! ./install [optional:target] -\
proxy=<proxy>   [may be needed to access remote se\
rvices]\\n\";\n    \n    print \"!!!!!!! mode:\";\\
n    foreach $m (keys(%MODE)){print \"$m \";}\n   \
 print \"\\n\";\n    print \"!!!!!!! Packages:\";\\
n    foreach $m (keys (%PG)){print \"$m \";}\n    \
print \"\\n\";\n    \n    print \"\\n\\n\";\n    e\
xit ($EXIT_FAILURE);\n  }\n\n\n\nmy (@argl)=($cl=~\
/(\\S+=[^=]+)\\s\\w+=/g);\npush (@argl, ($cl=~/(\\\
S+=[^=]+\\S)\\s*$/g));\n\nforeach $a (@argl)\n  {\\
n    if ( ($cl=~/CXX=(.*)/)){$CXX=$1;}\n    if ( (\
$cl=~/-CC=(.*)/    )){$CC=$1;}\n    if ( ($cl=~/-F\
C=(.*)/    )){$FC=$1;}\n    if ( ($cl=~/-CFLAGS=(.\
*)/)){$CFLAGS=$1;}\n    if ( ($cl=~/-CXXFLAGS=(.*)\
/)){$CXXFLAGS=$1;}\n  }\nour ($ROOT_INSTALL, $NO_Q\
UESTION, $default_update_action,$BINARIES_ONLY,$fo\
rce, $default_update_action, $INSTALL_DIR, $PLUGIN\
S_DIR, $DISTRIBUTIONS,$tclinkdb, $proxy, $clean);\\
nif ( ($cl=~/-root/)){$ROOT_INSTALL=1;}\nif ( ($cl\
=~/-no_question/)){$NO_QUESTION=1;}\nif ( ($cl=~/-\
update/)){$default_update_action=\"update\";}\nif \
( ($cl=~/-recompile/)){$recompile=1;}\n\n\n$BINARI\
ES_ONLY=1;\n\nif ( ($cl=~/-nobinaries/)){$BINARIES\
_ONLY=0;}\nif ( ($cl=~/-force/)){$force=1;$default\
_update_action=\"update\"}\nif ( ($cl=~/-exec=\\s*\
(\\S+)/)){$INSTALL_DIR=$1;}\nif ( ($cl=~/-plugins=\
\\s*(\\S+)/)){$PLUGINS_DIR=$1;}\nif ( ($cl=~/-dis=\
\\s*(\\S+)/)){$DISTRIBUTIONS=$1;}\n\nif ( ($cl=~/-\
tclinkdb=\\s*(\\S+)/)){$tclinkdb=$1;}\nif ( ($cl=~\
/-proxy=\\s*(\\S+)/)){$proxy=$1;}\nif ( ($cl=~/-cl\
ean/)){$clean=1;}\nif ( ($cl=~/-repo=\\s*(\\S+)/))\
{ $REPO_ROOT=$1; }\nif ( ($cl=~/-tcdir=\\s*(\\S+)/\
)){ $TCDIR=$1; }\n\nif ( ($cl=~/-email=\\s*(\\S+)/\
)){$email=$1;}\n\n\nif ($tclinkdb){&update_tclinkd\
b ($tclinkdb);}\n\n\nif( $REPO_ROOT ne \"\" ) {\n	\
if( $OSNAME eq \"\" ) { print \"You have specified\
 the repository folder but the required \\\"OSNAME\
\\\" enviroment variable is missing. \\n\"; exit 1\
; } \n	if( $OSARCH eq \"\" ) { print \"You have sp\
ecified the repository folder but the required \\\\
"OSARCH\\\" enviroment variable is missing. \\n\";\
 exit 1; } \n}\n\n\nif(!$TCDIR) { $TCDIR=\"$HOME/.\
t_coffee\"; }\n&add_dir ($TCDIR);\n&add_dir ($TCCA\
CHE=\"$TCDIR/cache\");\n&add_dir ($TCTMP=\"$CDIR/t\
mp\");\n&add_dir ($TCM=\"$TCDIR/mcoffee\");\n&add_\
dir ($TCMETHODS=\"$TCDIR/methods\");\n&add_dir ($T\
CPLUGINS=\"$TCDIR/plugins/$OS\");\n\n\nour $BASE=\\
"$CD/bin\";\nour $BIN=\"$BASE/cache/binaries/$OS\"\
;\nour $DOWNLOAD_DIR=\"$BASE/cache/download\";\nou\
r $DOWNLOAD_FILE=\"$DOWNLOAD_DIR/files\";\nour $TM\
P=\"$BASE/cache/tmp\";\n\n&add_dir($BASE);\n&add_d\
ir($BIN);\n&add_dir($DOWNLOAD_DIR);\n&add_dir($DOW\
NLOAD_FILE);\nif (!$DISTRIBUTIONS){$DISTRIBUTIONS=\
\"$DOWNLOAD_DIR/distributions\";}\n&add_dir ($DIST\
RIBUTIONS);\n&add_dir ($TMP);\n\n\nif    (!$PLUGIN\
S_DIR && !$ROOT_INSTALL){$PLUGINS_DIR=$TCPLUGINS;}\
\nelsif (!$PLUGINS_DIR &&  $ROOT_INSTALL){$PLUGINS\
_DIR=\"/usr/local/bin/\";}\n\nif    (!$INSTALL_DIR\
 && !$ROOT_INSTALL){$INSTALL_DIR=\"$TCDIR/bin/$OS\\
";mkpath ($INSTALL_DIR);}\nelsif (!$INSTALL_DIR &&\
  $ROOT_INSTALL){$INSTALL_DIR=\"/usr/local/bin/\";\
}\n\nif (-d \"mcoffee\"){`cp mcoffee/* $TCM`;}\n\n\
\nour $ENV_FILE=\"$TCDIR/.t_coffee_env\";\nunlink \
($ENV_FILE);\n&add2env_file ($ENV_FILE,\"EMAIL_4_T\
COFFEE\", $email);\n&add2env_file ($ENV_FILE,\"htt\
p_proxy_4_TCOFFEE\", $proxy);\n&env_file2putenv ($\
ENV_FILE);\n&set_proxy($proxy);\n\n\n\nmy ($target\
, $p, $r);\n$target=$p;\n\nforeach $p (  ((keys (%\
PG)),(keys(%MODE)),(@smode)) )\n  {\n    if ($ARGV\
[0] eq $p && $target eq \"\"){$target=$p;}\n  }\ni\
f ($target eq \"\"){exit ($EXIT_FAILURE);}\n\n\nfo\
reach $r (@required_applications)\n  {\n    my @ap\
p_list;\n    my $i;\n    $i=0;\n    \n    @app_lis\
t=split (/_OR_/, $r);\n    foreach my $pg (@app_li\
st)\n      {\n	$i+=&pg_is_installed ($pg);\n      \
}\n    if ($i==0)\n      {\n      print \"One of t\
he following packages must be installed to proceed\
: \";\n      foreach my $pg (@app_list)\n	{\n	  pr\
int (\"$pg \");\n	}\n      die;\n    }\n  }\n\n\n\\
n\n\n\n&sign_license_ni();\n\n\n$PG{C}{compiler}=g\
et_C_compiler($CC);\n$PG{Fortran}{compiler}=get_F_\
compiler($FC);\n$PG{CXX}{compiler}=$PG{CPP}{compil\
er}=$PG{GPP}{compiler}=get_CXX_compiler($CXX);\nif\
 ($CXXFLAGS){$PG{CPP}{options}=$PG{GPP}{options}=$\
PG{CXX}{options}=$CXXFLAGS;}\nif ($CFLAGS ne \"\" \
){$PG{C}{options}=$CFLAGS;}\nforeach my $c (keys(%\
PG))\n  {\n    my $arguments;\n    if ($PG{$c}{com\
piler})\n      {\n	$arguments=\"$PG{$c}{compiler_f\
lag}=$PG{$c}{compiler} \";\n	if ($PG{$c}{options})\
\n	  {\n	    $arguments.=\"$PG{$c}{options_flag}='\
\" . $PG{$c}{options} . \"' \";\n	  }\n	$PG{$c}{ar\
guments}=$arguments;\n      }\n  }\n\nif ($PG{$tar\
get}){$PG{$target}{install}=1;}\nelse\n  {\n    fo\
reach my $pg (keys(%PG))\n      {\n	if ( $target e\
q \"all\" || ($PG{$pg}{mode}=~/$target/))\n	  {\n	\
    $PG{$pg} {install}=1;\n	  }\n      }\n  }\n\nf\
oreach my $pg (keys(%PG))\n  {\n    if (!$PG{$pg}{\
update_action}){$PG{$pg}{update_action}=$default_u\
pdate_action;}\n    elsif ($PG{$pg}{update_action}\
 eq \"never\"){$PG{$pg}{install}=0;}\n    if ( $fo\
rce && $PG{$pg}{install})\n      {\n	`rm $BIN/$pg \
$BIN/$pg.exe $SILENT`;\n      }\n    if ($PG{$pg}{\
update_action} eq \"update\" && $PG{$pg}{install})\
{$PG{$pg}{update}=1;}\n  }\n\nif (($target=~/clean\
/))\n  {\n    print \"------- cleaning executables\
 -----\\n\";\n    `rm bin/* $SILENT`;\n    exit ($\
EXIT_SUCCESS);\n  }\n\nif ( !$PG{$target}){print \\
"------- Installing T-Coffee Modes\\n\";}\n\nforea\
ch my $m (keys(%MODE))\n  {\n    if ( $target eq \\
"all\" || $target eq $m)\n      {\n	print \"\\n---\
---- The installer will now install the $m compone\
nts $MODE{$m}{description}\\n\";\n	foreach my $pg \
(keys(%PG))\n	  {\n	    if ( $PG{$pg}{mode} =~/$m/\
 && $PG{$pg}{install})\n	      {\n		if ($PG{$pg}{t\
ouched}){print \"------- $PG{$pg}{dname}: already \
processed\\n\";}\n		else {$PG{$pg}{success}=&insta\
ll_pg($pg);$PG{$pg}{touched}=1;}\n	      }\n	  }\n\
      }\n  }\n\nif ( $PG{$target}){print \"-------\
 Installing Individual Package\\n\";}\nforeach my \
$pg (keys (%PG))\n  {\n    \n    if ( $PG{$pg}{ins\
tall} && !$PG{$pg}{touched})\n      {\n	print \"\\\
n------- Install $pg\\n\";\n	$PG{$pg}{success}=&in\
stall_pg($pg);$PG{$pg}{touched}=1;\n      }\n  }\n\
print \"------- Finishing The installation\\n\";\n\
my $final_report=&install ($INSTALL_DIR);\n\nprint\
 \"\\n\";\nprint \"*******************************\
**************************************\\n\";\nprin\
t \"********              INSTALLATION SUMMARY    \
      *****************\\n\";\nprint \"***********\
**************************************************\
********\\n\";\nprint \"------- SUMMARY package In\
stallation:\\n\";\nprint \"-------   Executable In\
stalled in: $PLUGINS_DIR\\n\";\n\nforeach my $pg (\
keys(%PG))\n  {\n    if ( $PG{$pg}{install})\n    \
  {\n	my $bin_status=($PG{$pg}{from_binary} && $PG\
{$pg}{success})?\"[from binary]\":\"\";\n	if     (\
 $PG{$pg}{new} && !$PG{$pg}{old})                 \
    {print \"*------        $PG{$pg}{dname}: insta\
lled $bin_status\\n\"; $PG{$pg}{status}=1;}\n	elsi\
f  ( $PG{$pg}{new} &&  $PG{$pg}{old})             \
        {print \"*------        $PG{$pg}{dname}: u\
pdated $bin_status\\n\"  ; $PG{$pg}{status}=1;} \n\
	elsif  (!$PG{$pg}{new} &&  $PG{$pg}{old} && !$PG{\
$pg}{update}){print \"*------        $PG{$pg}{dnam\
e}: previous\\n\" ; $PG{$pg}{status}=1;}\n	elsif  \
(!$PG{$pg}{new} &&  $PG{$pg}{old} &&  $PG{$pg}{upd\
ate}){print \"*------        $PG{$pg}{dname}: fail\
ed update (previous installation available)\\n\";$\
PG{$pg}{status}=0;}\n	else                        \
                                  {print \"*------\
        $PG{$pg}{dname}: failed installation\\n\";\
$PG{$pg}{status}=0;}\n      }\n  }\nmy $failure;\n\
\nif ( !$PG{$target}){print \"*------ SUMMARY mode\
 Installation:\\n\";}\nforeach my $m (keys(%MODE))\
\n  {\n  \n    if ( $target eq \"all\" || $target \
eq $m)\n      {\n	my $succesful=1;\n	foreach my $p\
g (keys(%PG))\n	  {\n	    if (($PG{$pg}{mode}=~/$m\
/) && $PG{$pg}{install} && $PG{$pg}{status}==0)\n	\
      {\n		$succesful=0;\n		print \"*!!!!!!       \
$PG{$pg}{dname}: Missing\\n\";\n	      }\n	  }\n	i\
f ( $succesful)\n	  {\n	    $MODE{$m}{status}=1;\n\
	    print \"*------       MODE $MODE{$m}{dname} S\
UCCESSFULLY installed\\n\";\n	  }\n	else\n	  {\n	 \
   $failure++;\n	    $MODE{$m}{status}=0;\n	    pr\
int \"*!!!!!!       MODE $MODE{$m}{dname} UNSUCCES\
SFULLY installed\\n\";\n	  }\n      }\n  }\n\n    \
\n      \nif ($clean==1 && ($BASE=~/install4tcoffe\
e/) ){print \"*------ Clean Installation Directory\
: $BASE\\n\";`rm -rf $BASE`;}\nforeach my $pg (key\
s(%PG)){if ($PG{$pg}{install} && $PG{$pg}{status}=\
=0){exit ($EXIT_FAILURE);}}\n\nif ($failure)\n  {\\
n    print \"*************************************\
********************************\\n\";\n    print \
\"********     SOME PACKAGES FAILED TO INSTALL    \
    *****************\\n\";\n    print \"*********\
**************************************************\
**********\\n\";\n    print \"\\nSome of the repor\
ted failures may be due to connectivity problems\"\
;\n    print \"\\nRerun the installation and the i\
nstaller will specifically try to install the miss\
ing packages\";\n    print \"\\nIf this Fails, go \
to the original website and install the package ma\
nually\";\n  }\n\nprint \"************************\
*********************************************\\n\"\
;\nprint \"********              FINALIZE YOUR INS\
TALLATION    *****************\\n\";\nprint \"****\
**************************************************\
***************\\n\";\nprint \"------- Your third \
party executables are in:\\n\"; \nprint \"------- \
      $PLUGINS_DIR:\\n\";\nprint \"------- Your t_\
coffee exccutable is in\\n\";\nprint \"-------    \
   $INSTALL_DIR:\\n\";\nprint \"------- In order t\
o make your installation permanent add these two l\
ines\\n\";\nprint \"export PATH=$INSTALL_DIR:\\$PA\
TH\\n\";\nprint \"export PLUGINS_4_TCOFFEE=$PLUGIN\
S_DIR:\\n\";\nif ($OS eq \"linux\")\n  {\n    prin\
t \"-------       to the file: $HOME/.bashrc\\n\";\
\n  }\nelse \n  {\n    print \"-------       to th\
e file: $HOME/.profile\\ OR $HOME/.basrc\";\n  }\n\
exit ($EXIT_SUCCESS);  \n  \nsub get_CXX_compiler\\
n  {\n    my $c=@_[0];\n    my (@clist)=(\"g++\");\
\n    \n    return get_compil ($c, @clist);\n }\ns\
ub get_C_compiler\n  {\n    my $c=@_[0];\n    my (\
@clist)=(\"gcc\", \"cc\", \"icc\");\n    \n    ret\
urn get_compil ($c, @clist);\n }\n\nsub get_F_comp\
iler\n  {\n    my ($c)=@_[0];\n    my @clist=(\"f7\
7\", \"g77\",\"g95\", \"gfortran\", \"ifort\");\n \
   return get_compil ($c, @clist);\n  } \n       \\
nsub get_compil\n  {\n    my ($fav,@clist)=(@_);\n\
    \n    #return the first compiler found install\
ed in the system. Check first the favorite\n    fo\
reach my $c ($fav,@clist)\n      {\n	if  (&pg_is_i\
nstalled ($c)){return $c;}\n      }\n    return \"\
\";\n  }\nsub exit_if_pg_not_installed\n  {\n    m\
y (@arg)=(@_);\n    \n    foreach my $p (@arg)\n  \
    {\n	if ( !&pg_is_installed ($p))\n	  {\n	    p\
rint \"!!!!!!!! The $p utility must be installed f\
or this installation to proceed [FATAL]\\n\";\n	  \
  die;\n	  }\n      }\n    return 1;\n  }\nsub set\
_proxy\n  {\n    my ($proxy)=(@_);\n    my (@list,\
$p);\n    \n    @list= (\"HTTP_proxy\", \"http_pro\
xy\", \"HTTP_PROXY\", \"ALL_proxy\", \"all_proxy\"\
,\"HTTP_proxy_4_TCOFFEE\",\"http_proxy_4_TCOFFEE\"\
);\n    \n    if (!$proxy)\n      {\n	foreach my $\
p (@list)\n	  {\n	    if ( ($ENV_SET{$p}) || $ENV{\
$p}){$proxy=$ENV{$p};}\n	  }\n      }\n    foreach\
 my $p(@list){$ENV{$p}=$proxy;}\n  }\n	\nsub check\
_internet_connection\n  {\n    my $internet;\n    \
\n    if ( -e \"x\"){unlink (\"x\");}\n    if     \
(&pg_is_installed    (\"wget\")){`wget www.google.\
com -Ox >/dev/null 2>/dev/null`;}\n    elsif  (&pg\
_is_installed    (\"curl\")){`curl www.google.com \
-ox >/dev/null 2>/dev/null`;}\n    else\n      {\n\
	printf stderr \"\\nERROR: No pg for remote file f\
etching [wget or curl][FATAL]\\n\";\n	exit ($EXIT_\
FAILURE);\n      }\n    \n    if ( !-e \"x\" || -s\
 \"x\" < 10){$internet=0;}\n    else {$internet=1;\
}\n    if (-e \"x\"){unlink \"x\";}\n    return $i\
nternet;\n  }\nsub url2file\n  {\n    my ($cmd, $f\
ile,$wget_arg, $curl_arg)=(@_);\n    my ($exit,$fl\
ag, $pg, $arg);\n    \n    if ($INTERNET || check_\
internet_connection ()){$INTERNET=1;}\n    else\n \
     {\n	print STDERR \"ERROR: No Internet Connect\
ion [FATAL:install.pl]\\n\";\n	exit ($EXIT_FAILURE\
);\n      }\n    \n    if     (&pg_is_installed   \
 (\"wget\")){$pg=\"wget\"; $flag=\"-O\";$arg=\"--t\
ries=2 --connect-timeout=10 --no-check-certificate\
 $wget_arg\";}\n    elsif  (&pg_is_installed    (\\
"curl\")){$pg=\"curl\"; $flag=\"-o\";$arg=$curl_ar\
g;}\n    else\n      {\n	printf stderr \"\\nERROR:\
 No pg for remote file fetching [wget or curl][FAT\
AL]\\n\";\n	exit ($EXIT_FAILURE);\n      }\n    \n\
    \n    if (-e $file){unlink($file);}\n    $exit\
=system \"$pg $cmd $flag$file $arg\";\n    return \
$exit;\n  }\n\nsub pg_is_installed\n  {\n    my ($\
p, $dir)=(@_);\n    my ($r,$m, $ret);\n    my ($su\
pported, $language, $compil);\n    \n  \n    if ( \
$PG{$p})\n      {\n	$language=$PG{$p}{language2};\\
n	$compil=$PG{$language}{compiler};\n      }\n    \
\n    if ( $compil eq \"CPAN\")\n      {\n	if ( sy\
stem (\"perl -M$p -e 1\")==$EXIT_SUCCESS){$ret=1;}\
\n	else {$ret=0;}\n      }\n    elsif ($dir)\n    \
  {\n	if (-e \"$dir/$p\" || -e \"$dir/$p\\.exe\"){\
$ret=1;}\n	else {$ret=0;}\n      }\n    elsif (-e \
\"$PLUGINS_DIR/$p\" || -e \"$PLUGINS_DIR/$p.exe\")\
{$ret=1;}\n    else\n      {\n	$r=`which $p 2>/dev\
/null`;\n	if ($r eq \"\"){$ret=0;}\n	else {$ret=1;\
}\n      }\n   \n    return $ret;\n  }\nsub instal\
l\n  {\n    my ($new_bin)=(@_);\n    my ($copied, \
$report);\n\n    \n    if (!$ROOT_INSTALL)\n      \
{\n	`$CP $BIN/* $PLUGINS_DIR`;\n	if (-e \"$BIN/t_c\
offee\")\n	  {\n	    `$CP $BIN/t_coffee $INSTALL_D\
IR`;\n	      unlink(\"$PLUGINS_DIR/t_coffee\");\n	\
  }\n	$copied=1;\n      }\n    else\n      {\n	$co\
pied=&root_run (\"You must be root to finalize the\
 installation\", \"$CP $BIN/* $PLUGINS_DIR $SILENT\
\");\n	if (-e \"$BIN/t_coffee\")\n	  {\n	    &root\
_run (\"You must be root to finalize the installat\
ion\", \"$CP $BIN/t_coffee $INSTALL_DIR\");\n	    \
&root_run (\"You must be root to finalize the inst\
allation\", \"rm  $PLUGINS_DIR/t_coffee\");\n	  }\\
n      }\n    \n     \n  if ( !$copied)\n    {\n  \
    $report=\"*!!!!!! Installation unsuccesful. Th\
e executables have been left in $BASE/bin\\n\";\n \
   }\n  elsif ( $copied && $ROOT)\n    {\n      $r\
eport=\"*------ Installation succesful. Your execu\
tables have been copied in $new_bin and are on you\
r PATH\\n\";\n    }\n  elsif ( $copied && !$ROOT)\\
n    {\n      $report= \"*!!!!!! T-Coffee has been\
 installed in $INSTALL_DIR\\n\";\n      $report= \\
"*!!!!!! T-Coffee and associated packages have bee\
n copied in: $PLUGINS_DIR\\n\";\n      $report.=\"\
*!!!!!! This T-Coffee location is NOT on your PATH\
 sytem variable\\n\";\n      if ( $OS eq \"linux\"\
)\n	{\n	  $report.=\"*!!!!!! You can do so by addi\
ng the following line in your ~/.bashrc file:\\n\"\
;\n	}\n      else\n	{\n	  $report.=\"*!!!!!! You c\
an do so by adding the following line in your ~/.p\
rofile file:\\n\";\n	}\n      $report.=\"*!!!!!! e\
xport PATH=$INSTALL_DIR:\\$PATH\\n\";\n    }\n  re\
turn $report;\n}\n\nsub sign_license_ni\n  {\n    \
my $F=new FileHandle;\n    open ($F, \"license.txt\
\");\n    while (<$F>)\n      {\n	print \"$_\";\n \
     }\n    close ($F);\n    \n    return;\n  }\n\\
nsub install_pg\n  {\n    my ($pg)=(@_);\n    my (\
$report, $previous, $language, $compiler, $return)\
;\n    \n    if (!$PG{$pg}{install}){return 1;}\n \
   \n    $previous=&pg_is_installed ($pg);\n    \n\
    if ($PG{$pg}{update_action} eq \"no_update\" &\
& $previous)\n      {\n	$PG{$pg}{old}=1;\n	$PG{$pg\
}{new}=0;\n	$return=1;\n      }\n    else\n      {\
\n	$PG{$pg}{old}=$previous;\n	\n	if ($PG{$pg} {lan\
guage2} eq \"Perl\"){&install_perl_package ($pg);}\
\n	elsif ($BINARIES_ONLY && &install_binary_packag\
e ($pg)){$PG{$pg}{from_binary}=1;}\n	elsif (&insta\
ll_source_package ($pg)){;}\n	else \n	  {\n	    \n\
	    if (!&supported_os($OS))\n	      {\n		print \\
"!!!!!!!! $pg compilation failed, binary unsupport\
ed for $OS\\n\"; \n	      }\n	    elsif (!($PG{$pg\
}{from_binary}=&install_binary_package ($pg)))\n	 \
     {\n		print \"!!!!!!!! $pg compilation and  bi\
nary installation failed\\n\";\n	      }\n	  }\n	$\
PG{$pg}{new}=$return=&pg_is_installed ($pg,$BIN);\\
n      }\n\n    \n    return $return;\n  }\nsub in\
stall_perl_package\n  {\n    my ($pg)=(@_);\n    m\
y ($report, $language, $compiler);\n    \n    $lan\
guage=$PG{$pg} {language2};\n    $compiler=$PG{$la\
nguage}{compiler};\n    \n    if (!&pg_is_installe\
d ($pg))\n      {\n	if ( $OS eq \"windows\"){`perl\
 -M$compiler -e 'install $pg'`;}\n	elsif ( $ROOT e\
q \"sudo\"){system (\"sudo perl -M$compiler -e 'in\
stall $pg'\");}\n	else {system (\"su root -c perl \
-M$compiler -e 'install $pg'\");}\n      }\n    re\
turn &pg_is_installed ($pg);\n  }\n\n\n\nsub insta\
ll_source_package\n  {\n    my ($pg)=(@_);\n    my\
 ($report, $download, $arguments, $language, $addr\
ess, $name, $ext, $main_dir, $distrib);\n    my $w\
get_tmp=\"$TMP/wget.tmp\";\n    my (@fl);\n    if \
( $default_update_action ne \"update\" && (-e \"$B\
IN/$pg\" || -e \"$BIN/$pg.exe\" )  ){return 1;}\n \
   \n    #\n    # check if the module exists in th\
e repository cache \n    #\n	if( repo_load($pg) ) \
{\n		return 1;\n	}\n    \n    if ($pg eq \"t_coffe\
e\")  {return   &install_t_coffee_source ($pg);}\n\
    elsif ($pg eq \"TMalign\"){return   &install_T\
Malign ($pg);}\n    \n    chdir $DISTRIBUTIONS;\n \
   \n    $download=$PG{$pg}{source};\n    \n    if\
 (($download =~/tgz/))\n      {\n	($address,$name,\
$ext)=($download=~/(.+\\/)([^\\/]+)(\\.tgz).*/);\n\
      }\n    elsif (($download=~/tar\\.gz/))\n    \
  {\n	($address,$name,$ext)=($download=~/(.+\\/)([\
^\\/]+)(\\.tar\\.gz).*/);\n      }\n    elsif (($d\
ownload=~/tar/))\n      {\n	($address,$name,$ext)=\
($download=~/(.+\\/)([^\\/]+)(\\.tar).*/);\n      \
}\n    else\n      {\n	($address,$name)=($download\
=~/(.+\\/)([^\\/]+)/);\n	$ext=\"\";\n      }\n    \
$distrib=\"$name$ext\";\n    \n    if ( !-d $pg){m\
kdir $pg;}\n    chdir $pg;\n   \n    #get the dist\
ribution if available\n    if ( -e \"$DOWNLOAD_DIR\
/$distrib\")\n      {\n	`$CP $DOWNLOAD_DIR/$distri\
b .`;\n      }\n    #UNTAR and Prepare everything\\
n    if (!-e \"$name.tar\" && !-e \"$name\")\n    \
  {\n	&check_rm ($wget_tmp);\n	print \"\\n------- \
Downloading/Installing $pg\\n\";\n	\n	if (!-e $dis\
trib && &url2file (\"$download\", \"$wget_tmp\")==\
$EXIT_SUCCESS)\n	  {\n	    \n	    `mv $wget_tmp $d\
istrib`;\n	    `$CP $distrib $DOWNLOAD_DIR/`;\n	  \
}\n\n	if (!-e $distrib)\n	  {\n	    print \"!!!!!!\
! Download of $pg distribution failed\\n\";\n	    \
print \"!!!!!!! Check Address: $PG{$pg}{source}\\n\
\";\n	    return 0;\n	  }\n	print \"\\n------- unz\
ipping/untaring $name\\n\";\n	if (($ext =~/z/))\n	\
  { \n	    &flush_command (\"gunzip -f $name$ext\"\
);\n	    \n	  }\n	if (($ext =~/tar/) || ($ext =~/t\
gz/))\n	  {\n	    &flush_command(\"tar -xvf $name.\
tar\");\n	  }\n      }\n    #Guess and enter the d\
istribution directory\n    @fl=ls($p);\n    foreac\
h my $f (@fl)\n      {\n	if (-d $f)\n	  {\n	    $m\
ain_dir=$f;\n	  }\n      }\n    if (-d $main_dir)\\
n	  \n      {\n	chdir $main_dir;}\n    else\n     \
 {\n	print \"Error: $main_dir does not exist\";\n \
     }\n    print \"\\n------- Compiling/Installin\
g $pg\\n\";\n    `make clean $SILENT`;\n    \n    \
\n    #\n    # SAP module\n    #\n    if ($pg eq \\
"sap\")\n      {\n	if (-e \"./configure\")\n	  {\n\
	    #new sap distribution\n	    \n	    &flush_com\
mand (\"./configure\");\n	    &flush_command (\"ma\
ke clean\");\n	    &flush_command (\"make\");\n	  \
  &check_cp (\"./src/$pg\", \"$BIN\");\n	    repo_\
store(\"./src/$pg\");\n	  }\n	else\n	  {\n	    #ol\
d style distribution\n	    `rm *.o sap  sap.exe ./\
util/aa/*.o  ./util/wt/.o $SILENT`;\n	    &flush_c\
ommand (\"make $arguments sap\");\n	    &check_cp \
($pg, \"$BIN\");\n	    repo_store($pg);\n	  }\n   \
   }\n    \n    #\n    # CLUSTALW2 module\n    #\n\
    elsif ($pg eq \"clustalw2\")\n      {\n	&flush\
_command(\"./configure\");\n	&flush_command(\"make\
 $arguments\");\n	&check_cp (\"./src/$pg\", \"$BIN\
\");\n	repo_store(\"./src/$pg\");\n      }\n\n    \
#\n    # CLUSTAL-OMEGA module\n    #\n    elsif ($\
pg eq \"clustalo\")\n      {\n	&flush_command(\"./\
configure\");\n	&flush_command(\"make $arguments\"\
);\n	&check_cp (\"./src/$pg\", \"$BIN\");\n	repo_s\
tore(\"./src/$pg\");\n      }\n\n    #\n    # STRI\
KE module\n    #\n    elsif ($pg eq \"strike\")\n \
     {\n	&flush_command(\"make $arguments\");\n	&c\
heck_cp (\"./bin/$pg\", \"$BIN\");\n	repo_store(\"\
./bin/$pg\");\n      }\n    \n    #\n    # FSA mod\
ule\n    # \n    elsif ($pg eq \"fsa\")\n      {\n\
	&flush_command(\"./configure --prefix=$BIN\");\n	\
&flush_command(\"make $arguments\");\n	&flush_comm\
and (\"make install\");\n\n	repo_store(\"fsa\", \"\
$BIN/bin\");\n	`mv $BIN/bin/* $BIN`;\n	`rmdir $BIN\
/bin`;\n      }\n    \n    #\n    # CLUSTALW modul\
e\n    #\n    elsif ($pg eq \"clustalw\")\n      {\
\n	&flush_command(\"make $arguments clustalw\");\n\
	`$CP $pg $BIN $SILENT`;\n	repo_store($pg);\n     \
 }\n    \n    #\n    # MAFFT module\n    #\n    el\
sif ($pg eq \"mafft\")\n      {\n	my $base=cwd();\\
n	my $c;\n	\n	#compile core\n	mkpath (\"./mafft/bi\
n\");\n	mkpath (\"./mafft/lib\");\n	chdir \"$base/\
core\";\n	`make clean $SILENT`;\n	&flush_command (\
\"make $arguments\");\n	&flush_command (\"make ins\
tall LIBDIR=../mafft/lib BINDIR=../mafft/bin\");\n\
	\n	#compile extension\n	chdir \"$base/extensions\\
";\n	`make clean $SILENT`;\n	&flush_command (\"mak\
e $arguments\");\n	&flush_command (\"make install \
LIBDIR=../mafft/lib BINDIR=../mafft/bin\");\n	\n	#\
put everything in mafft and copy the compiled stuf\
f in bin\n	chdir \"$base\";\n	if ($ROOT_INSTALL)\n\
	  {\n	    &root_run (\"You Must be Root to Instal\
l MAFFT\\n\", \"mkdir /usr/local/mafft/;$CP mafft/\
lib/* /usr/local/mafft;$CP mafft/lib/mafft* /usr/l\
ocal/bin ;$CP mafft/bin/mafft /usr/local/bin/; \")\
;\n	  }\n	else\n	  {\n	    `$CP mafft/lib/*  $BIN`\
;\n	    `$CP mafft/bin/mafft  $BIN`;\n	  }\n	`tar \
-cvf mafft.tar mafft`;\n	`gzip mafft.tar`;\n	`mv m\
afft.tar.gz $BIN`;\n	\n	repo_store(\"mafft/bin/maf\
ft\", \"mafft/lib/\", \"$BIN/mafft.tar.gz\");\n   \
   }\n      \n    #\n    # DIALIGN-TX module\n    \
#\n    elsif ( $pg eq \"dialign-tx\" )\n      {\n	\
my $f;\n	my $base=cwd();\n\n	chdir \"./source\";\n\
	if ($OS eq \"macosx\"){&flush_command (\"cp makef\
ile.MAC_OS makefile\");}\n\n	&flush_command (\" ma\
ke CPPFLAGS='-O3 -funroll-loops' all\");\n	\n	chdi\
r \"..\";\n	&check_cp (\"./source/$pg\", \"$BIN\")\
;\n	repo_store(\"./source/$pg\");\n      }\n      \
\n    #\n    # DIALIGN-T module \n    # (is the sa\
me as dialign-tx, but it is mantained for backward\
 name compatibility with tcoffee)\n    #\n    elsi\
f ( $pg eq \"dialign-t\" )\n      {\n	my $f;\n	my \
$base=cwd();\n\n	chdir \"./source\";\n	if ($OS eq \
\"macosx\"){&flush_command (\"cp makefile.MAC_OS m\
akefile\");}\n\n	&flush_command (\" make CPPFLAGS=\
'-O3 -funroll-loops' all\");\n	\n	chdir \"..\";\n	\
&check_cp (\"./source/dialign-tx\", \"$BIN/dialign\
-t\");\n	repo_store(\"$BIN/dialign-t\");	\n      }\
      \n      \n    #\n    # POA module\n    #\n  \
  elsif ($pg eq \"poa\")\n      {\n	&flush_command\
 (\"make $arguments poa\");\n	&check_cp (\"$pg\", \
\"$BIN\");\n	repo_store(\"$pg\");\n      }\n     \\
n     \n    #\n    # PROBCONS module\n    #\n    e\
lsif ( $pg eq \"probcons\")\n      {\n	&add_C_libr\
aries(\"./ProbabilisticModel.h\", \"list\", \"cstr\
ing\");\n	\n	`rm *.exe $SILENT`;\n	&flush_command \
(\"make $arguments probcons\");\n	&check_cp(\"$pg\\
", \"$BIN/$pg\");\n	repo_store(\"$pg\");\n      }\\
n      \n    #\n    # PROBCONS RNA module\n    #\n\
    elsif ( $pg eq \"probconsRNA\")\n      {\n	&ad\
d_C_libraries(\"./ProbabilisticModel.h\", \"list\"\
, \"cstring\");\n	&add_C_libraries(\"./Main.cc\", \
\"iomanip\", \"cstring\",\"climits\");\n	`rm *.exe\
 $SILENT`;\n	&flush_command (\"make $arguments pro\
bcons\");\n	&check_cp(\"probcons\", \"$BIN/$pg\");\
\n	repo_store(\"$BIN/$pg\");\n      }\n\n	#\n	# MU\
SCLE module\n	#\n    elsif (  $pg eq \"muscle\")\n\
      {	\n	`rm *.o muscle muscle.exe $SILENT`;\n	i\
f ($OS eq \"macosx\" || $OS eq \"linux\")\n	  {\n	\
    &replace_line_in_file (\"./Makefile\", \"LDLIB\
S = -lm -static\",  \"LDLIBS = -lm\");\n	  }\n	els\
if ($OS eq \"windows\")\n	  {\n	    &replace_line_\
in_file (\"./intmath.cpp\",  \"double log2e\",    \
  \"double cedric_log\");\n	    &replace_line_in_f\
ile (\"./intmath.cpp\",  \"double log2\",       \"\
double log_notuse\");\n	    &replace_line_in_file \
(\"./intmath.cpp\",  \"double cedric_log\", \"doub\
le log2e\");\n	  }\n	&flush_command (\"make $argum\
ents all\");\n	&check_cp(\"$pg\", \"$BIN\");\n	rep\
o_store(\"$pg\");	\n      }\n      \n     #\n     \
# MUS4 module\n     #\n     elsif (  $pg eq \"mus4\
\")\n      {\n	`rm *.o muscle muscle.exe $SILENT`;\
\n	&flush_command (\"./mk\");\n	&check_cp(\"$pg\",\
 \"$BIN\");\n	repo_store(\"$pg\");	\n      }\n    \
  \n    #\n    # PCMA module\n    #\n    elsif ( $\
pg eq \"pcma\")\n      {\n	if ($OS eq \"macosx\")\\
n	  {\n	    &replace_line_in_file (\"./alcomp2.c\"\
, \"malloc.h\",  \"\");\n	  }\n	&flush_command (\"\
make $arguments pcma\");\n	&check_cp(\"$pg\", \"$B\
IN\");\n	repo_store(\"$pg\");	\n      }\n      \n \
   #\n    # KALIGN module\n    #\n    elsif ($pg e\
q \"kalign\")\n      {\n	&flush_command (\"./confi\
gure\");\n	&flush_command(\"make $arguments\");\n	\
&check_cp (\"$pg\",$BIN);\n	repo_store(\"$pg\");	\\
n      }\n      \n    #\n    # AMAP module\n    #\\
n    elsif ( $pg eq \"amap\")\n      {\n	&add_C_li\
braries(\"./Amap.cc\", \"iomanip\", \"cstring\",\"\
climits\");	\n	`make clean $SILENT`;\n	&flush_comm\
and (\"make $arguments all\");\n	&check_cp (\"$pg\\
", $BIN);\n	repo_store(\"$pg\");	\n      }\n      \
\n    #\n    # PRODA module\n    #\n    elsif ( $p\
g eq \"proda\")\n      {\n	`sed -i '' 's/int errno\
 = 0;/int errno; errno = 0;/' Main.cc`;\n	&add_C_l\
ibraries(\"AlignedFragment.h\", \"vector\", \"iost\
ream\", \"cstring\",\"cstdlib\");\n	&add_C_librari\
es(\"Main.cc\", \"vector\", \"climits\");	\n	&add_\
C_libraries(\"Sequence.cc\", \"stdlib.h\", \"cstdi\
o\");	\n	&flush_command (\"make $arguments all\");\
\n	&check_cp (\"$pg\", $BIN);\n	repo_store(\"$pg\"\
);	\n      }\n      \n    #\n    # PRANK module\n \
   #\n    elsif ( $pg eq \"prank\")\n      {\n	&fl\
ush_command (\"make $arguments all\");\n	&check_cp\
 (\"$pg\", $BIN);\n	repo_store(\"$pg\");	\n      }\
\n      \n    #\n    # !!!! MUSTANG module\n    #\\
n     elsif ( $pg eq \"mustang\")\n      {\n	&flus\
h_command (\"rm ./bin/*\");\n	&flush_command (\"ma\
ke $arguments all\");\n\n	if ( $OS=~/windows/){&fl\
ush_command(\"cp ./bin/* $BIN/mustang.exe\");}\n	e\
lse {&flush_command(\"cp ./bin/* $BIN/mustang\");}\
\n	\n	repo_store(\"$BIN/mustang\");\n      }\n\n	#\
\n	# RNAplfold module\n	#\n    elsif ( $pg eq \"RN\
Aplfold\")\n      {\n	&flush_command(\"./configure\
\");\n	&flush_command (\"make $arguments all\");\n\
	&check_cp(\"./Progs/RNAplfold\", \"$BIN\");\n	&ch\
eck_cp(\"./Progs/RNAalifold\", \"$BIN\");\n	&check\
_cp(\"./Progs/RNAfold\", \"$BIN\");\n	\n	repo_stor\
e(\"./Progs/RNAplfold\", \"./Progs/RNAalifold\", \\
"./Progs/RNAfold\");\n      }\n      \n    #\n    \
# !!! RETREE module\n    #\n    elsif ( $pg eq \"r\
etree\")\n      {\n	chdir \"src\";\n	&flush_comman\
d (\"cp Makefile.unx Makefile\");\n	&flush_command\
 (\"make $arguments all\");\n	&flush_command (\"ma\
ke put\");\n	system \"cp ../exe/* $BIN\";\n	\n	rep\
o_store(\"retree\", \"../exe\");\n      }\n	\n    \
chdir $CDIR;\n    return &pg_is_installed ($pg, $B\
IN);\n  }\n\nsub install_t_coffee_source\n  {\n   \
 my ($pg)=(@_);\n    my ($report,$cflags, $argumen\
ts, $language, $compiler) ;\n\n    #1-Install T-Co\
ffee\n    chdir \"t_coffee_source\";\n    &flush_c\
ommand (\"make clean\");\n    print \"\\n------- C\
ompiling T-Coffee\\n\";\n    $language=$PG{$pg} {l\
anguage2};\n    $arguments=$PG{$language}{argument\
s};\n    \n    if ( $CC ne \"\")\n      {\n	print \
\"make -i $arguments t_coffee \\n\";\n	&flush_comm\
and (\"make -i $arguments t_coffee\");\n      }\n \
   &check_cp ($pg, $BIN);\n    \n    chdir $CDIR;\\
n    return &pg_is_installed ($pg, $BIN);\n  }\nsu\
b install_TMalign\n  {\n    my ($pg)=(@_);\n    my\
 $report;\n    chdir \"t_coffee_source\";\n    pri\
nt \"\\n------- Compiling TMalign\\n\";\n    `rm T\
Malign TMalign.exe $SILENT`;\n    if ( $FC ne \"\"\
){&flush_command (\"make -i $PG{Fortran}{arguments\
} TMalign\");}\n    &check_cp ($pg, $BIN);\n    re\
po_store($pg);\n\n    if ( !-e \"$BIN/$pg\" && pg_\
has_binary_distrib ($pg))\n      {\n	print \"!!!!!\
!! Compilation of $pg impossible. Will try to inst\
all from binary\\n\";\n	return &install_binary_pac\
kage ($pg);\n      }\n    chdir $CDIR;\n    return\
 &pg_is_installed ($pg, $BIN);\n  }\n\nsub pg_has_\
binary_distrib\n  {\n    my ($pg)=(@_);\n    if ($\
PG{$pg}{windows}){return 1;}\n    elsif ($PG{$pg}{\
osx}){return 1;}\n    elsif ($PG{$pg}{macosx}){ret\
urn 1;}\n\n    elsif ($PG{$pg}{linux}){return 1;}\\
n    return 0;\n  }\nsub install_binary_package\n \
 {\n    my ($pg)=(@_);\n    my ($base,$report,$nam\
e, $download, $arguments, $language, $dir);\n    m\
y $isdir;\n    &input_os();\n    \n    #\n    # - \
paolodt - Check if the module exists in the reposi\
tory cache \n    #\n	if( repo_load($pg) ) {\n	    \
$PG{$pg}{from_binary}=1;\n		return 1;\n	}\n    # -\
 paolodt - end \n    \n    if (!&supported_os($OS)\
){return 0;}\n    if ( $PG{$pg}{binary}){$name=$PG\
{$pg}{binary};}\n    else {$name=$pg;}\n    if ($n\
ame eq \"t_coffee\")\n      {\n	#check if local bi\
n is there\n	if (-e \"./bin/$OS/t_coffee\")\n	  {\\
n	    print \"\\n------- Installing  T-Coffee from\
 Pre-Compiled/Pre-Downloaded $OS binary\\n\";\n	  \
  print \"\\n------- If you want to trigger a fres\
h compilation use -recompile\\n\";\n	    &check_cp\
 (\"./bin/$OS/t_coffee\", $BIN);\n	    return &pg_\
is_installed ($pg, $BIN);\n	  }\n	#try to get prec\
ompiled binary -- available from MAC is distributi\
on from MAC\n	else\n	  {\n	    $download=\"$WEB_BA\
SE/Packages/Binaries/tcoffee/$OS/$name.$VERSION\";\
\n	  }\n      }\n    else\n      {\n	$download=\"$\
WEB_BASE/Packages/Binaries/plugins/$OS/$name\";\n \
     }\n    \n    $base=cwd();\n    chdir $TMP;\n \
   \n    if (!-e $name)\n      {\n	`rm x $SILENT`;\
\n	if ( url2file(\"$download\",\"x\")==$EXIT_SUCCE\
SS)\n	  {\n	    `mv x $name`;\n	  }\n      }\n    \
\n    if (!-e $name)\n      {\n	print \"!!!!!!! $P\
G{$pg}{dname}: Download of $pg binary failed\\n\";\
\n	print \"!!!!!!! $PG{$pg}{dname}: Check Address:\
 $download\\n\";\n	chdir $base;\n	return 0;\n     \
 }\n    print \"\\n------- Installing $pg\\n\";\n \
   \n    if ($name =~/tar\\.gz/)\n      {\n	`gunzi\
p  -f $name`;\n	`tar -xvf $pg.tar`;\n	chdir $pg;\n\
	`chmod u+x *`;\n 	`mv * $BIN`;\n	#if (!($pg=~/\\*\
/)){`rm -rf $pg`;}\n      }\n    else\n      {\n	&\
check_cp (\"$pg\", \"$BIN\");\n	`chmod u+x $BIN/$p\
g`; \n	unlink ($pg);\n      }\n    chdir $base;\n \
   $PG{$pg}{from_binary}=1;\n\n    return &pg_is_i\
nstalled ($pg, $BIN);\n  }\n\n	\nsub add_dir\n  {\\
n    my $dir=@_[0];\n    \n    if (!-e $dir && !-d\
 $dir)\n      {\n	my @l;\n	umask (0000);\n	@l=mkpa\
th ($dir,{mode => 0777});\n	\n      }\n    else\n \
     {\n	return 0;\n      }\n  }\nsub check_rm \n \
 {\n    my ($file)=(@_);\n    \n    if ( -e $file)\
\n      {\n	return unlink($file);\n      }\n    re\
turn 0;\n  }\nsub check_cp\n  {\n    my ($from, $t\
o)=(@_);\n    if ( !-e $from && -e \"$from\\.exe\"\
){$from=\"$from\\.exe\";}\n    if ( !-e $from){ret\
urn 0;}\n        \n    `$CP $from $to`;\n    retur\
n 1;\n  }\n\nsub repo_store \n{\n   # check that a\
ll required data are available\n   if( $REPO_ROOT \
eq \"\" ) { return; }\n\n\n    # extract the packa\
ge name from the specified path\n    my $pg =`base\
name $_[0]`;\n    chomp($pg);\n	\n    my $VER = $P\
G{$pg}{version};\n    my $CACHE = \"$REPO_ROOT/$pg\
/$VER/$OSNAME-$OSARCH\"; \n    \n    print \"-----\
--- Storing package: \\\"$pg\\\" to path: $CACHE\\\
n\";\n    \n    # clean the cache path if exists a\
nd create it again\n    `rm -rf $CACHE`;\n    `mkd\
ir -p $CACHE`;\n    \n 	for my $path (@_) {\n\n	  \
  # check if it is a single file \n	 	if( -f $path\
 ) {\n	    	`cp $path $CACHE`;\n		}\n		# .. or a d\
irectory, in this case copy all the content \n		el\
sif( -d $path ) {\n			opendir(IMD, $path);\n			my \
@thefiles= readdir(IMD);\n			closedir(IMD);\n			\n\
			for my $_file (@thefiles) {\n				if( $_file ne \
\".\" && $_file ne \"..\") {\n	    			`cp $path/$_\
file $CACHE`;\n				}\n			}\n		} \n	}	   \n    \n	\\
n}   \n\nsub repo_load \n{\n    my ($pg)=(@_);\n\n\
    #Bypass the Repository Cache\n    return 0;\n \
   # check that all required data are available\n \
   if( $REPO_ROOT eq \"\" ) { return 0; }\n\n    m\
y $VER = $PG{$pg}{version};\n    my $CACHE = \"$RE\
PO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n    if( !-e \
\"$CACHE/$pg\" ) {\n   	 	print \"-------- Module \
\\\"$pg\\\" NOT found on repository cache.\\n\";\n\
    	return 0;\n    }\n    \n    print \"-------- \
Module \\\"$pg\\\" found on repository cache. Usin\
g copy on path: $CACHE\\n\";\n    `cp $CACHE/* $BI\
N`;\n    return 1;\n}\n\nsub check_file_list_exist\
s \n  {\n    my ($base, @flist)=(@_);\n    my $f;\\
n\n    foreach $f (@flist)\n      {\n	if ( !-e \"$\
base/$f\"){return 0;}\n      }\n    return 1;\n  }\
\nsub ls\n  {\n    my $f=@_[0];\n    my @fl;\n    \
chomp(@fl=`ls -1 $f`);\n    return @fl;\n  }\nsub \
flush_command\n  {\n    my $command=@_[0];\n    my\
 $F=new FileHandle;\n    open ($F, \"$command|\");\
\n    while (<$F>){print \"    --- $_\";}\n    clo\
se ($F);\n  }    \n\nsub input_installation_direct\
ory\n  {\n    my $dir=@_[0];\n    my $new;\n    \n\
    print \"------- The current installation direc\
tory is: [$dir]\\n\";\n    print \"??????? Return \
to keep the default or new value:\";\n   \n    if \
($NO_QUESTION==0)\n      {\n	chomp ($new=<stdin>);\
\n	while ( $new ne \"\" && !input_yes (\"You have \
entered $new. Is this correct? ([y]/n):\"))\n	  {\\
n	    print \"???????New installation directory:\"\
;\n	    chomp ($new=<stdin>);\n	  }\n	$dir=($new e\
q \"\")?$dir:$new;\n	$dir=~s/\\/$//;\n      }\n   \
 \n    if ( -d $dir){return $dir;}\n    elsif (&ro\
ot_run (\"You must be root to create $dir\",\"mkdi\
r $dir\")==$EXIT_SUCCESS){return $dir;}\n    else\\
n      {\n	print \"!!!!!!! $dir could not be creat\
ed\\n\";\n	if ( $NO_QUESTION)\n	  {\n	    return \\
"\";\n	  }\n	elsif ( &input_yes (\"??????? Do you \
want to provide a new directory([y]/n)?:\"))\n	  {\
\n	    return input_installation_directory ($dir);\
\n	  }\n	else\n	  {\n	    return \"\";\n	  }\n    \
  }\n    \n  }\nsub input_yes\n  {\n    my $questi\
on =@_[0];\n    my $answer;\n\n    if ($NO_QUESTIO\
N==1){return 1;}\n    \n    if ($question eq \"\")\
{$question=\"??????? Do you wish to proceed ([y]/n\
)?:\";}\n    print $question;\n    chomp($answer=l\
c(<STDIN>));\n    if (($answer=~/^y/) || $answer e\
q \"\"){return 1;}\n    elsif ( ($answer=~/^n/)){r\
eturn 0;}\n    else\n      {\n	return input_yes($q\
uestion);\n      }\n  }\nsub root_run\n  {\n    my\
 ($txt, $cmd)=(@_);\n    \n    if ( system ($cmd)=\
=$EXIT_SUCCESS){return $EXIT_SUCCESS;}\n    else \\
n      {\n	print \"------- $txt\\n\";\n	if ( $ROOT\
 eq \"sudo\"){return system (\"sudo $cmd\");}\n	el\
se {return system (\"su root -c \\\"$cmd\\\"\");}\\
n      }\n  }\nsub get_root\n  {\n    if (&pg_is_i\
nstalled (\"sudo\")){return \"sudo\";}\n    else {\
return \"su\";}\n  }\n\nsub get_os\n  {\n    my $r\
aw_os=`uname`;\n    my $os;\n\n    $raw_os=lc ($ra\
w_os);\n    \n    if ($raw_os =~/cygwin/){$os=\"wi\
ndows\";}\n    elsif ($raw_os =~/linux/){$os=\"lin\
ux\";}\n    elsif ($raw_os =~/osx/){$os=\"macosx\"\
;}\n    elsif ($raw_os =~/darwin/){$os=\"macosx\";\
}\n    else\n      {\n	$os=$raw_os;\n      }\n    \
return $os;\n  }\nsub input_os\n  {\n    my $answe\
r;\n    if ($OS) {return $OS;}\n    \n    print \"\
??????? which os do you use: [w]indows, [l]inux, [\
m]acosx:?\";\n    $answer=lc(<STDIN>);\n\n    if (\
($answer=~/^m/)){$OS=\"macosx\";}\n    elsif ( ($a\
nswer=~/^w/)){$OS=\"windows\";}\n    elsif ( ($ans\
wer=~/^linux/)){$OS=\"linux\";}\n    \n    else\n \
     {\n	return &input_os();\n      }\n    return \
$OS;\n  }\n\nsub supported_os\n  {\n    my ($os)=(\
@_[0]);\n    return $SUPPORTED_OS{$os};\n  }\n\nsu\
b add2env_file\n  {\n    my ($env, $var, $value)=(\
@_);\n    my $F = new FileHandle;\n    my $t;\n   \
 if (!$value){return;}\n    #make sure new variabl\
es do not get duplicated\n    if ( -e $env)\n     \
 {\n	open ($F, \"$env\");\n	while (<$F>)\n	  {\n	 \
   my $line=$_;\n	    if (!($line=~/$var/)){$t.=$l\
ine;}\n	  }\n	close ($F);\n      }\n    $t.=\"$var\
=$value\\n\";\n    open ($F, \">$env\");\n    prin\
t $F \"$t\";\n    $ENV{$var}=$value;\n    close ($\
F);\n  }    \n    \n\n\nsub update_tclinkdb \n  {\\
n    my $file =@_[0];\n    my $name;\n    my $F=ne\
w FileHandle;\n    my ($download, $address, $name,\
 $l, $db);\n    \n    if ( $file eq \"update\"){$f\
ile=$TCLINKDB_ADDRESS;}\n    \n    if ( $file =~/h\
ttp:\\/\\// || $file =~/ftp:\\/\\//)\n      {\n	($\
address, $name)=($download=~/(.*)\\/([^\\/]+)$/);\\
n	`rm x $SILENT`;\n	if (&url2file ($file,\"x\")==$\
EXIT_SUCCESS)\n	  {\n	    print \"------- Susscess\
ful upload of $name\";\n	    `mv x $name`;\n	    $\
file=$name;\n	  }\n      }\n    open ($F, \"$file\\
");\n    while (<$F>)\n      {\n	my $l=$_;\n	if ((\
$l =~/^\\/\\//) || ($db=~/^#/)){;}\n	elsif ( !($l \
=~/\\w/)){;}\n	else\n	  {\n	    my @v=split (/\\s+\
/, $l);\n	    if ( $l=~/^MODE/)\n	      {\n		$MODE\
{$v[1]}{$v[2]}=$v[3];\n	      }\n	    elsif ($l=~/\
^PG/)\n	      {\n		$PG{$v[1]}{$v[2]}=$v[3];\n	    \
  }\n	  }\n      }\n    close ($F);\n    &post_pro\
cess_PG();\n    return;\n  }\n\n\n\nsub initialize\
_PG\n  {\n\n$PG{\"t_coffee\"}{\"4_TCOFFEE\"}=\"TCO\
FFEE\";\n$PG{\"t_coffee\"}{\"type\"}=\"sequence_mu\
ltiple_aligner\";\n$PG{\"t_coffee\"}{\"ADDRESS\"}=\
\"http://www.tcoffee.org\";\n$PG{\"t_coffee\"}{\"l\
anguage\"}=\"C++\";\n$PG{\"t_coffee\"}{\"language2\
\"}=\"CXX\";\n$PG{\"t_coffee\"}{\"source\"}=\"http\
://www.tcoffee.org/Packages/sources/tcoffee/stable\
/T-COFFEE_distribution.tar.gz\";\n$PG{\"t_coffee\"\
}{\"update_action\"}=\"always\";\n$PG{\"t_coffee\"\
}{\"binary\"}=\"t_coffee\";\n$PG{\"t_coffee\"}{\"m\
ode\"}=\"tcoffee,mcoffee,rcoffee,expresso,3dcoffee\
\";\n$PG{\"clustalo\"}{\"4_TCOFFEE\"}=\"CLUSTALO\"\
;\n$PG{\"clustalo\"}{\"type\"}=\"sequence_multiple\
_aligner\";\n$PG{\"clustalo\"}{\"ADDRESS\"}=\"http\
://www.clustal.org/omega/\";\n$PG{\"clustalo\"}{\"\
language\"}=\"C++\";\n$PG{\"clustalo\"}{\"language\
2\"}=\"C++\";\n$PG{\"clustalo\"}{\"source\"}=\"htt\
p://www.clustal.org/omega/clustal-omega-1.2.4.tar.\
gz\";\n$PG{\"clustalo\"}{\"mode\"}=\"mcoffee\";\n$\
PG{\"clustalo\"}{\"binary\"}=\"clustalo\";\n$PG{\"\
clustalo\"}{\"version\"}=\"1.2.4\";\n$PG{\"strike\\
"}{\"4_TCOFFEE\"}=\"STRIKE\";\n$PG{\"strike\"}{\"t\
ype\"}=\"sequence_alignment_scoring\";\n$PG{\"stri\
ke\"}{\"ADDRESS\"}=\"http://www.tcoffee.org/Projec\
ts/strike/index.html\";\n$PG{\"strike\"}{\"languag\
e\"}=\"C++\";\n$PG{\"strike\"}{\"language2\"}=\"CX\
X\";\n$PG{\"strike\"}{\"source\"}=\"http://www.tco\
ffee.org/Projects/strike/strike_v1.2.tar.bz2\";\n$\
PG{\"strike\"}{\"mode\"}=\"tcoffee,expresso\";\n$P\
G{\"strike\"}{\"version\"}=\"1.2\";\n$PG{\"strike\\
"}{\"binary\"}=\"strike\";\n$PG{\"clustalw2\"}{\"4\
_TCOFFEE\"}=\"CLUSTALW2\";\n$PG{\"clustalw2\"}{\"t\
ype\"}=\"sequence_multiple_aligner\";\n$PG{\"clust\
alw2\"}{\"ADDRESS\"}=\"http://www.clustal.org\";\n\
$PG{\"clustalw2\"}{\"language\"}=\"C++\";\n$PG{\"c\
lustalw2\"}{\"language2\"}=\"CXX\";\n$PG{\"clustal\
w2\"}{\"source\"}=\"http://www.clustal.org/downloa\
d/2.0.10/clustalw-2.0.10-src.tar.gz\";\n$PG{\"clus\
talw2\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"clu\
stalw2\"}{\"binary\"}=\"clustalw2\";\n$PG{\"clusta\
lw2\"}{\"version\"}=\"2.0.10\";\n$PG{\"clustalw\"}\
{\"4_TCOFFEE\"}=\"CLUSTALW\";\n$PG{\"clustalw\"}{\\
"type\"}=\"sequence_multiple_aligner\";\n$PG{\"clu\
stalw\"}{\"ADDRESS\"}=\"http://www.clustal.org\";\\
n$PG{\"clustalw\"}{\"language\"}=\"C\";\n$PG{\"clu\
stalw\"}{\"language2\"}=\"C\";\n$PG{\"clustalw\"}{\
\"source\"}=\"http://www.clustal.org/download/1.X/\
ftp-igbmc.u-strasbg.fr/pub/ClustalW/clustalw1.82.U\
NIX.tar.gz\";\n$PG{\"clustalw\"}{\"mode\"}=\"mcoff\
ee,rcoffee\";\n$PG{\"clustalw\"}{\"version\"}=\"1.\
82\";\n$PG{\"clustalw\"}{\"binary\"}=\"clustalw\";\
\n$PG{\"dialign-t\"}{\"4_TCOFFEE\"}=\"DIALIGNT\";\\
n$PG{\"dialign-t\"}{\"type\"}=\"sequence_multiple_\
aligner\";\n$PG{\"dialign-t\"}{\"ADDRESS\"}=\"http\
://dialign-tx.gobics.de/\";\n$PG{\"dialign-t\"}{\"\
DIR\"}=\"/usr/share/dialign-tx/\";\n$PG{\"dialign-\
t\"}{\"language\"}=\"C\";\n$PG{\"dialign-t\"}{\"la\
nguage2\"}=\"C\";\n$PG{\"dialign-t\"}{\"source\"}=\
\"http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar\
.gz\";\n$PG{\"dialign-t\"}{\"mode\"}=\"mcoffee\";\\
n$PG{\"dialign-t\"}{\"binary\"}=\"dialign-t\";\n$P\
G{\"dialign-t\"}{\"version\"}=\"1.0.2\";\n$PG{\"di\
align-tx\"}{\"4_TCOFFEE\"}=\"DIALIGNTX\";\n$PG{\"d\
ialign-tx\"}{\"type\"}=\"sequence_multiple_aligner\
\";\n$PG{\"dialign-tx\"}{\"ADDRESS\"}=\"http://dia\
lign-tx.gobics.de/\";\n$PG{\"dialign-tx\"}{\"DIR\"\
}=\"/usr/share/dialign-tx/\";\n$PG{\"dialign-tx\"}\
{\"language\"}=\"C\";\n$PG{\"dialign-tx\"}{\"langu\
age2\"}=\"C\";\n$PG{\"dialign-tx\"}{\"source\"}=\"\
http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.g\
z\";\n$PG{\"dialign-tx\"}{\"mode\"}=\"mcoffee\";\n\
$PG{\"dialign-tx\"}{\"binary\"}=\"dialign-tx\";\n$\
PG{\"dialign-tx\"}{\"version\"}=\"1.0.2\";\n$PG{\"\
poa\"}{\"4_TCOFFEE\"}=\"POA\";\n$PG{\"poa\"}{\"typ\
e\"}=\"sequence_multiple_aligner\";\n$PG{\"poa\"}{\
\"ADDRESS\"}=\"http://www.bioinformatics.ucla.edu/\
poa/\";\n$PG{\"poa\"}{\"language\"}=\"C\";\n$PG{\"\
poa\"}{\"language2\"}=\"C\";\n$PG{\"poa\"}{\"sourc\
e\"}=\"http://downloads.sourceforge.net/poamsa/poa\
V2.tar.gz\";\n$PG{\"poa\"}{\"DIR\"}=\"/usr/share/\\
";\n$PG{\"poa\"}{\"FILE1\"}=\"blosum80.mat\";\n$PG\
{\"poa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"poa\"}{\"b\
inary\"}=\"poa\";\n$PG{\"poa\"}{\"version\"}=\"2.0\
\";\n$PG{\"probcons\"}{\"4_TCOFFEE\"}=\"PROBCONS\"\
;\n$PG{\"probcons\"}{\"type\"}=\"sequence_multiple\
_aligner\";\n$PG{\"probcons\"}{\"ADDRESS\"}=\"http\
://probcons.stanford.edu/\";\n$PG{\"probcons\"}{\"\
language2\"}=\"CXX\";\n$PG{\"probcons\"}{\"languag\
e\"}=\"C++\";\n$PG{\"probcons\"}{\"source\"}=\"htt\
p://probcons.stanford.edu/probcons_v1_12.tar.gz\";\
\n$PG{\"probcons\"}{\"mode\"}=\"mcoffee\";\n$PG{\"\
probcons\"}{\"binary\"}=\"probcons\";\n$PG{\"probc\
ons\"}{\"version\"}=\"1.12\";\n$PG{\"msaprobs\"}{\\
"4_TCOFFEE\"}=\"MSAPROBS\";\n$PG{\"msaprobs\"}{\"t\
ype\"}=\"sequence_multiple_aligner\";\n$PG{\"msapr\
obs\"}{\"ADDRESS\"}=\"http://msaprobs.sourceforge.\
net/homepage.htm#latest\";\n$PG{\"msaprobs\"}{\"la\
nguage2\"}=\"CXX\";\n$PG{\"msaprobs\"}{\"language\\
"}=\"C++\";\n$PG{\"msaprobs\"}{\"source\"}=\"https\
://sourceforge.net/projects/msaprobs/files/MSAProb\
s-MPI/MSAProbs-MPI_rel1.0.5.tar.gz\";\n$PG{\"msapr\
obs\"}{\"mode\"}=\"mcoffee\";\n$PG{\"msaprobs\"}{\\
"binary\"}=\"msaprobs\";\n$PG{\"msaprobs\"}{\"vers\
ion\"}=\"1.05\";\n$PG{\"msaprobs\"}{\"update_actio\
n\"}=\"never\";\n$PG{\"upp\"}{\"4_TCOFFEE\"}=\"UPP\
\";\n$PG{\"upp\"}{\"type\"}=\"sequence_multiple_al\
igner\";\n$PG{\"upp\"}{\"ADDRESS\"}=\"http://www.c\
s.utexas.edu/users/phylo/software/upp/\";\n$PG{\"u\
pp\"}{\"language2\"}=\"CXX\";\n$PG{\"upp\"}{\"lang\
uage\"}=\"C++\";\n$PG{\"upp\"}{\"source\"}=\"https\
://github.com/smirarab/pasta/archive/upp.zip\";\n$\
PG{\"upp\"}{\"mode\"}=\"mcoffee\";\n$PG{\"upp\"}{\\
"binary\"}=\"upp\";\n$PG{\"upp\"}{\"version\"}=\"1\
\";\n$PG{\"upp\"}{\"update_action\"}=\"never\";\n$\
PG{\"famsa\"}{\"4_TCOFFEE\"}=\"FAMSA\";\n$PG{\"fam\
sa\"}{\"type\"}=\"sequence_multiple_aligner\";\n$P\
G{\"famsa\"}{\"ADDRESS\"}=\"https://github.com/ref\
resh-bio/FAMSA\";\n$PG{\"famsa\"}{\"language\"}=\"\
C++\";\n$PG{\"famsa\"}{\"language\"}=\"C++\";\n$PG\
{\"famsa\"}{\"source\"}=\"https://github.com/refre\
sh-bio/FAMSA.git\";\n$PG{\"famsa\"}{\"mode\"}=\"mc\
offee,rcoffee\";\n$PG{\"famsa\"}{\"binary\"}=\"fam\
sa\";\n$PG{\"famsa\"}{\"version\"}=\"1.1\";\n$PG{\\
"mafft\"}{\"4_TCOFFEE\"}=\"MAFFT\";\n$PG{\"mafft\"\
}{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"\
mafft\"}{\"ADDRESS\"}=\"http://align.bmr.kyushu-u.\
ac.jp/mafft/online/server/\";\n$PG{\"mafft\"}{\"la\
nguage\"}=\"C\";\n$PG{\"mafft\"}{\"language\"}=\"C\
\";\n$PG{\"mafft\"}{\"source\"}=\"http://mafft.cbr\
c.jp/alignment/software/mafft-7.310-with-extension\
s-src.tgz\";\n$PG{\"mafft\"}{\"mode\"}=\"mcoffee,r\
coffee\";\n$PG{\"mafft\"}{\"binary\"}=\"mafft.tar.\
gz\";\n$PG{\"mafft\"}{\"version\"}=\"7.310\";\n$PG\
{\"msa\"}{\"4_TCOFFEE\"}=\"MSA\";\n$PG{\"msa\"}{\"\
type\"}=\"sequence_multiple_aligner\";\n$PG{\"msa\\
"}{\"ADDRESS\"}=\"https://www.ncbi.nlm.nih.gov/CBB\
research/Schaffer/msa.html\";\n$PG{\"msa\"}{\"lang\
uage\"}=\"C\";\n$PG{\"msa\"}{\"language\"}=\"C\";\\
n$PG{\"msa\"}{\"source\"}=\"ftp://ftp.ncbi.nih.gov\
/pub/msa/msa.tar.Z\";\n$PG{\"msa\"}{\"mode\"}=\"mc\
offee\";\n$PG{\"msa\"}{\"binary\"}=\"msa.pl\";\n$P\
G{\"msa\"}{\"version\"}=\"1.0\";\n$PG{\"msa\"}{\"u\
pdate_action\"}=\"never\";\n$PG{\"dca\"}{\"4_TCOFF\
EE\"}=\"DCA\";\n$PG{\"dca\"}{\"type\"}=\"sequence_\
multiple_aligner\";\n$PG{\"dca\"}{\"ADDRESS\"}=\"h\
ttps://bibiserv2.cebitec.uni-bielefeld.de/dca\";\n\
$PG{\"dca\"}{\"language\"}=\"C\";\n$PG{\"dca\"}{\"\
language\"}=\"C\";\n$PG{\"dca\"}{\"source\"}=\"htt\
ps://bibiserv2.cebitec.uni-bielefeld.de/applicatio\
ns/dca/resources/downloads/dca-1.1-src.tar.gz\";\n\
$PG{\"dca\"}{\"mode\"}=\"mcoffee\";\n$PG{\"dca\"}{\
\"binary\"}=\"dca.pl\";\n$PG{\"dca\"}{\"version\"}\
=\"1.1\";\n$PG{\"dca\"}{\"update_action\"}=\"never\
\";\n$PG{\"muscle\"}{\"4_TCOFFEE\"}=\"MUSCLE\";\n$\
PG{\"muscle\"}{\"type\"}=\"sequence_multiple_align\
er\";\n$PG{\"muscle\"}{\"ADDRESS\"}=\"http://www.d\
rive5.com/muscle/\";\n$PG{\"muscle\"}{\"language\"\
}=\"C++\";\n$PG{\"muscle\"}{\"language2\"}=\"GPP\"\
;\n$PG{\"muscle\"}{\"source\"}=\"http://www.drive5\
.com/muscle/downloads3.7/muscle3.7_src.tar.gz\";\n\
$PG{\"muscle\"}{\"windows\"}=\"http://www.drive5.c\
om/muscle/downloads3.7/muscle3.7_win32.zip\";\n$PG\
{\"muscle\"}{\"linux\"}=\"http://www.drive5.com/mu\
scle/downloads3.7/muscle3.7_linux_ia32.tar.gz\";\n\
$PG{\"muscle\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$P\
G{\"muscle\"}{\"version\"}=\"3.7\";\n$PG{\"pcma\"}\
{\"4_TCOFFEE\"}=\"PCMA\";\n$PG{\"pcma\"}{\"type\"}\
=\"sequence_multiple_aligner\";\n$PG{\"pcma\"}{\"A\
DDRESS\"}=\"http://prodata.swmed.edu/pcma/pcma.php\
\";\n$PG{\"pcma\"}{\"language\"}=\"C\";\n$PG{\"pcm\
a\"}{\"language2\"}=\"C\";\n$PG{\"pcma\"}{\"source\
\"}=\"http://prodata.swmed.edu/download/pub/PCMA/p\
cma.tar.gz\";\n$PG{\"pcma\"}{\"mode\"}=\"mcoffee\"\
;\n$PG{\"pcma\"}{\"version\"}=\"1.0\";\n$PG{\"kali\
gn\"}{\"4_TCOFFEE\"}=\"KALIGN\";\n$PG{\"kalign\"}{\
\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"ka\
lign\"}{\"ADDRESS\"}=\"http://msa.cgb.ki.se\";\n$P\
G{\"kalign\"}{\"language\"}=\"C\";\n$PG{\"kalign\"\
}{\"language2\"}=\"C\";\n$PG{\"kalign\"}{\"source\\
"}=\"http://msa.cgb.ki.se/downloads/kalign/current\
.tar.gz\";\n$PG{\"kalign\"}{\"mode\"}=\"mcoffee\";\
\n$PG{\"kalign\"}{\"version\"}=\"1.0\";\n$PG{\"ama\
p\"}{\"4_TCOFFEE\"}=\"AMAP\";\n$PG{\"amap\"}{\"typ\
e\"}=\"sequence_multiple_aligner\";\n$PG{\"amap\"}\
{\"ADDRESS\"}=\"http://bio.math.berkeley.edu/amap/\
\";\n$PG{\"amap\"}{\"language\"}=\"C++\";\n$PG{\"a\
map\"}{\"language2\"}=\"CXX\";\n$PG{\"amap\"}{\"so\
urce\"}=\"https://github.com/mes5k/amap-align/arch\
ive/amap.zip\";\n$PG{\"amap\"}{\"mode\"}=\"mcoffee\
\";\n$PG{\"amap\"}{\"version\"}=\"2.0\";\n$PG{\"am\
ap\"}{\"update_action\"}=\"never\";\n$PG{\"proda\"\
}{\"4_TCOFFEE\"}=\"PRODA\";\n$PG{\"proda\"}{\"type\
\"}=\"sequence_multiple_aligner\";\n$PG{\"proda\"}\
{\"ADDRESS\"}=\"http://proda.stanford.edu\";\n$PG{\
\"proda\"}{\"language\"}=\"C++\";\n$PG{\"proda\"}{\
\"language2\"}=\"CXX\";\n$PG{\"proda\"}{\"source\"\
}=\"http://proda.stanford.edu/proda_1_0.tar.gz\";\\
n$PG{\"proda\"}{\"mode\"}=\"mcoffee\";\n$PG{\"prod\
a\"}{\"version\"}=\"1.0\";\n$PG{\"prank\"}{\"4_TCO\
FFEE\"}=\"PRANK\";\n$PG{\"prank\"}{\"type\"}=\"seq\
uence_multiple_aligner\";\n$PG{\"prank\"}{\"ADDRES\
S\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/\";\\
n$PG{\"prank\"}{\"language\"}=\"C++\";\n$PG{\"pran\
k\"}{\"language2\"}=\"CXX\";\n$PG{\"prank\"}{\"sou\
rce\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/sr\
c/prank/prank.src.100802.tgz\";\n$PG{\"prank\"}{\"\
mode\"}=\"mcoffee\";\n$PG{\"prank\"}{\"version\"}=\
\"100303\";\n$PG{\"sap\"}{\"4_TCOFFEE\"}=\"SAP\";\\
n$PG{\"sap\"}{\"type\"}=\"structure_pairwise_align\
er\";\n$PG{\"sap\"}{\"ADDRESS\"}=\"https://mathbio\
.crick.ac.uk/wiki/Software#SAP\";\n$PG{\"sap\"}{\"\
language\"}=\"C\";\n$PG{\"sap\"}{\"language2\"}=\"\
C\";\n$PG{\"sap\"}{\"source\"}=\"https://github.co\
m/jkleinj/SAP/archive/v.1.1.3.tar.gz\";\n$PG{\"sap\
\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"sap\"}\
{\"version\"}=\"1.1.3\";\n$PG{\"sap\"}{\"binary\"}\
=\"sap\";\n$PG{\"TMalign\"}{\"4_TCOFFEE\"}=\"TMALI\
GN\";\n$PG{\"TMalign\"}{\"type\"}=\"structure_pair\
wise_aligner\";\n$PG{\"TMalign\"}{\"ADDRESS\"}=\"h\
ttp://zhanglab.ccmb.med.umich.edu/TM-align/TMalign\
.f\";\n$PG{\"TMalign\"}{\"language\"}=\"Fortran\";\
\n$PG{\"TMalign\"}{\"language2\"}=\"Fortran\";\n$P\
G{\"TMalign\"}{\"source\"}=\"http://zhanglab.ccmb.\
med.umich.edu/TM-align/TMalign.f\";\n$PG{\"TMalign\
\"}{\"linux\"}=\"http://zhanglab.ccmb.med.umich.ed\
u/TM-align/TMalign_32.gz\";\n$PG{\"TMalign\"}{\"mo\
de\"}=\"expresso,3dcoffee\";\n$PG{\"TMalign\"}{\"v\
ersion\"}=\"2013.05.11\";\n$PG{\"mustang\"}{\"4_TC\
OFFEE\"}=\"MUSTANG\";\n$PG{\"mustang\"}{\"type\"}=\
\"structure_pairwise_aligner\";\n$PG{\"mustang\"}{\
\"ADDRESS\"}=\"http://lcb.infotech.monash.edu.au/m\
ustang/\";\n$PG{\"mustang\"}{\"language\"}=\"C++\"\
;\n$PG{\"mustang\"}{\"language2\"}=\"CXX\";\n$PG{\\
"mustang\"}{\"source\"}=\"http://lcb.infotech.mona\
sh.edu.au/mustang/mustang_v3.2.3.tgz\";\n$PG{\"mus\
tang\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"mu\
stang\"}{\"version\"}=\"3.2.3\";\n$PG{\"lsqman\"}{\
\"4_TCOFFEE\"}=\"LSQMAN\";\n$PG{\"lsqman\"}{\"type\
\"}=\"structure_pairwise_aligner\";\n$PG{\"lsqman\\
"}{\"ADDRESS\"}=\"empty\";\n$PG{\"lsqman\"}{\"lang\
uage\"}=\"empty\";\n$PG{\"lsqman\"}{\"language2\"}\
=\"empty\";\n$PG{\"lsqman\"}{\"source\"}=\"empty\"\
;\n$PG{\"lsqman\"}{\"update_action\"}=\"never\";\n\
$PG{\"lsqman\"}{\"mode\"}=\"expresso,3dcoffee\";\n\
$PG{\"align_pdb\"}{\"4_TCOFFEE\"}=\"ALIGN_PDB\";\n\
$PG{\"align_pdb\"}{\"type\"}=\"structure_pairwise_\
aligner\";\n$PG{\"align_pdb\"}{\"ADDRESS\"}=\"empt\
y\";\n$PG{\"align_pdb\"}{\"language\"}=\"empty\";\\
n$PG{\"align_pdb\"}{\"language2\"}=\"empty\";\n$PG\
{\"align_pdb\"}{\"source\"}=\"empty\";\n$PG{\"alig\
n_pdb\"}{\"update_action\"}=\"never\";\n$PG{\"alig\
n_pdb\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"f\
ugueali\"}{\"4_TCOFFEE\"}=\"FUGUE\";\n$PG{\"fuguea\
li\"}{\"type\"}=\"structure_pairwise_aligner\";\n$\
PG{\"fugueali\"}{\"ADDRESS\"}=\"http://mizuguchila\
b.org/fugue/\";\n$PG{\"fugueali\"}{\"language\"}=\\
"empty\";\n$PG{\"fugueali\"}{\"language2\"}=\"empt\
y\";\n$PG{\"fugueali\"}{\"source\"}=\"empty\";\n$P\
G{\"fugueali\"}{\"update_action\"}=\"never\";\n$PG\
{\"fugueali\"}{\"mode\"}=\"expresso,3dcoffee\";\n$\
PG{\"dalilite.pl\"}{\"4_TCOFFEE\"}=\"DALILITEc\";\\
n$PG{\"dalilite.pl\"}{\"type\"}=\"structure_pairwi\
se_aligner\";\n$PG{\"dalilite.pl\"}{\"ADDRESS\"}=\\
"built_in\";\n$PG{\"dalilite.pl\"}{\"ADDRESS2\"}=\\
"http://www.ebi.ac.uk/Tools/webservices/services/d\
alilite\";\n$PG{\"dalilite.pl\"}{\"language\"}=\"P\
erl\";\n$PG{\"dalilite.pl\"}{\"language2\"}=\"Perl\
\";\n$PG{\"dalilite.pl\"}{\"source\"}=\"empty\";\n\
$PG{\"dalilite.pl\"}{\"update_action\"}=\"never\";\
\n$PG{\"dalilite.pl\"}{\"mode\"}=\"expresso,3dcoff\
ee\";\n$PG{\"probconsRNA\"}{\"4_TCOFFEE\"}=\"PROBC\
ONSRNA\";\n$PG{\"probconsRNA\"}{\"type\"}=\"RNA_mu\
ltiple_aligner\";\n$PG{\"probconsRNA\"}{\"ADDRESS\\
"}=\"http://probcons.stanford.edu/\";\n$PG{\"probc\
onsRNA\"}{\"language\"}=\"C++\";\n$PG{\"probconsRN\
A\"}{\"language2\"}=\"CXX\";\n$PG{\"probconsRNA\"}\
{\"source\"}=\"http://probcons.stanford.edu/probco\
nsRNA.tar.gz\";\n$PG{\"probconsRNA\"}{\"mode\"}=\"\
mcoffee,rcoffee\";\n$PG{\"probconsRNA\"}{\"version\
\"}=\"1.0\";\n$PG{\"sfold\"}{\"4_TCOFFEE\"}=\"CONS\
AN\";\n$PG{\"sfold\"}{\"type\"}=\"RNA_pairwise_ali\
gner\";\n$PG{\"sfold\"}{\"ADDRESS\"}=\"http://sela\
b.janelia.org/software/consan/\";\n$PG{\"sfold\"}{\
\"language\"}=\"empty\";\n$PG{\"sfold\"}{\"languag\
e2\"}=\"empty\";\n$PG{\"sfold\"}{\"source\"}=\"emp\
ty\";\n$PG{\"sfold\"}{\"update_action\"}=\"never\"\
;\n$PG{\"sfold\"}{\"mode\"}=\"rcoffee\";\n$PG{\"RN\
Aplfold\"}{\"4_TCOFFEE\"}=\"RNAPLFOLD\";\n$PG{\"RN\
Aplfold\"}{\"type\"}=\"RNA_secondarystructure_pred\
ictor\";\n$PG{\"RNAplfold\"}{\"ADDRESS\"}=\"http:/\
/www.tbi.univie.ac.at/RNA/\";\n$PG{\"RNAplfold\"}{\
\"language\"}=\"C\";\n$PG{\"RNAplfold\"}{\"languag\
e2\"}=\"C\";\n$PG{\"RNAplfold\"}{\"source\"}=\"htt\
p://www.tbi.univie.ac.at/RNA/packages/source/Vienn\
aRNA-2.1.9.tar.gz\";\n$PG{\"RNAplfold\"}{\"mode\"}\
=\"rcoffee,\";\n$PG{\"RNAplfold\"}{\"binary\"}=\"R\
NAplfold.tar.gz\";\n$PG{\"RNAplfold\"}{\"version\"\
}=\"2.1.9\";\n$PG{\"retree\"}{\"4_TCOFFEE\"}=\"PHY\
LIP\";\n$PG{\"retree\"}{\"type\"}=\"Phylogeny\";\n\
$PG{\"retree\"}{\"ADDRESS\"}=\"http://evolution.gs\
.washington.edu/phylip/\";\n$PG{\"retree\"}{\"lang\
uage\"}=\"C\";\n$PG{\"retree\"}{\"language2\"}=\"C\
\";\n$PG{\"retree\"}{\"source\"}=\"http://www.tcof\
fee.org/Packages/mirrors/source/phylip-3.66.tar.gz\
\";\n$PG{\"retree\"}{\"mode\"}=\"trmsd,\";\n$PG{\"\
retree\"}{\"binary\"}=\"retree.tar.gz\";\n$PG{\"re\
tree\"}{\"version\"}=\"3.66\";\n$PG{\"hmmtop\"}{\"\
4_TCOFFEE\"}=\"HMMTOP\";\n$PG{\"hmmtop\"}{\"type\"\
}=\"protein_secondarystructure_predictor\";\n$PG{\\
"hmmtop\"}{\"ADDRESS\"}=\"www.enzim.hu/hmmtop/\";\\
n$PG{\"hmmtop\"}{\"language\"}=\"C\";\n$PG{\"hmmto\
p\"}{\"language2\"}=\"C\";\n$PG{\"hmmtop\"}{\"sour\
ce\"}=\"http://www.tcoffee.org/Packages/mirrors/hm\
mtop2.1.tgz\";\n$PG{\"hmmtop\"}{\"binary\"}=\"hmmt\
op\";\n$PG{\"hmmtop\"}{\"update_action\"}=\"never\\
";\n$PG{\"hmmtop\"}{\"mode\"}=\"psicoffee\";\n$PG{\
\"hmmtop\"}{\"version\"}=\"2.1\";\n$PG{\"gorIV\"}{\
\"4_TCOFFEE\"}=\"GOR4\";\n$PG{\"gorIV\"}{\"type\"}\
=\"protein_secondarystructure_predictor\";\n$PG{\"\
gorIV\"}{\"ADDRESS\"}=\"http://mig.jouy.inra.fr/lo\
giciels/gorIV/\";\n$PG{\"gorIV\"}{\"language\"}=\"\
C\";\n$PG{\"gorIV\"}{\"language2\"}=\"C\";\n$PG{\"\
gorIV\"}{\"source\"}=\"http://www.tcoffee.org/Pack\
ages/mirrors/GOR_IV.tar.gz\";\n$PG{\"gorIV\"}{\"up\
date_action\"}=\"never\";\n$PG{\"gorIV\"}{\"mode\"\
}=\"tcoffee\";\n$PG{\"wublast.pl\"}{\"4_TCOFFEE\"}\
=\"EBIWUBLASTc\";\n$PG{\"wublast.pl\"}{\"type\"}=\\
"protein_homology_predictor\";\n$PG{\"wublast.pl\"\
}{\"ADDRESS\"}=\"built_in\";\n$PG{\"wublast.pl\"}{\
\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/webserv\
ices/services/wublast\";\n$PG{\"wublast.pl\"}{\"la\
nguage\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"languag\
e2\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"source\"}=\\
"empty\";\n$PG{\"wublast.pl\"}{\"update_action\"}=\
\"never\";\n$PG{\"wublast.pl\"}{\"mode\"}=\"psicof\
fee,expresso,accurate\";\n$PG{\"blastpgp.pl\"}{\"4\
_TCOFFEE\"}=\"EBIBLASTPGPc\";\n$PG{\"blastpgp.pl\"\
}{\"type\"}=\"protein_homology_predictor\";\n$PG{\\
"blastpgp.pl\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"\
blastpgp.pl\"}{\"ADDRESS2\"}=\"http://www.ebi.ac.u\
k/Tools/webservices/services/blastpgp\";\n$PG{\"bl\
astpgp.pl\"}{\"language\"}=\"Perl\";\n$PG{\"blastp\
gp.pl\"}{\"language2\"}=\"Perl\";\n$PG{\"blastpgp.\
pl\"}{\"source\"}=\"empty\";\n$PG{\"blastpgp.pl\"}\
{\"update_action\"}=\"never\";\n$PG{\"blastpgp.pl\\
"}{\"mode\"}=\"psicoffee,expresso,accurate\";\n$PG\
{\"blastall\"}{\"4_TCOFFEE\"}=\"blastall\";\n$PG{\\
"blastall\"}{\"type\"}=\"protein_homology_predicto\
r\";\n$PG{\"blastall\"}{\"ADDRESS\"}=\"ftp://ftp.n\
cbi.nih.gov/blast/executables/LATEST\";\n$PG{\"bla\
stall\"}{\"language\"}=\"C\";\n$PG{\"blastall\"}{\\
"language2\"}=\"C\";\n$PG{\"blastall\"}{\"source\"\
}=\"ftp://ftp.ncbi.nlm.nih.gov/blast/executables/b\
last+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz\";\n$PG{\\
"blastall\"}{\"update_action\"}=\"never\";\n$PG{\"\
blastall\"}{\"mode\"}=\"psicoffee,expresso,3dcoffe\
e\";\n$PG{\"legacy_blast.pl\"}{\"4_TCOFFEE\"}=\"NC\
BIBLAST\";\n$PG{\"legacy_blast.pl\"}{\"type\"}=\"p\
rotein_homology_predictor\";\n$PG{\"legacy_blast.p\
l\"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.gov/blast/e\
xecutables/LATEST\";\n$PG{\"legacy_blast.pl\"}{\"l\
anguage\"}=\"C\";\n$PG{\"legacy_blast.pl\"}{\"lang\
uage2\"}=\"C\";\n$PG{\"legacy_blast.pl\"}{\"source\
\"}=\"ftp://ftp.ncbi.nlm.nih.gov/blast/executables\
/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz\";\n$PG\
{\"legacy_blast.pl\"}{\"update_action\"}=\"never\"\
;\n$PG{\"legacy_blast.pl\"}{\"mode\"}=\"psicoffee,\
expresso,3dcoffee\";\n$PG{\"SOAP::Lite\"}{\"4_TCOF\
FEE\"}=\"SOAPLITE\";\n$PG{\"SOAP::Lite\"}{\"type\"\
}=\"library\";\n$PG{\"SOAP::Lite\"}{\"ADDRESS\"}=\\
"http://cpansearch.perl.org/src/MKUTTER/SOAP-Lite-\
0.710.08/Makefile.PL\";\n$PG{\"SOAP::Lite\"}{\"lan\
guage\"}=\"Perl\";\n$PG{\"SOAP::Lite\"}{\"language\
2\"}=\"Perl\";\n$PG{\"SOAP::Lite\"}{\"source\"}=\"\
empty\";\n$PG{\"SOAP::Lite\"}{\"update_action\"}=\\
"never\";\n$PG{\"SOAP::Lite\"}{\"mode\"}=\"none\";\
\n$PG{\"XML::Simple\"}{\"4_TCOFFEE\"}=\"XMLSIMPLE\\
";\n$PG{\"XML::Simple\"}{\"type\"}=\"library\";\n$\
PG{\"XML::Simple\"}{\"ADDRESS\"}=\"http://search.c\
pan.org/~grantm/XML-Simple-2.18/lib/XML/Simple.pm\\
";\n$PG{\"XML::Simple\"}{\"language\"}=\"Perl\";\n\
$PG{\"XML::Simple\"}{\"language2\"}=\"Perl\";\n$PG\
{\"XML::Simple\"}{\"source\"}=\"empty\";\n$PG{\"XM\
L::Simple\"}{\"mode\"}=\"psicoffee,expresso,accura\
te\";\n$PG{\"x3dna\"}{\"4_TCOFFEE\"}=\"x3dna\";\n$\
PG{\"x3dna\"}{\"type\"}=\"RNA_secondarystructure_p\
redictor\";\n$PG{\"x3dna\"}{\"ADDRESS\"}=\"http://\
x3dna.bio.columbia.edu/\";\n$PG{\"x3dna\"}{\"sourc\
e\"}=\"http://www.tcoffee.org/Packages/mirrors/sou\
rce/x3dna-v2.3-linux-64bit.tar.gz\";\n$PG{\"x3dna\\
"}{\"mode\"}=\"saracoffee\";\n$PG{\"x3dna\"}{\"upd\
ate_action\"}=\"never\";\n$PG{\"fsa\"}{\"4_TCOFFEE\
\"}=\"FSA\";\n$PG{\"fsa\"}{\"type\"}=\"sequence_mu\
ltiple_aligner\";\n$PG{\"fsa\"}{\"ADDRESS\"}=\"htt\
p://fsa.sourceforge.net/\";\n$PG{\"fsa\"}{\"langua\
ge\"}=\"C++\";\n$PG{\"fsa\"}{\"language2\"}=\"CXX\\
";\n$PG{\"fsa\"}{\"source\"}=\"http://sourceforge.\
net/projects/fsa/files/fsa-1.15.3.tar.gz/download/\
\";\n$PG{\"fsa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"fs\
a\"}{\"version\"}=\"1.15.3\";\n$PG{\"fsa\"}{\"upda\
te_action\"}=\"never\";\n$PG{\"mus4\"}{\"4_TCOFFEE\
\"}=\"MUS4\";\n$PG{\"mus4\"}{\"type\"}=\"sequence_\
multiple_aligner\";\n$PG{\"mus4\"}{\"ADDRESS\"}=\"\
http://www.drive5.com/muscle/\";\n$PG{\"mus4\"}{\"\
language\"}=\"C++\";\n$PG{\"mus4\"}{\"language2\"}\
=\"GPP\";\n$PG{\"mus4\"}{\"source\"}=\"http://www.\
drive5.com/muscle/muscle4.0_src.tar.gz\";\n$PG{\"m\
us4\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"mus4\\
"}{\"version\"}=\"4.0\";\n$PG{\"mus4\"}{\"update_a\
ction\"}=\"never\";\n$MODE{\"tcoffee\"}{\"name\"}=\
\"tcoffee\";\n$MODE{\"rcoffee\"}{\"name\"}=\"rcoff\
ee\";\n$MODE{\"3dcoffee\"}{\"name\"}=\"3dcoffee\";\
\n$MODE{\"mcoffee\"}{\"name\"}=\"mcoffee\";\n$MODE\
{\"expresso\"}{\"name\"}=\"expresso\";\n$MODE{\"tr\
msd\"}{\"name\"}=\"trmsd\";\n$MODE{\"accurate\"}{\\
"name\"}=\"accurate\";\n$MODE{\"seq_reformat\"}{\"\
name\"}=\"seq_reformat\";\n\n\n$PG{C}{compiler}=\"\
gcc\";\n$PG{C}{compiler_flag}=\"CC\";\n$PG{C}{opti\
ons}=\"\";\n$PG{C}{options_flag}=\"CFLAGS\";\n$PG{\
C}{type}=\"compiler\";\n\n$PG{\"CXX\"}{compiler}=\\
"g++\";\n$PG{\"CXX\"}{compiler_flag}=\"CXX\";\n$PG\
{\"CXX\"}{options}=\"\";\n$PG{\"CXX\"}{options_fla\
g}=\"CXXFLAGS\";\n$PG{CXX}{type}=\"compiler\";\n\n\
$PG{\"CPP\"}{compiler}=\"g++\";\n$PG{\"CPP\"}{comp\
iler_flag}=\"CPP\";\n$PG{\"CPP\"}{options}=\"\";\n\
$PG{\"CPP\"}{options_flag}=\"CPPFLAGS\";\n$PG{CPP}\
{type}=\"compiler\";\n\n$PG{\"GPP\"}{compiler}=\"g\
++\";\n$PG{\"GPP\"}{compiler_flag}=\"GPP\";\n$PG{\\
"GPP\"}{options}=\"\";\n$PG{\"GPP\"}{options_flag}\
=\"CFLAGS\";\n$PG{GPP}{type}=\"compiler\";\n\n$PG{\
Fortran}{compiler}=\"g77\";\n$PG{Fortran}{compiler\
_flag}=\"FCC\";\n$PG{Fortran}{type}=\"compiler\";\\
n\n$PG{Perl}{compiler}=\"CPAN\";\n$PG{Perl}{type}=\
\"compiler\";\n\n$SUPPORTED_OS{macosx}=\"Macintosh\
\";\n$SUPPORTED_OS{linux}=\"Linux\";\n$SUPPORTED_O\
S{windows}=\"Cygwin\";\n\n\n\n$MODE{t_coffee}{desc\
ription}=\" for regular multiple sequence alignmen\
ts\";\n$MODE{rcoffee} {description}=\" for RNA mul\
tiple sequence alignments\";\n\n$MODE{psicoffee} {\
description}=\" for Homology Extended multiple seq\
uence alignments\";\n$MODE{expresso}{description}=\
\" for very accurate structure based multiple sequ\
ence alignments\";\n$MODE{\"3dcoffee\"}{descriptio\
n}=\" for multiple structure alignments\";\n$MODE{\
mcoffee} {description}=\" for combining alternativ\
e multiple sequence alignment packages\\n------- i\
nto a unique meta-package. The installer will uplo\
ad several MSA packages and compile them\\n\n\";\n\
\n\n&post_process_PG();\nreturn;\n}\n\nsub post_pr\
ocess_PG\n  {\n    my $p;\n    \n    %PG=&name2dna\
me (%PG);\n    %MODE=&name2dname(%MODE);\n    fore\
ach $p (keys(%PG)){if ( $PG{$p}{type} eq \"compile\
r\"){$PG{$p}{update_action}=\"never\";}}\n    \n  \
}\n\nsub name2dname\n  {\n    my (%L)=(@_);\n    m\
y ($l, $ml);\n    \n    foreach my $pg (keys(%L))\\
n      {\n	$l=length ($pg);\n	if ( $l>$ml){$ml=$l;\
}\n      }\n    $ml+=1;\n    foreach my $pg (keys(\
%L))\n      {\n	my $name;\n	$l=$ml-length ($pg);\n\
	$name=$pg;\n	for ( $b=0; $b<$l; $b++)\n	  {\n	   \
 $name .=\" \";\n	  }\n	$L{$pg}{dname}=$name;\n   \
   }\n    return %L;\n  }\n\nsub env_file2putenv\n\
  {\n    my $f=@_[0];\n    my $F=new FileHandle;\n\
    my $n;\n    \n    open ($F, \"$f\");\n    whil\
e (<$F>)\n      {\n	my $line=$_;\n	my($var, $value\
)=($_=~/(\\S+)\\=(\\S*)/);\n	$ENV{$var}=$value;\n	\
$ENV_SET{$var}=1;\n	$n++;\n      }\n    close ($F)\
;\n    return $n;\n  }\n\nsub replace_line_in_file\
\n  {\n    my ($file, $wordin, $wordout)=@_;\n    \
my $O=new FileHandle;\n    my $I=new FileHandle;\n\
    my $l;\n    if (!-e $file){return;}\n    \n   \
 system (\"mv $file $file.old\");\n    open ($O, \\
">$file\");\n    open ($I, \"$file.old\");\n    wh\
ile (<$I>)\n      {\n	$l=$_;\n	if (!($l=~/$wordin/\
)){print $O \"$l\";}\n	elsif ( $wordout ne \"\"){$\
l=~s/$wordin/$wordout/g;print $O \"$l\";}\n      }\
\n    close ($O);\n    close ($I);\n    return;\n \
 }\n\nsub add_C_libraries\n  {\n   my ($file,$firs\
t,@list)=@_;\n   \n    my $O=new FileHandle;\n    \
my $I=new FileHandle;\n    my ($l,$anchor);\n    i\
f (!-e $file){return;}\n   \n    $anchor=\"#includ\
e <$first>\";\n	 \n    system (\"mv $file $file.ol\
d\");\n    open ($O, \">$file\");\n    open ($I, \\
"$file.old\");\n    while (<$I>)\n      {\n	$l=$_;\
\n	print $O \"$l\";\n	if (!($l=~/$anchor/))\n	   {\
\n	    \n	    foreach my $lib (@list)\n	       {\n\
                  print $O \"#include <$lib>\\n\";\
\n	       }\n           }\n      }\n    close ($O)\
;\n    close ($I);\n    return;\n    }\n","use Env\
;\nuse Cwd;\n@suffix=(\"tmp\", \"temp\", \"cache\"\
, \"t_coffee\", \"core\", \"tcoffee\");\n\nif ($#A\
RGV==-1)\n  {\n    print \"clean_cache.pl -file <f\
ile to add in -dir> -dir=<dir> -size=<value in Mb>\
\\n0: unlimited -1 always.\\nWill only clean direc\
tories matching:[\";\n    foreach $k(@suffix){prin\
t \"*$k* \";}\n    print \"]\\n\";\n    exit (EXIT\
_FAILURE);\n  }\n\n$cl=join (\" \",@ARGV);\nif (($\
cl=~/\\-no_action/))\n  {\n    exit (EXIT_SUCCESS)\
;\n  }\n\nif (($cl=~/\\-debug/))\n  {\n    $DEBUG=\
1;\n  }\nelse\n  {\n    $DEBUG=0;\n  }\n\nif (($cl\
=~/\\-dir=(\\S+)/))\n  {\n    $dir=$1;\n  }\nelse\\
n  {\n    $dir=\"./\";\n  }\n\nif ($cl=~/\\-file=(\
\\S+)/)\n  {\n    $file=$1;\n  }\nelse\n  {\n    $\
file=0;\n  }\n\nif ($cl=~/\\-size=(\\S+)/)\n  {\n \
   $max_size=$1;\n  }\nelse\n  {\n    $max_size=0;\
#unlimited\n  }\nif ($cl=~/\\-force/)\n  {\n    $f\
orce=1;\n  }\nelse\n  {\n    $force=0;\n  }\n\nif \
($cl=~/\\-age=(\\S+)/)\n  {\n    $max_age=$1;\n  }\
\nelse\n  {\n    $max_age=0;#unlimited\n  }\n\n$ma\
x_size*=1000000;\nif ( ! -d $dir)\n  {\n    print \
STDERR \"\\nCannot process $dir: does not exist \\\
n\";\n    exit (EXIT_FAILURE);\n  }\n\nif ( !($dir\
=~/^\\//))\n  {\n    $base=cwd();\n    $dir=\"$bas\
e/$dir\";\n  }\n\n$proceed=0;\nforeach $s (@suffix\
)\n  {\n    \n    if (($dir=~/$s/)){$proceed=1;}\n\
    $s=uc ($s);\n    if (($dir=~/$s/)){$proceed=1;\
}\n  }\nif ( $proceed==0)\n  {\n    print STDERR \\
"Clean_cache.pl can only clean directories whose a\
bsolute path name contains the following strings:\\
";\n    foreach $w (@suffix) {print STDERR \"$w \"\
;$w=lc($w); print STDERR \"$w \";}\n    print STDE\
RR \"\\nCannot process $dir\\n\";\n    exit (EXIT_\
FAILURE);\n  }\n\n$name_file=\"$dir/name_file.txt\\
";\n$size_file=\"$dir/size_file.txt\";\nif ( $forc\
e){&create_ref_file ($dir,$name_file,$size_file);}\
\nif ($file){&add_file ($dir, $name_file, $size_fi\
le, $file);}\n&clean_dir ($dir, $name_file, $size_\
file, $max_size,$max_age);\nexit (EXIT_SUCCESS);\n\
\nsub clean_dir \n  {\n    my ($dir, $name_file, $\
size_file, $max_size, $max_age)=@_;\n    my ($tot_\
size, $size, $f, $s);\n\n  \n    $tot_size=&get_to\
t_size ($dir, $name_file, $size_file);\n\n    if (\
 $tot_size<=$max_size){return ;}\n    else {$max_s\
ize/=2;}\n    \n    #recreate the name file in cas\
e some temprary files have not been properly regis\
tered\n    &create_ref_file ($dir, $name_file, $si\
ze_file, $max_age);\n  \n    $new_name_file=&vtmpn\
am();\n    open (R, \"$name_file\");\n    open (W,\
 \">$new_name_file\");\n    while (<R>)\n      {\n\
	my $line=$_;\n	\n	($f, $s)=($line=~/(\\S+) (\\S+)\
/);\n	if ( !($f=~/\\S/)){next;}\n	\n	elsif ($max_s\
ize && $tot_size>=$max_size && !($f=~/name_file/))\
\n	  {\n	    remove ( \"$dir/$f\");\n	    $tot_siz\
e-=$s;\n	  }\n	elsif ( $max_age && -M(\"$dir/$f\")\
>=$max_age)\n	  {\n	    remove ( \"$dir/$f\");\n	 \
   $tot_size-=$s;\n	  }\n	else\n	  {\n	    print W\
 \"$f $s\\n\";\n	  }\n      }\n    close (R);\n   \
 close (W);\n    open (F, \">$size_file\");\n    p\
rint F \"$tot_size\";\n    if ( -e $new_name_file)\
{`mv $new_name_file $name_file`;}\n    close (F);\\
n  }\nsub get_tot_size\n  {\n    my ($dir, $name_f\
ile, $size_file)=@_;\n    my $size;\n    \n    if \
( !-d $dir){return 0;}\n    if ( !-e $name_file)\n\
      {\n	\n	&create_ref_file ($dir, $name_file, $\
size_file);\n      }\n    open (F, \"$size_file\")\
;\n    $size=<F>;\n    close (F);\n    chomp ($siz\
e);\n    return $size;\n  }\nsub size \n  {\n    m\
y $f=@_[0];\n\n    if ( !-d $f){return -s($f);}\n \
   else {return &dir2size($f);}\n  }\nsub dir2size\
\n  {\n    my $d=@_[0];\n    my ($s, $f);\n    \n \
   if ( !-d $d) {return 0;}\n    \n    foreach $f \
(&dir2list ($d))\n      {\n	if ( -d $f){$s+=&dir2s\
ize (\"$d/$f\");}\n	else {$s+= -s \"$dir/$f\";}\n \
     }\n    return $s;\n  }\n\nsub remove \n  {\n \
   my $file=@_[0];\n    my ($f);\n    \n    debug_\
print( \"--- $file ---\\n\");\n    if (($file eq \\
".\") || ($file eq \"..\") || ($file=~/\\*/)){retu\
rn EXIT_FAILURE;}\n    elsif ( !-d $file)\n      {\
\n	debug_print (\"unlink $file\\n\");\n	if (-e $fi\
le){unlink ($file);}\n      }\n    elsif ( -d $fil\
e)\n      {\n	debug_print (\"++++++++ $file ++++++\
+\\n\");\n	foreach $f (&dir2list($file))\n	  {\n	 \
   &remove (\"$file/$f\");\n	  }\n	debug_print (\"\
rmdir $file\\n\");\n	rmdir $file;\n      }\n    el\
se\n      {\n	debug_print (\"????????? $file ?????\
???\\n\");\n      }\n    return EXIT_SUCCESS;\n  }\
\n\nsub dir2list\n  {\n    my $dir=@_[0];\n    my \
(@list1, @list2,@list3, $l);\n\n    opendir (DIR,$\
dir);\n    @list1=readdir (DIR);\n    closedir (DI\
R);\n    \n    foreach $l (@list1)\n      {\n	if (\
 $l ne \".\" && $l ne \"..\"){@list2=(@list2, $l);\
}\n      }\n    @list3 = sort { (-M \"$dir/$list2[\
$b]\") <=> (-M \"$dir/$list2[$a]\")} @list2;\n    \
return @list3;\n    \n  }\n\nsub debug_print\n  {\\
n    \n    if ($DEBUG==1){print @_;}\n    \n  }\ns\
ub create_ref_file\n  {\n    my ($dir,$name_file,$\
size_file)=@_;\n    my ($f, $s, $tot_size, @l);\n \
   \n    if ( !-d $dir){return;}\n    \n    @l=&di\
r2list ($dir);\n    open (F, \">$name_file\");\n  \
  foreach $f (@l)\n      {\n	$s=&size(\"$dir/$f\")\
;\n	$tot_size+=$s;\n	print F \"$f $s\\n\";\n      \
}\n    &myecho ($tot_size, \">$size_file\");\n    \
close (F);\n  }\nsub add_file \n  {\n    my ($dir,\
$name_file,$size_file,$file)=@_;\n    my ($s, $tot\
_size);\n    \n    if ( !-d $dir)   {return;}\n   \
 if ( !-e \"$dir/$file\" ) {return;}\n    if ( !-e\
 $name_file){&create_ref_file ($dir,$name_file,$si\
ze_file);}\n					    \n    $s=&size(\"$dir/$file\"\
);\n    open (F, \">>$name_file\");\n    print F \\
"$file\\n\";\n    close (F);\n\n    $tot_size=&get\
_tot_size ($dir,$name_file,$size_file);\n    $tot_\
size+=$s;\n    &myecho ($tot_size, \">$size_file\"\
);\n    \n  }\n	\nsub myecho\n  {\n    my ($string\
, $file)=@_;\n    open (ECHO, $file) || die;\n    \
print ECHO \"$string\";\n    close (ECHO);\n  }\n \
   \n		\n	\nsub vtmpnam\n  {\n    my $tmp_file_nam\
e;\n    $tmp_name_counter++;\n    $tmp_file_name=\\
"tmp_file_for_clean_cache_pdb$$.$tmp_name_counter\\
";\n    $tmp_file_list[$ntmp_file++]=$tmp_file_nam\
e;\n    if ( -e $tmp_file_name) {return &vtmpnam (\
);}\n    else {return $tmp_file_name;}\n  }\n","\n\
my $address=\"http://www.tcoffee.org/Data/Datasets\
/NatureProtocolsDataset.tar.gz\";\nmy $out=\"Natur\
eProtocolsDataset.tar.gz\";\n&url2file ($address,$\
out);\n\nif ( -e $out)\n  {\n    \n    system (\"g\
unzip NatureProtocolsDataset.tar.gz\");\n    syste\
m (\"tar -xvf NatureProtocolsDataset.tar\");\n  	s\
ystem (\"rm -rf NatureProtocolsDataset.tar\");  \n\
    print \"Your Data Set is in the Folder 'Nature\
ProtocolsDataset'\\n\";\n  }\nelse \n  {\n    prin\
t \"Could not Download Dataset --- Web site may be\
 down -- Try again later\\n\";\n  }\n\n\n\n\nsub u\
rl2file\n{\n    my ($address, $out, $wget_arg, $cu\
rl_arg)=(@_);\n    my ($pg, $flag, $r, $arg, $coun\
t);\n    \n    if (!$CONFIGURATION){&check_configu\
ration (\"wget\", \"INTERNET\", \"gzip\");$CONFIGU\
RATION=1;}\n    \n    if (&pg_is_installed (\"wget\
\"))   {$pg=\"wget\"; $flag=\"-O\";$arg=$wget_arg;\
}\n    elsif (&pg_is_installed (\"curl\")){$pg=\"c\
url\"; $flag=\"-o\";$arg=$curl_arg;}\n    return s\
ystem (\"$pg $address $flag $out>/dev/null 2>/dev/\
null\");\n\n}\n\nsub pg_is_installed\n  {\n    my \
@ml=@_;\n    my $r, $p, $m;\n    my $supported=0;\\
n    \n    my $p=shift (@ml);\n    if ($p=~/::/)\n\
      {\n	if (system (\"perl -M$p -e 1\")==$EXIT_S\
UCCESS){return 1;}\n	else {return 0;}\n      }\n  \
  else\n      {\n	$r=`which $p 2>/dev/null`;\n	if \
($r eq \"\"){return 0;}\n	else {return 1;}\n      \
}\n  }\nsub check_configuration \n    {\n      my \
@l=@_;\n      my $v;\n      foreach my $p (@l)\n	{\
\n	  \n	  if   ( $p eq \"EMAIL\")\n	    { \n	     \
 if ( !($EMAIL=~/@/))\n		{\n		  exit (EXIT_FAILURE\
);\n		}\n	    }\n	  elsif( $p eq \"INTERNET\")\n	 \
   {\n	      if ( !&check_internet_connection())\n\
		{\n		  exit (EXIT_FAILURE);\n		}\n	    }\n	  els\
if( $p eq \"wget\")\n	    {\n	      if (!&pg_is_in\
stalled (\"wget\") && !&pg_is_installed (\"curl\")\
)\n		{\n		  exit (EXIT_FAILURE);\n		}\n	    }\n	  \
elsif( !(&pg_is_installed ($p)))\n	    {\n	      e\
xit (EXIT_FAILURE);\n	    }\n	}\n      return 1;\n\
    }\nsub check_internet_connection\n  {\n    my \
$internet;\n    my $tmp;\n    &check_configuration\
 ( \"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \\
n    if     (&pg_is_installed    (\"wget\")){`wget\
 www.google.com -O$tmp >/dev/null 2>/dev/null`;}\n\
    elsif  (&pg_is_installed    (\"curl\")){`curl \
www.google.com -o$tmp >/dev/null 2>/dev/null`;}\n \
   \n    if ( !-e $tmp || -s $tmp < 10){$internet=\
0;}\n    else {$internet=1;}\n    if (-e $tmp){unl\
ink $tmp;}\n\n    return $internet;\n  }\n\nsub vt\
mpnam\n      {\n	my $r=rand(100000);\n	my $f=\"fil\
e.$r.$$\";\n	while (-e $f)\n	  {\n	    $f=vtmpnam(\
);\n	  }\n	push (@TMPFILE_LIST, $f);\n	return $f;\\
n      }\n\n","\n$t_coffee=\"t_coffee\";\n\nforeac\
h $value ( @ARGV)\n  {\n    $seq_file=$seq_file.\"\
 \".$value;\n  }\n\n$name=$ARGV[0];\n$name=~s/\\.[\
^\\.]*$//;\n$lib_name=\"$name.mocca_lib\";\n$type=\
`t_coffee $seq_file -get_type -quiet`;\nchop ($typ\
e);\n\nif ( $type eq \"PROTEIN\"){$lib_mode=\"lali\
gn_rs_s_pair -lalign_n_top 20\";}\nelsif ( $type e\
q\"DNA\"){$lib_mode=\"lalign_rs_s_dna_pair -lalign\
_n_top 40\";}\n\nif ( !(-e $lib_name))\n  {\n	  \n\
  $command=\"$t_coffee -mocca -seq_weight=no -cosm\
etic_penalty=0 -mocca_interactive -in $lib_mode -o\
ut_lib $lib_name -infile $seq_file\";\n  \n  }\nel\
sif ( (-e $lib_name))\n  {\n  $command=\"$t_coffee\
 -mocca -seq_weight=no -cosmetic_penalty=0 -mocca_\
interactive -in $lib_name -infile $seq_file\";\n  \
\n  }\n\nsystem ($command);\n\nexit;\n\n","my $WSD\
L = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/W\
SDaliLite.wsdl';\n\nuse SOAP::Lite;\nuse Data::Dum\
per;\nuse Getopt::Long qw(:config no_ignore_case b\
undling);\nuse File::Basename;\n\nmy $checkInterva\
l = 5;\n\nmy %params=(\n	    'async' => '1', # Use\
 async mode and simulate sync mode in client\n	   \
 );\nGetOptions(\n    'pdb1=s'     => \\$params{'s\
equence1'},\n    'chainid1=s' => \\$params{'chaini\
d1'},\n    'pdb2=s'     => \\$params{'sequence2'},\
\n    'chainid2=s' => \\$params{'chainid2'},\n    \
\"help|h\"	 => \\$help, # Usage info\n    \"async|\
a\"	 => \\$async, # Asynchronous submission\n    \\
"polljob\"	 => \\$polljob, # Get results\n    \"st\
atus\"	 => \\$status, # Get status\n    \"jobid|j=\
s\"  => \\$jobid, # JobId\n    \"email|S=s\"  => \\
\$params{email}, # E-mail address\n    \"trace\"  \
    => \\$trace, # SOAP messages\n    \"sequence=s\
\" => \\$sequence, # Input PDB\n    );\n\nmy $scri\
ptName = basename($0, ());\nif($help) {\n    &usag\
e();\n    exit(0);\n}\n\nif($trace) {\n    print \\
"Tracing active\\n\";\n    SOAP::Lite->import(+tra\
ce => 'debug');\n}\n\nmy $soap = SOAP::Lite\n    -\
>service($WSDL)\n    ->on_fault(sub {\n        my \
$soap = shift;\n        my $res = shift;\n        \
# Throw an exception for all faults\n        if(re\
f($res) eq '') {\n            die($res);\n        \
} else {\n            die($res->faultstring);\n   \
     }\n        return new SOAP::SOM;\n    }\n    \
           );\n\nif( !($polljob || $status) &&\n  \
  !( defined($params{'sequence1'}) && defined($par\
ams{'sequence2'}) )\n    ) {\n    print STDERR 'Er\
ror: bad option combination', \"\\n\";\n    &usage\
();\n    exit(1);\n}\nelsif($polljob && defined($j\
obid)) {\n    print \"Getting results for job $job\
id\\n\";\n    getResults($jobid);\n}\nelsif($statu\
s && defined($jobid)) {\n    print STDERR \"Gettin\
g status for job $jobid\\n\";\n    my $result = $s\
oap->checkStatus($jobid);\n    print STDOUT \"$res\
ult\", \"\\n\";\n    if($result eq 'DONE') {\n	pri\
nt STDERR \"To get results: $scriptName --polljob \
--jobid $jobid\\n\";\n    }\n}\nelse {\n    if(-f \
$params{'sequence1'}) {\n	$params{'sequence1'} = r\
ead_file($params{'sequence1'});\n    }\n    if(-f \
$params{'sequence2'}) {\n	$params{'sequence2'} = r\
ead_file($params{'sequence2'});\n    }\n\n    my $\
jobid;\n    my $paramsData = SOAP::Data->name('par\
ams')->type(map=>\\%params);\n    # For SOAP::Lite\
 0.60 and earlier parameters are passed directly\n\
    if($SOAP::Lite::VERSION eq '0.60' || $SOAP::Li\
te::VERSION =~ /0\\.[1-5]/) {\n        $jobid = $s\
oap->runDaliLite($paramsData);\n    }\n    # For S\
OAP::Lite 0.69 and later parameter handling is dif\
ferent, so pass\n    # undef's for templated param\
s, and then pass the formatted args.\n    else {\n\
        $jobid = $soap->runDaliLite(undef,\n				  \
   $paramsData);\n    }\n\n    if (defined($async)\
) {\n	print STDOUT $jobid, \"\\n\";\n        print\
 STDERR \"To check status: $scriptName --status --\
jobid $jobid\\n\";\n    } else { # Synchronous mod\
e\n        print STDERR \"JobId: $jobid\\n\";\n   \
     sleep 1;\n        getResults($jobid);\n    }\\
n}\n\nsub clientPoll($) {\n    my $jobid = shift;\\
n    my $result = 'PENDING';\n    # Check status a\
nd wait if not finished\n    #print STDERR \"Check\
ing status: $jobid\\n\";\n    while($result eq 'RU\
NNING' || $result eq 'PENDING') {\n        $result\
 = $soap->checkStatus($jobid);\n        print STDE\
RR \"$result\\n\";\n        if($result eq 'RUNNING\
' || $result eq 'PENDING') {\n            # Wait b\
efore polling again.\n            sleep $checkInte\
rval;\n        }\n    }\n}\n\nsub getResults($) {\\
n    $jobid = shift;\n    # Check status, and wait\
 if not finished\n    clientPoll($jobid);\n    # U\
se JobId if output file name is not defined\n    u\
nless(defined($outfile)) {\n        $outfile=$jobi\
d;\n    }\n    # Get list of data types\n    my $r\
esultTypes = $soap->getResults($jobid);\n    # Get\
 the data and write it to a file\n    if(defined($\
outformat)) { # Specified data type\n        my $s\
elResultType;\n        foreach my $resultType (@$r\
esultTypes) {\n            if($resultType->{type} \
eq $outformat) {\n                $selResultType =\
 $resultType;\n            }\n        }\n        $\
res=$soap->poll($jobid, $selResultType->{type});\n\
        write_file($outfile.'.'.$selResultType->{e\
xt}, $res);\n    } else { # Data types available\n\
        # Write a file for each output type\n     \
   for my $resultType (@$resultTypes){\n          \
  #print \"Getting $resultType->{type}\\n\";\n    \
        $res=$soap->poll($jobid, $resultType->{typ\
e});\n            write_file($outfile.'.'.$resultT\
ype->{ext}, $res);\n        }\n    }\n}\n\nsub rea\
d_file($) {\n    my $filename = shift;\n    open(F\
ILE, $filename);\n    my $content;\n    my $buffer\
;\n    while(sysread(FILE, $buffer, 1024)) {\n	$co\
ntent.= $buffer;\n    }\n    close(FILE);\n    ret\
urn $content;\n}\n\nsub write_file($$) {\n    my (\
$tmp,$entity) = @_;\n    print STDERR \"Creating r\
esult file: \".$tmp.\"\\n\";\n    unless(open (FIL\
E, \">$tmp\")) {\n	return 0;\n    }\n    syswrite(\
FILE, $entity);\n    close (FILE);\n    return 1;\\
n}\n\nsub usage {\n    print STDERR <<EOF\nDaliLit\
e\n========\n\nPairwise comparison of protein stru\
ctures\n\n[Required]\n\n  --pdb1                : \
str  : PDB ID for structure 1\n  --pdb2           \
     : str  : PDB ID for structure 2\n\n[Optional]\
\n\n  --chain1              : str  : Chain identif\
er in structure 1\n  --chain2              : str  \
: Chain identifer in structure 2\n\n[General]\n\n \
 -h, --help            :      : prints this help t\
ext\n  -S, --email           : str  : user email a\
ddress\n  -a, --async           :      : asynchron\
ous submission\n      --status          :      : p\
oll for the status of a job\n      --polljob      \
   :      : poll for the results of a job\n  -j, -\
-jobid           : str  : jobid for an asynchronou\
s job\n  -O, --outfile         : str  : file name \
for results (default is jobid)\n      --trace	    \
    :      : show SOAP messages being interchanged\
 \n\nSynchronous job:\n\n  The results/errors are \
returned as soon as the job is finished.\n  Usage:\
 $scriptName --email <your\\@email> [options] pdbF\
ile [--outfile string]\n  Returns: saves the resul\
ts to disk\n\nAsynchronous job:\n\n  Use this if y\
ou want to retrieve the results at a later time. T\
he results \n  are stored for up to 24 hours. \n  \
The asynchronous submission mode is recommended wh\
en users are submitting \n  batch jobs or large da\
tabase searches	\n  Usage: $scriptName --email <yo\
ur\\@email> --async [options] pdbFile\n  Returns: \
jobid\n\n  Use the jobid to query for the status o\
f the job. \n  Usage: $scriptName --status --jobid\
 <jobId>\n  Returns: string indicating the status \
of the job:\n    DONE - job has finished\n    RUNN\
ING - job is running\n    NOT_FOUND - job cannot b\
e found\n    ERROR - the jobs has encountered an e\
rror\n\n  When done, use the jobid to retrieve the\
 status of the job. \n  Usage: $scriptName --pollj\
ob --jobid <jobId> [--outfile string]\n\n[Help]\n\\
n  For more detailed help information refer to\n  \
http://www.ebi.ac.uk/DaliLite/\nEOF\n;\n}\n","my $\
WSDL = 'http://www.ebi.ac.uk/Tools/webservices/wsd\
l/WSWUBlast.wsdl';\n\nuse strict;\nuse SOAP::Lite;\
\nuse Getopt::Long qw(:config no_ignore_case bundl\
ing);\nuse File::Basename;\n\nmy $checkInterval = \
15;\n\nmy $numOpts = scalar(@ARGV);\nmy ($outfile,\
 $outformat, $help, $async, $polljob, $status, $id\
s, $jobid, $trace, $sequence);\nmy %params= ( # De\
faults\n	      'async' => 1, # Force into async mo\
de\n	      'exp' => 10.0, # E-value threshold\n	  \
    'numal' => 50, # Maximum number of alignments\\
n	      'scores' => 100, # Maximum number of score\
s\n            );\nGetOptions( # Map the options i\
nto variables\n    \"program|p=s\"     => \\$param\
s{program}, # BLAST program\n    \"database|D=s\" \
   => \\$params{database}, # Search database\n    \
\"matrix|m=s\"      => \\$params{matrix}, # Scorin\
g matrix\n    \"exp|E=f\"         => \\$params{exp\
}, # E-value threshold\n    \"echofilter|e\"    =>\
 \\$params{echofilter}, # Display filtered sequenc\
e\n    \"filter|f=s\"      => \\$params{filter}, #\
 Low complexity filter name\n    \"alignments|b=i\\
"  => \\$params{numal}, # Number of alignments\n  \
  \"scores|s=i\"      => \\$params{scores}, # Numb\
er of scores\n    \"sensitivity|S=s\" => \\$params\
{sensitivity}, # Search sensitivity\n    \"sort|t=\
s\"	      => \\$params{sort}, # Sort hits by...\n \
   \"stats|T=s\"       => \\$params{stats}, # Scor\
ing statistic to use\n    \"strand|d=s\"      => \\
\$params{strand}, # Strand to use in DNA vs. DNA s\
earch\n    \"topcombon|c=i\"   => \\$params{topcom\
bon}, # Consistent sets of HSPs\n    \"outfile=s\"\
       => \\$outfile, # Output file\n    \"outform\
at|o=s\"   => \\$outformat, # Output format\n    \\
"help|h\"	      => \\$help, # Usage info\n    \"as\
ync|a\"	      => \\$async, # Asynchronous mode\n  \
  \"polljob\"	      => \\$polljob, # Get results\n\
    \"status\"	      => \\$status, # Get job statu\
s\n    \"ids\"             => \\$ids, # Get ids fr\
om result\n    \"jobid|j=s\"       => \\$jobid, # \
JobId\n    \"email=s\"         => \\$params{email}\
, # E-mail address\n    \"trace\"           => \\$\
trace, # SOAP trace\n    \"sequence=s\"      => \\\
$sequence, # Query sequence\n    );\n\nmy $scriptN\
ame = basename($0, ());\nif($help || $numOpts == 0\
) {\n    &usage();\n    exit(0);\n}\n\nif($trace){\
\n    print STDERR \"Tracing active\\n\";\n    SOA\
P::Lite->import(+trace => 'debug');\n}\n\nmy $soap\
 = SOAP::Lite\n    ->service($WSDL)\n    ->proxy('\
http://localhost/',\n    #proxy => ['http' => 'htt\
p://your.proxy.server/'], # HTTP proxy\n    timeou\
t => 600, # HTTP connection timeout\n    )\n    ->\
on_fault(sub { # SOAP fault handler\n        my $s\
oap = shift;\n        my $res = shift;\n        # \
Throw an exception for all faults\n        if(ref(\
$res) eq '') {\n            die($res);\n        } \
else {\n            die($res->faultstring);\n     \
   }\n        return new SOAP::SOM;\n    }\n      \
         );\n\nif( !($polljob || $status || $ids) \
&&\n    !( defined($ARGV[0]) || defined($sequence)\
 )\n    ) {\n    print STDERR 'Error: bad option c\
ombination', \"\\n\";\n    &usage();\n    exit(1);\
\n}\nelsif($polljob && defined($jobid)) {\n    pri\
nt \"Getting results for job $jobid\\n\";\n    get\
Results($jobid);\n}\nelsif($status && defined($job\
id)) {\n    print STDERR \"Getting status for job \
$jobid\\n\";\n    my $result = $soap->checkStatus(\
$jobid);\n    print STDOUT \"$result\\n\";\n    if\
($result eq 'DONE') {\n	print STDERR \"To get resu\
lts: $scriptName --polljob --jobid $jobid\\n\";\n \
   }\n}  \nelsif($ids && defined($jobid)) {\n    p\
rint STDERR \"Getting ids from job $jobid\\n\";\n \
   getIds($jobid);\n}\nelse {\n    # Prepare input\
 data\n    my $content;\n    my (@contents) = ();\\
n    if(-f $ARGV[0] || $ARGV[0] eq '-') {	\n	$cont\
ent={type=>'sequence',content=>read_file($ARGV[0])\
};	\n    }\n    if($sequence) {	\n	if(-f $sequence\
 || $sequence eq '-') {	\n	    $content={type=>'se\
quence',content=>read_file($ARGV[0])};	\n	} else {\
\n	    $content={type=>'sequence',content=>$sequen\
ce};\n	}\n    }\n    push @contents, $content;\n\n\
    # Submit the job\n    my $paramsData = SOAP::D\
ata->name('params')->type(map=>\\%params);\n    my\
 $contentData = SOAP::Data->name('content')->value\
(\\@contents);\n    # For SOAP::Lite 0.60 and earl\
ier parameters are passed directly\n    if($SOAP::\
Lite::VERSION eq '0.60' || $SOAP::Lite::VERSION =~\
 /0\\.[1-5]/) {\n        $jobid = $soap->runWUBlas\
t($paramsData, $contentData);\n    }\n    # For SO\
AP::Lite 0.69 and later parameter handling is diff\
erent, so pass\n    # undef's for templated params\
, and then pass the formatted args.\n    else {\n \
       $jobid = $soap->runWUBlast(undef, undef,\n	\
			   $paramsData, $contentData);\n    }\n\n    # \
Asynchronous mode: output jobid and exit.\n    if \
(defined($async)) {\n	print STDOUT $jobid, \"\\n\"\
;\n        print STDERR \"To check status: $script\
Name --status --jobid $jobid\\n\";\n    }\n    # S\
ynchronous mode: try to get results\n    else {\n \
       print STDERR \"JobId: $jobid\\n\";\n       \
 sleep 1;\n        getResults($jobid);\n    }\n}\n\
\nsub getIds($) {\n    my $jobid = shift;\n    my \
$results = $soap->getIds($jobid);\n    for my $res\
ult (@$results){\n	print \"$result\\n\";\n    }\n}\
\n\nsub clientPoll($) {\n    my $jobid = shift;\n \
   my $result = 'PENDING';\n    # Check status and\
 wait if not finished\n    while($result eq 'RUNNI\
NG' || $result eq 'PENDING') {\n        $result = \
$soap->checkStatus($jobid);\n        print STDERR \
\"$result\\n\";\n        if($result eq 'RUNNING' |\
| $result eq 'PENDING') {\n            # Wait befo\
re polling again.\n            sleep $checkInterva\
l;\n        }\n    }\n}\n\nsub getResults($) {\n  \
  my $jobid = shift;\n    my $res;\n    # Check st\
atus, and wait if not finished\n    clientPoll($jo\
bid);\n    # Use JobId if output file name is not \
defined\n    unless(defined($outfile)) {\n        \
$outfile=$jobid;\n    }\n    # Get list of data ty\
pes\n    my $resultTypes = $soap->getResults($jobi\
d);\n    # Get the data and write it to a file\n  \
  if(defined($outformat)) { # Specified data type\\
n	if($outformat eq 'xml') {$outformat = 'toolxml';\
}\n	if($outformat eq 'txt') {$outformat = 'toolout\
put';}\n        my $selResultType;\n        foreac\
h my $resultType (@$resultTypes) {\n            if\
($resultType->{type} eq $outformat) {\n           \
     $selResultType = $resultType;\n            }\\
n        }\n        $res=$soap->poll($jobid, $selR\
esultType->{type});\n	if($outfile eq '-') {\n	    \
 write_file($outfile, $res);\n	} else {\n	    writ\
e_file($outfile.'.'.$selResultType->{ext}, $res);\\
n	}\n    } else { # Data types available\n        \
# Write a file for each output type\n        for m\
y $resultType (@$resultTypes){\n            #print\
 STDERR \"Getting $resultType->{type}\\n\";\n     \
       $res=$soap->poll($jobid, $resultType->{type\
});\n	    if($outfile eq '-') {\n		write_file($out\
file, $res);\n	    } else {\n		write_file($outfile\
.'.'.$resultType->{ext}, $res);\n	    }\n        }\
\n    }\n}\n\nsub read_file($) {\n    my $filename\
 = shift;\n    my ($content, $buffer);\n    if($fi\
lename eq '-') {\n	while(sysread(STDIN, $buffer, 1\
024)) {\n	    $content .= $buffer;\n	}\n    }\n   \
 else { # File\n	open(FILE, $filename) or die \"Er\
ror: unable to open input file\";\n	while(sysread(\
FILE, $buffer, 1024)) {\n	    $content .= $buffer;\
\n	}\n	close(FILE);\n    }\n    return $content;\n\
}\n\nsub write_file($$) {\n    my ($filename, $dat\
a) = @_;\n    print STDERR 'Creating result file: \
' . $filename . \"\\n\";\n    if($filename eq '-')\
 {\n	print STDOUT $data;\n    }\n    else {\n	open\
(FILE, \">$filename\") or die \"Error: unable to o\
pen output file\";\n	syswrite(FILE, $data);\n	clos\
e(FILE);\n    }\n}\n\nsub usage {\n    print STDER\
R <<EOF\nWU-BLAST\n========\n\nRapid sequence data\
base search programs utilizing the BLAST algorithm\
.\n   \n[Required]\n\n      --email       : str  :\
 user email address \n  -p, --program	    : str  :\
 BLAST program to use: blastn, blastp, blastx, \n \
                            tblastn or tblastx\n  \
-D, --database    : str  : database to search\n  s\
eqFile           : file : query sequence data file\
 (\"-\" for STDIN)\n\n[Optional]\n\n  -m, --matrix\
	    : str  : scoring matrix\n  -E, --exp	    : re\
al : 0<E<= 1000. Statistical significance threshol\
d\n                             for reporting data\
base sequence matches.\n  -e, --echofilter  :     \
 : display the filtered query sequence in the outp\
ut\n  -f, --filter	    : str  : activates filterin\
g of the query sequence\n  -b, --alignments  : int\
  : number of alignments to be reported\n  -s, --s\
cores	    : int  : number of scores to be reported\
\n  -S, --sensitivity : str  :\n  -t, --sort	    :\
 str  :\n  -T, --stats       : str  :\n  -d, --str\
and      : str  : DNA strand to search with in DNA\
 vs. DNA searches \n  -c, --topcombon   :      :\n\
\n[General]	\n\n  -h, --help       :      : prints\
 this help text\n  -a, --async      :      : force\
s to make an asynchronous query\n      --status   \
  :      : poll for the status of a job\n      --p\
olljob    :      : poll for the results of a job\n\
  -j, --jobid      : str  : jobid that was returne\
d when an asynchronous job \n                     \
       was submitted.\n  -O, --outfile    : str  :\
 name of the file results should be written to \n \
                           (default is based on th\
e jobid; \"-\" for STDOUT)\n  -o, --outformat  : s\
tr  : txt or xml output (no file is written)\n    \
  --trace	   :      : show SOAP messages being int\
erchanged \n\nSynchronous job:\n\n  The results/er\
rors are returned as soon as the job is finished.\\
n  Usage: $scriptName --email <your\\@email> [opti\
ons...] seqFile\n  Returns: saves the results to d\
isk\n\nAsynchronous job:\n\n  Use this if you want\
 to retrieve the results at a later time. The resu\
lts \n  are stored for up to 24 hours. \n  The asy\
nchronous submission mode is recommended when user\
s are submitting \n  batch jobs or large database \
searches	\n  Usage: $scriptName --async --email <y\
our\\@email> [options...] seqFile\n  Returns : job\
id\n\n  Use the jobid to query for the status of t\
he job. \n  Usage: $scriptName --status --jobid <j\
obId>\n  Returns : string indicating the status of\
 the job:\n    DONE - job has finished\n    RUNNIN\
G - job is running\n    NOT_FOUND - job cannot be \
found\n    ERROR - the jobs has encountered an err\
or\n\n  When done, use the jobid to retrieve the s\
tatus of the job. \n  Usage: $scriptName --polljob\
 --jobid <jobId> [--outfile string]\n  Returns: sa\
ves the results to disk\n\n[Help]\n\nFor more deta\
iled help information refer to \nhttp://www.ebi.ac\
.uk/blast2/WU-Blast2_Help_frame.html\n \nEOF\n;\n}\
\n","\nmy $WSDL = 'http://www.ebi.ac.uk/Tools/webs\
ervices/wsdl/WSBlastpgp.wsdl';\n\nuse SOAP::Lite;\\
nuse Getopt::Long qw(:config no_ignore_case bundli\
ng);\nuse File::Basename;\n\nmy $checkInterval = 1\
5;\n\nmy %params=(\n	    'async' => '1', # Use asy\
nc mode and simulate sync mode in client\n	    );\\
nGetOptions(\n    \"mode=s\"           => \\$param\
s{mode}, # Search mode: PSI-Blast or PHI-Blast\n  \
  \"database|d=s\"     => \\$params{database}, # D\
atabase to search\n    \"matrix|M=s\"       => \\$\
params{matrix},# Scoring maxtrix\n    \"exp|e=f\" \
         => \\$params{exp}, # E-value\n    \"expmu\
lti|h=f\"     => \\$params{expmulti}, # E-value\n \
   \"filter|F=s\"       => \\$params{filter}, # Lo\
w complexity filter\n    \"dropoff|X=i\"      => \\
\$params{dropoff}, # Dropoff score\n    \"finaldro\
poff|Z=i\" => \\$params{finaldropoff}, # Final dro\
poff score\n    \"scores|v=i\"       => \\$params{\
scores}, # Max number of scores\n    \"align=i\"  \
        => \\$params{align}, # Alignment view\n   \
 \"startregion|S=i\"  => \\$params{startregion}, #\
 Start of region in query\n    \"endregion|H=i\"  \
  => \\$params{endregion}, # End of region in quer\
y\n    \"maxpasses|j=i\"    => \\$params{maxpasses\
}, # Number of PSI iterations\n    \"opengap|G=i\"\
      => \\$params{opengap}, # Gap open penalty\n \
   \"extendgap|E=i\"    => \\$params{extendgap}, #\
 Gap extension penalty\n    \"pattern=s\"        =\
> \\$params{pattern}, # PHI-BLAST pattern\n    \"u\
sagemode|p=s\"    => \\$params{usagemode}, # PHI-B\
LAST program\n    \"appxml=s\"         => \\$param\
s{appxml}, # Application XML\n    \"sequence=s\"  \
     => \\$sequence, # Query sequence\n    \"help\\
"	       => \\$help, # Usage info\n    \"polljob\"\
	       => \\$polljob, # Get results\n    \"status\
\"	       => \\$status, # Get status\n    \"ids\" \
     	       => \\$ids, # Get ids from result\n   \
 \"jobid=s\"          => \\$jobid, # JobId\n    \"\
outfile=s\"        => \\$outfile, # Output filenam\
e\n    \"outformat|o=s\"    => \\$outformat, # Out\
put file format\n    \"async|a\"	       => \\$asyn\
c, # Async submission\n    \"email=s\"          =>\
 \\$params{email}, # User e-mail address\n    \"tr\
ace\"            => \\$trace, # Show SOAP messages\
\n    );\n\nmy $scriptName = basename($0, ());\nif\
($help) {\n    &usage();\n    exit(0);\n}\n\nif ($\
trace){\n    print \"Tracing active\\n\";\n    SOA\
P::Lite->import(+trace => 'debug');\n}\n\nmy $soap\
 = SOAP::Lite\n    ->service($WSDL)\n    ->on_faul\
t(sub {\n        my $soap = shift;\n        my $re\
s = shift;\n        # Throw an exception for all f\
aults\n        if(ref($res) eq '') {\n            \
die($res);\n        } else {\n            die($res\
->faultstring);\n        }\n        return new SOA\
P::SOM;\n    }\n               );\n\nif( !($polljo\
b || $status || $ids) &&\n    !( (defined($ARGV[0]\
) && -f $ARGV[0]) || defined($sequence) )\n    ) {\
\n    print STDERR 'Error: bad option combination'\
, \"\\n\";\n    &usage();\n    exit(1);\n}\nelsif(\
$polljob && defined($jobid)) {\n    print \"Gettin\
g results for job $jobid\\n\";\n    getResults($jo\
bid);\n}\nelsif($status && defined($jobid)) {\n   \
 print STDERR \"Getting status for job $jobid\\n\"\
;\n    my $result = $soap->checkStatus($jobid);\n \
   print STDOUT $result, \"\\n\";\n    if($result \
eq 'DONE') {\n	print STDERR \"To get results: $scr\
iptName --polljob --jobid $jobid\\n\";\n    }\n}  \
\nelsif($ids && defined($jobid)) {\n    print STDE\
RR \"Getting ids from job $jobid\\n\";\n    getIds\
($jobid);\n}\nelse {\n    if(-f $ARGV[0]) {	\n	$co\
ntent={type=>'sequence', content=>read_file($ARGV[\
0])};	\n    }\n    if($sequence) {	\n	if(-f $seque\
nce) {\n	    $content={type=>'sequence', content=>\
read_file($sequence)};	\n	} else {\n	    $content=\
{type=>'sequence', content=>$sequence};\n	}\n    }\
\n    push @content, $content;\n\n    my $jobid;\n\
    my $paramsData = SOAP::Data->name('params')->t\
ype(map=>\\%params);\n    my $contentData = SOAP::\
Data->name('content')->value(\\@content);\n    # F\
or SOAP::Lite 0.60 and earlier parameters are pass\
ed directly\n    if($SOAP::Lite::VERSION eq '0.60'\
 || $SOAP::Lite::VERSION =~ /0\\.[1-5]/) {\n      \
  $jobid = $soap->runBlastpgp($paramsData, $conten\
tData);\n    }\n    # For SOAP::Lite 0.69 and late\
r parameter handling is different, so pass\n    # \
undef's for templated params, and then pass the fo\
rmatted args.\n    else {\n        $jobid = $soap-\
>runBlastpgp(undef, undef,\n				    $paramsData, $\
contentData);\n    }\n\n    if (defined($async)) {\
\n	print STDOUT $jobid, \"\\n\";\n        print ST\
DERR \"To check status: $scriptName --status --job\
id $jobid\\n\";\n    } else { # Synchronous mode\n\
        print STDERR \"JobId: $jobid\\n\";\n      \
  sleep 1;\n        getResults($jobid);\n    }\n}\\
n\nsub getIds($) {\n    $jobid = shift;\n    my $r\
esults = $soap->getIds($jobid);\n    for $result (\
@$results){\n	print \"$result\\n\";\n    }\n}\n\ns\
ub clientPoll($) {\n    my $jobid = shift;\n    my\
 $result = 'PENDING';\n    # Check status and wait\
 if not finished\n    #print STDERR \"Checking sta\
tus: $jobid\\n\";\n    while($result eq 'RUNNING' \
|| $result eq 'PENDING') {\n        $result = $soa\
p->checkStatus($jobid);\n        print STDERR \"$r\
esult\\n\";\n        if($result eq 'RUNNING' || $r\
esult eq 'PENDING') {\n            # Wait before p\
olling again.\n            sleep $checkInterval;\n\
        }\n    }\n}\n\nsub getResults($) {\n    $j\
obid = shift;\n    # Check status, and wait if not\
 finished\n    clientPoll($jobid);\n    # Use JobI\
d if output file name is not defined\n    unless(d\
efined($outfile)) {\n        $outfile=$jobid;\n   \
 }\n    # Get list of data types\n    my $resultTy\
pes = $soap->getResults($jobid);\n    # Get the da\
ta and write it to a file\n    if(defined($outform\
at)) { # Specified data type\n        my $selResul\
tType;\n        foreach my $resultType (@$resultTy\
pes) {\n            if($resultType->{type} eq $out\
format) {\n                $selResultType = $resul\
tType;\n            }\n        }\n        $res=$so\
ap->poll($jobid, $selResultType->{type});\n       \
 write_file($outfile.'.'.$selResultType->{ext}, $r\
es);\n    } else { # Data types available\n       \
 # Write a file for each output type\n        for \
my $resultType (@$resultTypes){\n            #prin\
t \"Getting $resultType->{type}\\n\";\n           \
 $res=$soap->poll($jobid, $resultType->{type});\n \
           write_file($outfile.'.'.$resultType->{e\
xt}, $res);\n        }\n    }\n}\n\nsub read_file(\
$) {\n    my $filename = shift;\n    open(FILE, $f\
ilename);\n    my $content;\n    my $buffer;\n    \
while(sysread(FILE, $buffer, 1024)) {\n	$content.=\
 $buffer;\n    }\n    close(FILE);  \n    return $\
content;\n}\n\nsub write_file($$) {\n    my ($tmp,\
$entity) = @_;\n    print STDERR \"Creating result\
 file: \".$tmp.\"\\n\";\n    unless(open (FILE, \"\
>$tmp\")) {\n	return 0;\n    }\n    syswrite(FILE,\
 $entity);\n    close (FILE);\n    return 1;\n}\n\\
nsub usage {\n    print STDERR <<EOF\nBlastpgp\n==\
======\n   \nThe blastpgp program implements the P\
SI-BLAST and PHI-BLAST variations\nof NCBI BLAST.\\
n\nFor more detailed help information refer to\nht\
tp://www.ebi.ac.uk/blastpgp/blastpsi_help_frame.ht\
ml\n \nBlastpgp specific options:\n\n[Required]\n\\
n      --mode            : str  : search mode to u\
se: PSI-Blast or PHI-Blast\n  -d, --database      \
  : str  : protein database to search\n  seqFile  \
             : file : query sequence\n\n[Optional]\
\n\n  -M, --matrix          : str  : scoring matri\
x\n  -e, --exp             : real : Expectation va\
lue\n  -h, --expmulti        : real : threshold (m\
ultipass model)\n  -F, --filter          : str  : \
filter query sequence with SEG [T,F]\n  -m, --alig\
n           : int  : alignment view option:\n     \
                            0 - pairwise, 1 - M/S \
identities,\n                                 2 - \
M/S non-identities, 3 - Flat identities,\n        \
                         4 - Flat non-identities\n\
  -G, --opengap         : int  : cost to open a ga\
p\n  -E, --extendgap       : int  : cost to extend\
 a gap\n  -g, --gapalign        : str  : Gapped [T\
,F]\n  -v, --scores          : int  : number of sc\
ores to be reported\n  -j, --maxpasses       : int\
  : number of iterations\n  -X, --dropoff         \
: int  : Dropoff score\n  -Z, --finaldropoff    : \
int  : Dropoff for final alignment\n  -S, --startr\
egion     : int  : Start of required region in que\
ry\n  -H, --endregion       : int  : End of requir\
ed region in query\n  -k, --pattern         : str \
 : Hit File (PHI-BLAST only)\n  -p, --usagemode   \
    : str  : Program option (PHI-BLAST only):\n   \
                              blastpgp, patseedp, \
seedp\n\n[General]\n\n      --help            :   \
   : prints this help text\n  -a, --async         \
  :      : forces to make an asynchronous query\n \
     --status          :      : poll for the statu\
s of a job\n      --polljob         :      : poll \
for the results of a job\n      --jobid           \
: str  : jobid of an asynchronous job\n      --ids\
             :      : get hit identifiers for resu\
lt \n  -O, --outfile         : str  : name of the \
file results should be written to\n               \
                  (default is based on the jobid)\\
n  -o, --outformat       : str  : txt or xml outpu\
t (no file is written)\n      --trace           : \
     : show SOAP messages being interchanged\n\nSy\
nchronous job:\n\n  The results/errors are returne\
d as soon as the job is finished.\n  Usage: blastp\
gp.pl --email <your@email> [options...] seqfile\n \
 Returns: saves the results to disk\n\nAsynchronou\
s job:\n\n  Use this if you want to retrieve the r\
esults at a later time. The results\n  are stored \
for up to 24 hours.\n  The asynchronous submission\
 mode is recommended when users are submitting\n  \
batch jobs or large database searches\n  Usage: bl\
astpgp.pl --email <your@email> --async [options...\
] seqFile\n  Returns: jobid\n\n  Use the jobid to \
query for the status of the job.\n  Usage: blastpg\
p.pl --status --jobid <jobId>\n  Returns: string i\
ndicating the status of the job\n    DONE - job ha\
s finished\n    RUNNING - job is running\n    NOT_\
FOUND - job cannot be found\n    ERROR - the jobs \
has encountered an error\n\n  When done, use the j\
obid to retrieve the results of the job.\n  Usage:\
 blastpgp.pl --polljob --jobid <jobId> [--outfile \
<fileName>]\n  Returns: saves the results to disk\\
nEOF\n;\n}\n","\n=head1 NAME\n\nncbiblast.pl\n\n=h\
ead1 DESCRIPTION\n\nNCBI Blast (REST) web service \
Perl client using L<LWP>.\n\nTested with:\n\n=over\
\n\n=item *\nL<LWP> 6.35, L<XML::Simple> 2.25 and \
Perl 5.22.0 (MacOS 10.13.6)\n\n=back\n\nFor furthe\
r information see:\n\n=over\n\n=item *\nL<https://\
www.ebi.ac.uk/Tools/webservices/>\n\n=back\n\n=hea\
d1 LICENSE\n\nCopyright 2012-2018 EMBL - European \
Bioinformatics Institute\n\nLicensed under the Apa\
che License, Version 2.0 (the \"License\");\nyou m\
ay not use this file except in compliance with the\
 License.\nYou may obtain a copy of the License at\
\n\n    http://www.apache.org/licenses/LICENSE-2.0\
\n\nUnless required by applicable law or agreed to\
 in writing, software\ndistributed under the Licen\
se is distributed on an \"AS IS\" BASIS,\nWITHOUT \
WARRANTIES OR CONDITIONS OF ANY KIND, either expre\
ss or implied.\nSee the License for the specific l\
anguage governing permissions and\nlimitations und\
er the License.\n\nPerl Client Automatically gener\
ated with:\nhttps://github.com/ebi-wp/webservice-c\
lients-generator\n\n=cut\n\nuse strict;\nuse warni\
ngs;\n\nuse English;\nuse LWP;\nuse XML::Simple;\n\
use Getopt::Long qw(:config no_ignore_case bundlin\
g);\nuse File::Basename;\nuse Data::Dumper;\nuse T\
ime::HiRes qw(usleep);\n\nmy $baseUrl = 'https://w\
ww.ebi.ac.uk/Tools/services/rest/ncbiblast';\nmy $\
version = '2020-02-28 14:43';\n\nmy $checkInterval\
 = 3;\n\nmy $maxErrorStatusCount = 3;\n\nmy $outpu\
tLevel = 1;\n\nmy $numOpts = scalar(@ARGV);\nmy %p\
arams = (\n    'debugLevel' => 0,\n    'maxJobs'  \
  => 1\n);\nmy @database;\nGetOptions(\n    # Tool\
 specific options\n    'program=s'       => \\$par\
ams{'program'},        # The BLAST program to be u\
sed for the Sequence Similarity Search.\n    'task\
=s'          => \\$params{'task'},           # Tas\
k option (only selectable for blastn)\n    'matrix\
=s'        => \\$params{'matrix'},         # (Prot\
ein searches) The substitution matrix used for sco\
ring alignments when searching the database.\n    \
'alignments=i'    => \\$params{'alignments'},     \
# Maximum number of match alignments reported in t\
he result output.\n    'scores=i'        => \\$par\
ams{'scores'},         # Maximum number of match s\
core summaries reported in the result output.\n   \
 'exp=s'           => \\$params{'exp'},           \
 # Limits the number of scores and alignments repo\
rted based on the expectation value. This is the m\
aximum number of times the match is expected to oc\
cur by chance.\n    'dropoff=i'       => \\$params\
{'dropoff'},        # The amount a score can drop \
before gapped extension of word hits is halted\n  \
  'match_scores=s'  => \\$params{'match_scores'}, \
  # (Nucleotide searches) The match score is the b\
onus to the alignment score when matching the same\
 base. The mismatch is the penalty when failing to\
 match.\n    'gapopen=i'       => \\$params{'gapop\
en'},        # Penalty taken away from the score w\
hen a gap is created in sequence. Increasing the g\
ap openning penalty will decrease the number of ga\
ps in the final alignment.\n    'gapext=i'        \
=> \\$params{'gapext'},         # Penalty taken aw\
ay from the score for each base or residue in the \
gap. Increasing the gap extension penalty favors s\
hort gaps in the final alignment, conversly decrea\
sing the gap extension penalty favors long gaps in\
 the final alignment.\n    'filter=s'        => \\\
$params{'filter'},         # Filter regions of low\
 sequence complexity. This can avoid issues with l\
ow complexity sequences where matches are found du\
e to composition rather than meaningful sequence s\
imilarity. However in some cases filtering also ma\
sks regions of interest and so should be used with\
 caution.\n    'seqrange=s'      => \\$params{'seq\
range'},       # Specify a range or section of the\
 input sequence to use in the search. Example: Spe\
cifying '34-89' in an input sequence of total leng\
th 100, will tell BLAST to only use residues 34 to\
 89, inclusive.\n    'gapalign'        => \\$param\
s{'gapalign'},       # This is a true/false settin\
g that tells the program the perform optimised ali\
gnments within regions involving gaps. If set to t\
rue, the program will perform an alignment using g\
aps. Otherwise, if it is set to false, it will rep\
ort only individual HSP where two sequence match e\
ach other, and thus will not produce alignments wi\
th gaps.\n    'wordsize=i'      => \\$params{'word\
size'},       # Word size for wordfinder algorithm\
\n    'taxids=s'        => \\$params{'taxids'},   \
      # Specify one or more TaxIDs so that the BLA\
ST search becomes taxonomically aware.\n    'comps\
tats=s'     => \\$params{'compstats'},      # Use \
composition-based statistics.\n    'align=i'      \
   => \\$params{'align'},          # Formating for\
 the alignments\n    'transltable=i'   => \\$param\
s{'transltable'},    # Query Genetic code to use i\
n translation\n    'stype=s'         => \\$params{\
'stype'},          # Indicates if the sequence is \
protein or DNA/RNA.\n    'sequence=s'      => \\$p\
arams{'sequence'},       # The query sequence can \
be entered directly into this form. The sequence c\
an be in GCG, FASTA, EMBL (Nucleotide only), GenBa\
nk, PIR, NBRF, PHYLIP or UniProtKB/Swiss-Prot (Pro\
tein only) format. A partially formatted sequence \
is not accepted. Adding a return to the end of the\
 sequence may help certain applications understand\
 the input. Note that directly using data from wor\
d processors may yield unpredictable results as hi\
dden/control characters may be present.\n    'data\
base=s'      => \\@database,                 # Dat\
abase\n    # Generic options\n    'email=s'       \
  => \\$params{'email'},          # User e-mail ad\
dress\n    'title=s'         => \\$params{'title'}\
,          # Job title\n    'outfile=s'       => \\
\$params{'outfile'},        # Output file name\n  \
  'outformat=s'     => \\$params{'outformat'},    \
  # Output file type\n    'jobid=s'         => \\$\
params{'jobid'},          # JobId\n    'help|h'   \
       => \\$params{'help'},           # Usage hel\
p\n    'asyncjob'        => \\$params{'asyncjob'},\
       # Asynchronous submission\n    'polljob'   \
      => \\$params{'polljob'},        # Get result\
s\n    'pollFreq=f'      => \\$params{'pollFreq'},\
       # Poll Frequency\n    'resultTypes'     => \
\\$params{'resultTypes'},    # Get result types\n \
   'status'          => \\$params{'status'},      \
   # Get status\n    'params'          => \\$param\
s{'params'},         # List input parameters\n    \
'paramDetail=s'   => \\$params{'paramDetail'},    \
# Get details for parameter\n    'multifasta'     \
 => \\$params{'multifasta'},     # Multiple fasta \
input\n    'useSeqId'        => \\$params{'useSeqI\
d'},       # Seq Id file name\n    'maxJobs=i'    \
   => \\$params{'maxJobs'},        # Max. parallel\
 jobs\n\n    'verbose'         => \\$params{'verbo\
se'},        # Increase output level\n    'version\
'         => \\$params{'version'},        # Prints\
 out the version of the Client and exit.\n    'qui\
et'           => \\$params{'quiet'},          # De\
crease output level\n    'debugLevel=i'    => \\$p\
arams{'debugLevel'},     # Debugging level\n    'b\
aseUrl=s'       => \\$baseUrl,                  # \
Base URL for service.\n);\nif ($params{'verbose'})\
 {$outputLevel++}\nif ($params{'quiet'}) {$outputL\
evel--}\nif ($params{'pollFreq'}) {$checkInterval \
= $params{'pollFreq'} * 1000 * 1000}\nif ($params{\
'baseUrl'}) {$baseUrl = $params{'baseUrl'}}\n$para\
ms{\"database\"} = [@database];\n\n&print_debug_me\
ssage('MAIN', 'LWP::VERSION: ' . $LWP::VERSION,\n \
   1);\n\n&print_debug_message('MAIN', \"params:\\\
n\" . Dumper(\\%params), 11);\n\nmy $ua;\n\nmy $sc\
riptName = basename($0, ());\n\nif ($params{'help'\
} || $numOpts == 0) {\n    &usage();\n    exit(0);\
\n}\n\n&print_debug_message('MAIN', 'baseUrl: ' . \
$baseUrl, 1);\nif (\n    !(\n        $params{'poll\
job'}\n            || $params{'resultTypes'}\n    \
        || $params{'status'}\n            || $para\
ms{'params'}\n            || $params{'paramDetail'\
}\n            || $params{'version'}\n    )\n     \
   && !(defined($ARGV[0]) || defined($params{'sequ\
ence'}))\n) {\n\n    # Bad argument combination, s\
o print error message and usage\n    print STDERR \
'Error: bad option combination', \"\\n\";\n    &us\
age();\n    exit(1);\n}\nelsif ($params{'params'})\
 {\n    &print_tool_params();\n}\n\nelsif ($params\
{'paramDetail'}) {\n    &print_param_details($para\
ms{'paramDetail'});\n}\n\nelsif ($params{'version'\
}) {\n  print STDOUT 'Revision: ' . $version, \"\\\
n\";\n  exit(1);\n}\n\nelsif ($params{'status'} &&\
 defined($params{'jobid'})) {\n    &print_job_stat\
us($params{'jobid'});\n}\n\nelsif ($params{'result\
Types'} && defined($params{'jobid'})) {\n    &prin\
t_result_types($params{'jobid'});\n}\n\nelsif ($pa\
rams{'polljob'} && defined($params{'jobid'})) {\n \
   &get_results($params{'jobid'});\n}\n\nelse {\n \
   # Multiple input sequence mode, assume fasta fo\
rmat.\n    if (defined($params{'multifasta'}) && $\
params{'multifasta'}) {\n        &multi_submit_job\
();\n    }\n\n    # Entry identifier list file.\n \
   elsif ((defined($params{'sequence'}) && $params\
{'sequence'} =~ m/^\\@/)\n        || (defined($ARG\
V[0]) && $ARGV[0] =~ m/^\\@/)) {\n        my $list\
_filename = $params{'sequence'} || $ARGV[0];\n    \
    $list_filename =~ s/^\\@//;\n        &list_fil\
e_submit_job($list_filename);\n    }\n    # Defaul\
t: single sequence/identifier.\n    else {\n      \
  # Warn for invalid batch only option use.\n     \
   if (defined($params{'useSeqId'}) && $params{'us\
eSeqId'}) {\n            print STDERR \"Warning: -\
-useSeqId option ignored.\\n\";\n            delet\
e $params{'useSeqId'};\n        }\n        if (def\
ined($params{'maxJobs'}) && $params{'maxJobs'} > 1\
) {\n            print STDERR \"Warning: --maxJobs\
 option ignored.\\n\";\n            $params{'maxJo\
bs'} = 1;\n        }\n        # Load the sequence \
data and submit.\n        &submit_job(&load_data()\
);\n    }\n}\n\n\n\n=head1 FUNCTIONS\n\n=cut\n\n\n\
=head2 rest_user_agent()\n\nGet a LWP UserAgent to\
 use to perform REST requests.\n\n  my $ua = &rest\
_user_agent();\n\n=cut\n\nsub rest_user_agent() {\\
n    print_debug_message('rest_user_agent', 'Begin\
', 21);\n    # Create an LWP UserAgent for making \
HTTP calls.\n    my $ua = LWP::UserAgent->new();\n\
    # Set 'User-Agent' HTTP header to identifiy th\
e client.\n    my $revisionNumber = 0;\n    $revis\
ionNumber = \"Revision: \" . $version;\n    $ua->a\
gent(\"EBI-Sample-Client/$revisionNumber ($scriptN\
ame; $OSNAME) \" . $ua->agent());\n    # Configure\
 HTTP proxy support from environment.\n    $ua->en\
v_proxy;\n    print_debug_message('rest_user_agent\
', 'End', 21);\n    return $ua;\n}\n\n=head2 rest_\
error()\n\nCheck a REST response for an error cond\
ition. An error is mapped to a die.\n\n  &rest_err\
or($response, $content_data);\n\n=cut\n\nsub rest_\
error() {\n    print_debug_message('rest_error', '\
Begin', 21);\n    my $response = shift;\n    my $c\
ontentdata;\n    if (scalar(@_) > 0) {\n        $c\
ontentdata = shift;\n    }\n    if (!defined($cont\
entdata) || $contentdata eq '') {\n        $conten\
tdata = $response->content();\n    }\n    # Check \
for HTTP error codes\n    if ($response->is_error)\
 {\n        my $error_message = '';\n        # HTM\
L response.\n        if ($contentdata =~ m/<h1>([^\
<]+)<\\/h1>/) {\n            $error_message = $1;\\
n        }\n        #  XML response.\n        elsi\
f ($contentdata =~ m/<description>([^<]+)<\\/descr\
iption>/) {\n            $error_message = $1;\n   \
     }\n        # die 'http status: ' . $response-\
>code . ' ' . $response->message . '  ' . $error_m\
essage;\n    }\n    print_debug_message('rest_erro\
r', 'End', 21);\n}\n\n=head2 rest_request()\n\nPer\
form a REST request (HTTP GET).\n\n  my $response_\
str = &rest_request($url);\n\n=cut\n\nsub rest_req\
uest {\n    print_debug_message('rest_request', 'B\
egin', 11);\n    my $requestUrl = shift;\n    prin\
t_debug_message('rest_request', 'URL: ' . $request\
Url, 11);\n\n    # Get an LWP UserAgent.\n    $ua \
= &rest_user_agent() unless defined($ua);\n    # A\
vailable HTTP compression methods.\n    my $can_ac\
cept;\n    eval {\n        $can_accept = HTTP::Mes\
sage::decodable();\n    };\n    $can_accept = '' u\
nless defined($can_accept);\n    # Perform the req\
uest\n    my $response = $ua->get($requestUrl,\n  \
      'Accept-Encoding' => $can_accept, # HTTP com\
pression.\n    );\n    print_debug_message('rest_r\
equest', 'HTTP status: ' . $response->code,\n     \
   11);\n    print_debug_message('rest_request',\n\
        'response length: ' . length($response->co\
ntent()), 11);\n    print_debug_message('rest_requ\
est',\n        'request:' . \"\\n\" . $response->r\
equest()->as_string(), 32);\n    print_debug_messa\
ge('rest_request',\n        'response: ' . \"\\n\"\
 . $response->as_string(), 32);\n    # Unpack poss\
ibly compressed response.\n    my $retVal;\n    if\
 (defined($can_accept) && $can_accept ne '') {\n  \
      $retVal = $response->decoded_content();\n   \
 }\n    # If unable to decode use orginal content.\
\n    $retVal = $response->content() unless define\
d($retVal);\n    # Check for an error.\n    &rest_\
error($response, $retVal);\n    print_debug_messag\
e('rest_request', 'retVal: ' . $retVal, 12);\n    \
print_debug_message('rest_request', 'End', 11);\n\\
n    # Return the response data\n    return $retVa\
l;\n}\n\n=head2 rest_get_parameters()\n\nGet list \
of tool parameter names.\n\n  my (@param_list) = &\
rest_get_parameters();\n\n=cut\n\nsub rest_get_par\
ameters {\n    print_debug_message('rest_get_param\
eters', 'Begin', 1);\n    my $url = $baseUrl . '/p\
arameters/';\n    my $param_list_xml_str = rest_re\
quest($url);\n    my $param_list_xml = XMLin($para\
m_list_xml_str);\n    my (@param_list) = @{$param_\
list_xml->{'id'}};\n    print_debug_message('rest_\
get_parameters', 'End', 1);\n    return(@param_lis\
t);\n}\n\n=head2 rest_get_parameter_details()\n\nG\
et details of a tool parameter.\n\n  my $paramDeta\
il = &rest_get_parameter_details($param_name);\n\n\
=cut\n\nsub rest_get_parameter_details {\n    prin\
t_debug_message('rest_get_parameter_details', 'Beg\
in', 1);\n    my $parameterId = shift;\n    print_\
debug_message('rest_get_parameter_details',\n     \
   'parameterId: ' . $parameterId, 1);\n    my $ur\
l = $baseUrl . '/parameterdetails/' . $parameterId\
;\n    my $param_detail_xml_str = rest_request($ur\
l);\n    my $param_detail_xml = XMLin($param_detai\
l_xml_str);\n    print_debug_message('rest_get_par\
ameter_details', 'End', 1);\n    return($param_det\
ail_xml);\n}\n\n=head2 rest_run()\n\nSubmit a job.\
\n\n  my $job_id = &rest_run($email, $title, \\%pa\
rams );\n\n=cut\n\nsub rest_run {\n    print_debug\
_message('rest_run', 'Begin', 1);\n    my $email =\
 shift;\n    my $title = shift;\n    my $params = \
shift;\n    $email = '' if (!$email);\n    print_d\
ebug_message('rest_run', 'email: ' . $email, 1);\n\
    if (defined($title)) {\n        print_debug_me\
ssage('rest_run', 'title: ' . $title, 1);\n    }\n\
    print_debug_message('rest_run', 'params: ' . D\
umper($params), 1);\n\n    # Get an LWP UserAgent.\
\n    $ua = &rest_user_agent() unless defined($ua)\
;\n\n    # Clean up parameters\n    my (%tmp_param\
s) = %{$params};\n    $tmp_params{'email'} = $emai\
l;\n    $tmp_params{'title'} = $title;\n    foreac\
h my $param_name (keys(%tmp_params)) {\n        if\
 (!defined($tmp_params{$param_name})) {\n         \
   delete $tmp_params{$param_name};\n        }\n  \
  }\n\n    # Submit the job as a POST\n    my $url\
 = $baseUrl . '/run';\n    my $response = $ua->pos\
t($url, \\%tmp_params);\n    print_debug_message('\
rest_run', 'HTTP status: ' . $response->code, 11);\
\n    print_debug_message('rest_run',\n        're\
quest:' . \"\\n\" . $response->request()->as_strin\
g(), 11);\n    print_debug_message('rest_run',\n  \
      'response: ' . length($response->as_string()\
) . \"\\n\" . $response->as_string(), 11);\n\n    \
# Check for an error.\n    &rest_error($response);\
\n\n    # The job id is returned\n    my $job_id =\
 $response->content();\n    print_debug_message('r\
est_run', 'End', 1);\n    return $job_id;\n}\n\n=h\
ead2 rest_get_status()\n\nCheck the status of a jo\
b.\n\n  my $status = &rest_get_status($job_id);\n\\
n=cut\n\nsub rest_get_status {\n    print_debug_me\
ssage('rest_get_status', 'Begin', 1);\n    my $job\
_id = shift;\n    print_debug_message('rest_get_st\
atus', 'jobid: ' . $job_id, 2);\n    my $status_st\
r = 'UNKNOWN';\n    my $url = $baseUrl . '/status/\
' . $job_id;\n    $status_str = &rest_request($url\
);\n    print_debug_message('rest_get_status', 'st\
atus_str: ' . $status_str, 2);\n    print_debug_me\
ssage('rest_get_status', 'End', 1);\n    return $s\
tatus_str;\n}\n\n=head2 rest_get_result_types()\n\\
nGet list of result types for finished job.\n\n  m\
y (@result_types) = &rest_get_result_types($job_id\
);\n\n=cut\n\nsub rest_get_result_types {\n    pri\
nt_debug_message('rest_get_result_types', 'Begin',\
 1);\n    my $job_id = shift;\n    print_debug_mes\
sage('rest_get_result_types', 'jobid: ' . $job_id,\
 2);\n    my (@resultTypes);\n    my $url = $baseU\
rl . '/resulttypes/' . $job_id;\n    my $result_ty\
pe_list_xml_str = &rest_request($url);\n    my $re\
sult_type_list_xml = XMLin($result_type_list_xml_s\
tr);\n    (@resultTypes) = @{$result_type_list_xml\
->{'type'}};\n    print_debug_message('rest_get_re\
sult_types',\n        scalar(@resultTypes) . ' res\
ult types', 2);\n    print_debug_message('rest_get\
_result_types', 'End', 1);\n    return(@resultType\
s);\n}\n\n=head2 rest_get_result()\n\nGet result d\
ata of a specified type for a finished job.\n\n  m\
y $result = rest_get_result($job_id, $result_type)\
;\n\n=cut\n\nsub rest_get_result {\n    print_debu\
g_message('rest_get_result', 'Begin', 1);\n    my \
$job_id = shift;\n    my $type = shift;\n    print\
_debug_message('rest_get_result', 'jobid: ' . $job\
_id, 1);\n    print_debug_message('rest_get_result\
', 'type: ' . $type, 1);\n    my $url = $baseUrl .\
 '/result/' . $job_id . '/' . $type;\n    my $resu\
lt = &rest_request($url);\n    print_debug_message\
('rest_get_result', length($result) . ' characters\
',\n        1);\n    print_debug_message('rest_get\
_result', 'End', 1);\n    return $result;\n}\n\n\n\
=head2 print_debug_message()\n\nPrint debug messag\
e at specified debug level.\n\n  &print_debug_mess\
age($method_name, $message, $level);\n\n=cut\n\nsu\
b print_debug_message {\n    my $function_name = s\
hift;\n    my $message = shift;\n    my $level = s\
hift;\n    if ($level <= $params{'debugLevel'}) {\\
n        print STDERR '[', $function_name, '()] ',\
 $message, \"\\n\";\n    }\n}\n\n=head2 print_tool\
_params()\n\nPrint list of tool parameters.\n\n  &\
print_tool_params();\n\n=cut\n\nsub print_tool_par\
ams {\n    print_debug_message('print_tool_params'\
, 'Begin', 1);\n    my (@param_list) = &rest_get_p\
arameters();\n    foreach my $param (sort (@param_\
list)) {\n        print $param, \"\\n\";\n    }\n \
   print_debug_message('print_tool_params', 'End',\
 1);\n}\n\n=head2 print_param_details()\n\nPrint d\
etails of a tool parameter.\n\n  &print_param_deta\
ils($param_name);\n\n=cut\n\nsub print_param_detai\
ls {\n    print_debug_message('print_param_details\
', 'Begin', 1);\n    my $paramName = shift;\n    p\
rint_debug_message('print_param_details', 'paramNa\
me: ' . $paramName, 2);\n    my $paramDetail = &re\
st_get_parameter_details($paramName);\n    print $\
paramDetail->{'name'}, \"\\t\", $paramDetail->{'ty\
pe'}, \"\\n\";\n    print $paramDetail->{'descript\
ion'}, \"\\n\";\n    if (defined($paramDetail->{'v\
alues'}->{'value'})) {\n        if (ref($paramDeta\
il->{'values'}->{'value'}) eq 'ARRAY') {\n        \
    foreach my $value (@{$paramDetail->{'values'}-\
>{'value'}}) {\n                &print_param_value\
($value);\n            }\n        }\n        else \
{\n            &print_param_value($paramDetail->{'\
values'}->{'value'});\n        }\n    }\n    print\
_debug_message('print_param_details', 'End', 1);\n\
}\n\n=head2 print_param_value()\n\nPrint details o\
f a tool parameter value.\n\n  &print_param_detail\
s($param_value);\n\nUsed by print_param_details() \
to handle both singluar and array values.\n\n=cut\\
n\nsub print_param_value {\n    my $value = shift;\
\n    print $value->{'value'};\n    if ($value->{'\
defaultValue'} eq 'true') {\n        print \"\\t\"\
, 'default';\n    }\n    print \"\\n\";\n    print\
 \"\\t\", $value->{'label'}, \"\\n\";\n    if (def\
ined($value->{'properties'})) {\n        foreach\n\
        my $key (sort ( keys(%{$value->{'propertie\
s'}{'property'}}) )) {\n            if (ref($value\
->{'properties'}{'property'}{$key}) eq 'HASH'\n   \
             && defined($value->{'properties'}{'pr\
operty'}{$key}{'value'})\n            ) {\n       \
         print \"\\t\", $key, \"\\t\",\n          \
          $value->{'properties'}{'property'}{$key}\
{'value'}, \"\\n\";\n            }\n            el\
se {\n                print \"\\t\", $value->{'pro\
perties'}{'property'}{'key'},\n                   \
 \"\\t\", $value->{'properties'}{'property'}{'valu\
e'}, \"\\n\";\n                last;\n            \
}\n        }\n    }\n}\n\n=head2 print_job_status(\
)\n\nPrint status of a job.\n\n  &print_job_status\
($job_id);\n\n=cut\n\nsub print_job_status {\n    \
print_debug_message('print_job_status', 'Begin', 1\
);\n    my $jobid = shift;\n    print_debug_messag\
e('print_job_status', 'jobid: ' . $jobid, 1);\n   \
 if ($outputLevel > 0) {\n        print STDERR 'Ge\
tting status for job ', $jobid, \"\\n\";\n    }\n \
   my $result = &rest_get_status($jobid);\n    pri\
nt \"$result\\n\";\n    if ($result eq 'FINISHED' \
&& $outputLevel > 0) {\n        print STDERR \"To \
get results: perl $scriptName --polljob --jobid \"\
 . $jobid\n            . \"\\n\";\n    }\n    prin\
t_debug_message('print_job_status', 'End', 1);\n}\\
n\n=head2 print_result_types()\n\nPrint available \
result types for a job.\n\n  &print_result_types($\
job_id);\n\n=cut\n\nsub print_result_types {\n    \
print_debug_message('result_types', 'Begin', 1);\n\
    my $jobid = shift;\n    print_debug_message('r\
esult_types', 'jobid: ' . $jobid, 1);\n    if ($ou\
tputLevel > 0) {\n        print STDERR 'Getting re\
sult types for job ', $jobid, \"\\n\";\n    }\n   \
 my $status = &rest_get_status($jobid);\n    if ($\
status eq 'PENDING' || $status eq 'RUNNING') {\n  \
      print STDERR 'Error: Job status is ', $statu\
s,\n            '. To get result types the job mus\
t be finished.', \"\\n\";\n    }\n    else {\n    \
    my (@resultTypes) = &rest_get_result_types($jo\
bid);\n        if ($outputLevel > 0) {\n          \
  print STDOUT 'Available result types:', \"\\n\";\
\n        }\n        foreach my $resultType (@resu\
ltTypes) {\n            print STDOUT $resultType->\
{'identifier'}, \"\\n\";\n            if (defined(\
$resultType->{'label'})) {\n                print \
STDOUT \"\\t\", $resultType->{'label'}, \"\\n\";\n\
            }\n            if (defined($resultType\
->{'description'})) {\n                print STDOU\
T \"\\t\", $resultType->{'description'}, \"\\n\";\\
n            }\n            if (defined($resultTyp\
e->{'mediaType'})) {\n                print STDOUT\
 \"\\t\", $resultType->{'mediaType'}, \"\\n\";\n  \
          }\n            if (defined($resultType->\
{'fileSuffix'})) {\n                print STDOUT \\
"\\t\", $resultType->{'fileSuffix'}, \"\\n\";\n   \
         }\n        }\n        if ($status eq 'FIN\
ISHED' && $outputLevel > 0) {\n            print S\
TDERR \"\\n\", 'To get results:', \"\\n\",\n      \
          \"  perl $scriptName --polljob --jobid \\
" . $params{'jobid'} . \"\\n\",\n                \\
"  perl $scriptName --polljob --outformat <type> -\
-jobid \"\n                    . $params{'jobid'} \
. \"\\n\";\n        }\n    }\n    print_debug_mess\
age('result_types', 'End', 1);\n}\n\n=head2 submit\
_job()\n\nSubmit a job to the service.\n\n  &submi\
t_job($seq);\n\n=cut\n\nsub submit_job {\n    prin\
t_debug_message('submit_job', 'Begin', 1);\n\n    \
# Set input sequence\n    $params{'sequence'} = sh\
ift;\n    my $seq_id = shift;\n\n    # Load parame\
ters\n    &load_params();\n\n    # Submit the job\\
n    my $jobid = &rest_run($params{'email'}, $para\
ms{'title'}, \\%params);\n\n    # Asynchronous sub\
mission.\n    if (defined($params{'asyncjob'})) {\\
n        print STDOUT $jobid, \"\\n\";\n        if\
 ($outputLevel > 0) {\n            print STDERR\n \
               \"To check status: perl $scriptName\
 --status --jobid $jobid\\n\";\n        }\n    }\n\
\n    # Simulate synchronous submission serial mod\
e.\n    else {\n        if ($outputLevel > 0) {\n \
           print STDERR \"JobId: $jobid\\n\";\n   \
     } else {\n            print STDERR \"$jobid\\\
n\";\n        }\n        usleep($checkInterval);\n\
        # Get results.\n        &get_results($jobi\
d, $seq_id);\n\n    }\n    print_debug_message('su\
bmit_job', 'End', 1);\n    return $jobid;\n}\n=hea\
d2 multi_submit_job()\n\nSubmit multiple jobs assu\
ming input is a collection of fasta formatted sequ\
ences.\n\n  &multi_submit_job();\n\n=cut\n\nsub mu\
lti_submit_job {\n    print_debug_message('multi_s\
ubmit_job', 'Begin', 1);\n    my (@filename_list) \
= ();\n\n    # Query sequence\n    if (defined($AR\
GV[0])) {                  # Bare option\n        \
if (-f $ARGV[0] || $ARGV[0] eq '-') { # File\n    \
        push(@filename_list, $ARGV[0]);\n        }\
\n        else {\n            warn 'Warning: Input\
 file \"' . $ARGV[0] . '\" does not exist';\n     \
   }\n    }\n    if ($params{'sequence'}) {       \
                               # Via --sequence\n \
       if (-f $params{'sequence'} || $params{'sequ\
ence'} eq '-') { # File\n            push(@filenam\
e_list, $params{'sequence'});\n        }\n        \
else {\n            warn 'Warning: Input file \"'\\
n                . $params{'sequence'}\n          \
      . '\" does not exist';\n        }\n    }\n\n\
    # Job identifier tracking for parallel executi\
on.\n    my @jobid_list = ();\n    my $job_number \
= 0;\n    $/ = '>';\n    foreach my $filename (@fi\
lename_list) {\n        my $INFILE;\n        if ($\
filename eq '-') { # STDIN.\n            open($INF\
ILE, '<-')\n                or die 'Error: unable \
to STDIN (' . $! . ')';\n        }\n        else {\
 # File.\n            open($INFILE, '<', $filename\
)\n                or die 'Error: unable to open f\
ile '\n                . $filename . ' ('\n       \
         . $! . ')';\n        }\n        while (<$\
INFILE>) {\n            my $seq = $_;\n           \
 $seq =~ s/>$//;\n            if ($seq =~ m/(\\S+)\
/) {\n                my $seq_id = $1;\n          \
      print STDERR \"Submitting job for: $seq_id\\\
n\"\n                    if ($outputLevel > 0);\n \
               $seq = '>' . $seq;\n               \
 &print_debug_message('multi_submit_job', $seq, 11\
);\n                $job_number++;\n              \
  my $job_id = &submit_job($seq, $seq_id);\n\n    \
            my $job_info_str = sprintf('%s %d %d',\
 $job_id, 0, $job_number);\n\n                push\
(@jobid_list, $job_info_str);\n            }\n\n  \
          # Parallel mode, wait for job(s) to fini\
sh to free slots.\n            while ($params{'max\
Jobs'} > 1\n                && scalar(@jobid_list)\
 >= $params{'maxJobs'}) {\n                &_job_l\
ist_poll(\\@jobid_list);\n                print_de\
bug_message('multi_submit_job',\n                 \
   'Remaining jobs: ' . scalar(@jobid_list), 1);\n\
            }\n        }\n        close $INFILE;\n\
    }\n\n    # Parallel mode, wait for remaining j\
obs to finish.\n    while ($params{'maxJobs'} > 1 \
&& scalar(@jobid_list) > 0) {\n        &_job_list_\
poll(\\@jobid_list);\n        print_debug_message(\
'multi_submit_job',\n            'Remaining jobs: \
' . scalar(@jobid_list), 1);\n    }\n    print_deb\
ug_message('multi_submit_job', 'End', 1);\n}\n\n\n\
=head2 _job_list_poll()\n\nPoll the status of a li\
st of jobs and fetch results for finished jobs.\n\\
n  while(scalar(@jobid_list) > 0) {\n    &_job_lis\
t_poll(\\@jobid_list);\n  }\n\n=cut\n\nsub _job_li\
st_poll {\n    print_debug_message('_job_list_poll\
', 'Begin', 1);\n    my $jobid_list = shift;\n    \
print_debug_message('_job_list_poll', 'Num jobs: '\
 . scalar(@$jobid_list),\n        11);\n\n    # Lo\
op though job Id list polling job status.\n    for\
 (my $jobNum = (scalar(@$jobid_list) - 1); $jobNum\
 > -1; $jobNum--) {\n        my ($jobid, $seq_id, \
$error_count, $job_number) =\n            split(/\\
\s+/, $jobid_list->[$jobNum]);\n        print_debu\
g_message('_job_list_poll', 'jobNum: ' . $jobNum, \
12);\n        print_debug_message('_job_list_poll'\
,\n            'Job info: ' . $jobid_list->[$jobNu\
m], 12);\n\n        # Get job status.\n        my \
$job_status = &rest_get_status($jobid);\n        p\
rint_debug_message('_job_list_poll', 'Status: ' . \
$job_status, 12);\n\n        # Fetch results and r\
emove finished/failed jobs from list.\n        if \
(\n            !(\n                $job_status eq \
'RUNNING'\n                    || $job_status eq '\
PENDING'\n                    || ($job_status eq '\
ERROR'\n                    && $error_count < $max\
ErrorStatusCount)\n            )\n        ) {\n   \
         if ($job_status eq 'ERROR' || $job_status\
 eq 'FAILED') {\n                print STDERR\n   \
                 \"Warning: job $jobid failed for \
sequence $job_number: $seq_id\\n\";\n            }\
\n            # Duplicated getting results.\n     \
       #&get_results($jobid, $seq_id);\n          \
  splice(@$jobid_list, $jobNum, 1);\n        }\n  \
      else {\n\n            # Update error count, \
increment for new error or clear old errors.\n    \
        if ($job_status eq 'ERROR') {\n           \
     $error_count++;\n            }\n            e\
lsif ($error_count > 0) {\n                $error_\
count--;\n            }\n\n            # Update jo\
b tracking info.\n            my $job_info_str = s\
printf('%s %s %d %d',\n                $jobid, $se\
q_id, $error_count, $job_number);\n            $jo\
bid_list->[$jobNum] = $job_info_str;\n        }\n \
   }\n    print_debug_message('_job_list_poll', 'N\
um jobs: ' . scalar(@$jobid_list),\n        11);\n\
    print_debug_message('_job_list_poll', 'End', 1\
);\n}\n\n=head2 list_file_submit_job()\n\nSubmit m\
ultiple jobs using a file containing a list of ent\
ry identifiers as\ninput.\n\n  &list_file_submit_j\
ob($list_filename)\n\n=cut\n\nsub list_file_submit\
_job {\n    print_debug_message('list_file_submit_\
job', 'Begin', 1);\n    my $filename = shift;\n\n \
   # Open the file of identifiers.\n    my $LISTFI\
LE;\n    if ($filename eq '-') { # STDIN.\n       \
 open($LISTFILE, '<-')\n            or die 'Error:\
 unable to STDIN (' . $! . ')';\n    }\n    else {\
 # File.\n        open($LISTFILE, '<', $filename)\\
n            or die 'Error: unable to open file ' \
. $filename . ' (' . $! . ')';\n    }\n\n    # Job\
 identifier tracking for parallel execution.\n    \
my @jobid_list = ();\n    my $job_number = 0;\n\n \
   # Iterate over identifiers, submitting each job\
\n    while (<$LISTFILE>) {\n        my $line = $_\
;\n        chomp($line);\n        if ($line ne '')\
 {\n            &print_debug_message('list_file_su\
bmit_job', 'line: ' . $line, 2);\n            if (\
$line =~ m/\\w:\\w/) {\n                # Check th\
is is an identifier\n                my $seq_id = \
$line;\n                print STDERR \"Submitting \
job for: $seq_id\\n\"\n                    if ($ou\
tputLevel > 0);\n                $job_number++;\n \
               my $job_id = &submit_job($seq_id, $\
seq_id);\n                my $job_info_str =\n    \
                sprintf('%s %s %d %d', $job_id, $s\
eq_id, 0, $job_number);\n                push(@job\
id_list, $job_info_str);\n            }\n         \
   else {\n                print STDERR\n         \
           \"Warning: line \\\"$line\\\" is not re\
cognised as an identifier\\n\";\n            }\n\n\
            # Parallel mode, wait for job(s) to fi\
nish to free slots.\n            while ($params{'m\
axJobs'} > 1\n                && scalar(@jobid_lis\
t) >= $params{'maxJobs'}) {\n                &_job\
_list_poll(\\@jobid_list);\n                print_\
debug_message('list_file_submit_job',\n           \
         'Remaining jobs: ' . scalar(@jobid_list),\
 1);\n            }\n        }\n    }\n    close $\
LISTFILE;\n\n    # Parallel mode, wait for remaini\
ng jobs to finish.\n    while ($params{'maxJobs'} \
> 1 && scalar(@jobid_list) > 0) {\n        &_job_l\
ist_poll(\\@jobid_list);\n        print_debug_mess\
age('list_file_submit_job',\n            'Remainin\
g jobs: ' . scalar(@jobid_list), 1);\n    }\n    p\
rint_debug_message('list_file_submit_job', 'End', \
1);\n}\n\n\n=head2 load_data()\n\nLoad sequence da\
ta from file or option specified on the command-li\
ne.\n\n  &load_data();\n\n=cut\n\nsub load_data {\\
n    print_debug_message('load_data', 'Begin', 1);\
\n    my $retSeq;\n\n    # Query sequence\n    if \
(defined($ARGV[0])) {                  # Bare opti\
on\n        if (-f $ARGV[0] || $ARGV[0] eq '-') { \
# File\n            $retSeq = &read_file($ARGV[0])\
;\n        }\n        else { # DB:ID or sequence\n\
            $retSeq = $ARGV[0];\n        }\n    }\\
n    if ($params{'sequence'}) {                   \
                   # Via --sequence\n        if (-\
f $params{'sequence'} || $params{'sequence'} eq '-\
') { # File\n            $retSeq = &read_file($par\
ams{'sequence'});\n        }\n        else { # DB:\
ID or sequence\n            $retSeq = $params{'seq\
uence'};\n        }\n    }\n    print_debug_messag\
e('load_data', 'End', 1);\n    return $retSeq;\n}\\
n\n=head2 load_params()\n\nLoad job parameters fro\
m command-line options.\n\n  &load_params();\n\n=c\
ut\n\nsub load_params {\n    print_debug_message('\
load_params', 'Begin', 1);\n\n    # Pass default v\
alues and fix bools (without default value)\n    i\
f ($params{'stype'} eq 'protein') {\n        if (!\
$params{'task'}) {\n            $params{'task'} = \
'blastp'\n        }\n    }\n    if ($params{'stype\
'} eq 'nucleotide') {\n        if (!$params{'task'\
}) {\n            $params{'task'} = 'blastn'\n    \
    }\n    }\n    if ($params{'stype'} eq 'vector'\
) {\n        if (!$params{'task'}) {\n            \
$params{'task'} = 'blastn'\n        }\n    }\n\n  \
  if ($params{'stype'} eq 'protein') {\n        if\
 (!$params{'matrix'}) {\n            $params{'matr\
ix'} = 'BLOSUM62'\n        }\n    }\n    if ($para\
ms{'stype'} eq 'nucleotide') {\n        if (!$para\
ms{'matrix'}) {\n            $params{'matrix'} = '\
NONE'\n        }\n    }\n    if ($params{'stype'} \
eq 'vector') {\n        if (!$params{'matrix'}) {\\
n            $params{'matrix'} = 'NONE'\n        }\
\n    }\n\n    if (!$params{'alignments'}) {\n    \
    $params{'alignments'} = 50\n    }\n\n    if (!\
$params{'scores'}) {\n        $params{'scores'} = \
50\n    }\n\n    if (!$params{'exp'}) {\n        $\
params{'exp'} = '10'\n    }\n\n    if (!$params{'d\
ropoff'}) {\n        $params{'dropoff'} = 0\n    }\
\n\n    if ($params{'stype'} eq 'nucleotide') {\n \
       if (!$params{'match_scores'}) {\n          \
  $params{'match_scores'} = '1,-3'\n        }\n   \
 }\n    if ($params{'stype'} eq 'vector') {\n     \
   if (!$params{'match_scores'}) {\n            $p\
arams{'match_scores'} = '1,-3'\n        }\n    }\n\
\n    if (!$params{'gapopen'}) {\n        $params{\
'gapopen'} = -1\n    }\n\n    if (!$params{'gapext\
'}) {\n        $params{'gapext'} = -1\n    }\n\n  \
  if ($params{'stype'} eq 'protein') {\n        if\
 (!$params{'filter'}) {\n            $params{'filt\
er'} = 'F'\n        }\n    }\n    if ($params{'sty\
pe'} eq 'nucleotide') {\n        if (!$params{'fil\
ter'}) {\n            $params{'filter'} = 'T'\n   \
     }\n    }\n    if ($params{'stype'} eq 'vector\
') {\n        if (!$params{'filter'}) {\n         \
   $params{'filter'} = 'T'\n        }\n    }\n\n  \
  if (!$params{'gapalign'}) {\n        $params{'ga\
palign'} = 'true'\n    }\n\n    if (!$params{'comp\
stats'}) {\n        $params{'compstats'} = 'F'\n  \
  }\n\n    if (!$params{'align'}) {\n        $para\
ms{'align'} = 0\n    }\n\n    if (!$params{'transl\
table'}) {\n        $params{'transltable'} = 1\n  \
  }\n\n    print_debug_message('load_params', 'End\
', 1);\n}\n\n=head2 client_poll()\n\nClient-side j\
ob polling.\n\n  &client_poll($job_id);\n\n=cut\n\\
nsub client_poll {\n    print_debug_message('clien\
t_poll', 'Begin', 1);\n    my $jobid = shift;\n   \
 my $status = 'PENDING';\n\n    # Check status and\
 wait if not finished. Terminate if three attempts\
 get \"ERROR\".\n    my $errorCount = 0;\n    whil\
e ($status eq 'RUNNING'\n        || $status eq 'PE\
NDING'\n        || ($status eq 'ERROR' && $errorCo\
unt < 2)) {\n        $status = rest_get_status($jo\
bid);\n        print STDERR \"$status\\n\" if ($ou\
tputLevel > 0);\n        if ($status eq 'ERROR') {\
\n            $errorCount++;\n        }\n        e\
lsif ($errorCount > 0) {\n            $errorCount-\
-;\n        }\n        if ($status eq 'RUNNING'\n \
           || $status eq 'PENDING'\n            ||\
 $status eq 'ERROR') {\n\n            # Wait befor\
e polling again.\n            usleep($checkInterva\
l);\n        }\n    }\n    print_debug_message('cl\
ient_poll', 'End', 1);\n    return $status;\n}\n\n\
=head2 get_results()\n\nGet the results for a job \
identifier.\n\n  &get_results($job_id);\n\n=cut\n\\
nsub get_results {\n    print_debug_message('get_r\
esults', 'Begin', 1);\n    my $jobid = shift;\n   \
 print_debug_message('get_results', 'jobid: ' . $j\
obid, 1);\n    my $seq_id = shift;\n    print_debu\
g_message('get_results', 'seq_id: ' . $seq_id, 1) \
if ($seq_id);\n\n    my $output_basename = $jobid;\
\n\n    # Verbose\n    if ($outputLevel > 1) {\n  \
      print 'Getting results for job ', $jobid, \"\
\\n\";\n    }\n\n    # Check status, and wait if n\
ot finished\n    client_poll($jobid);\n\n    # Def\
ault output file names use JobId, however the name\
 can be specified...\n    if (defined($params{'out\
file'})) {\n        $output_basename = $params{'ou\
tfile'};\n    }\n    # Or use sequence identifer.\\
n    elsif (defined($params{'useSeqId'} && defined\
($seq_id) && $seq_id ne '')) {\n        $output_ba\
sename = $seq_id;\n\n        # Make safe to use as\
 a file name.\n        $output_basename =~ s/\\W/_\
/g;\n    }\n\n    # Use JobId if output file name \
is not defined\n    else {\n        unless (define\
d($params{'outfile'})) {\n            #$params{'ou\
tfile'} = $jobid;\n            $output_basename = \
$jobid;\n        }\n    }\n\n    # Get list of dat\
a types\n    my (@resultTypes) = rest_get_result_t\
ypes($jobid);\n\n\n    # Get the data and write it\
 to a file\n    if (defined($params{'outformat'}))\
 {\n        # Specified data type\n        # check\
 to see if there are multiple formats (comma separ\
ated)\n        my $sep = \",\";\n        my (@mult\
ResultTypes);\n        if ($params{'outformat'} =~\
 /$sep/) {\n            @multResultTypes = split('\
,', $params{'outformat'});\n        }\n        els\
e {\n            $multResultTypes[0] = $params{'ou\
tformat'};\n        }\n        # check if the prov\
ided formats are recognised\n        foreach my $i\
nputType (@multResultTypes) {\n            my $exp\
ectation = 0;\n            foreach my $resultType \
(@resultTypes) {\n                if ($resultType-\
>{'identifier'} eq $inputType && $expectation eq 0\
) {\n                    $expectation = 1;\n      \
          }\n            }\n            if ($expec\
tation ne 1) {\n                die 'Error: unknow\
n result format \"' . $inputType . '\"';\n        \
    }\n        }\n        # if so get the files\n \
       my $selResultType;\n        foreach my $res\
ultType (@resultTypes) {\n            if (grep {$_\
 eq $resultType->{'identifier'}} @multResultTypes)\
 {\n                $selResultType = $resultType;\\
n                my $result = rest_get_result($job\
id, $selResultType->{'identifier'});\n            \
    if (defined($params{'outfile'}) && $params{'ou\
tfile'} eq '-') {\n                    write_file(\
$params{'outfile'}, $result);\n                }\n\
                else {\n                    write_\
file(\n                        $output_basename . \
'.'\n                            . $selResultType-\
>{'identifier'} . '.'\n                           \
 . $selResultType->{'fileSuffix'},\n              \
          $result\n                    );\n       \
         }\n            }\n        }\n    }\n    e\
lse { # Data types available\n        # Write a fi\
le for each output type\n        for my $resultTyp\
e (@resultTypes) {\n            if ($outputLevel >\
 1) {\n                print STDERR 'Getting ', $r\
esultType->{'identifier'}, \"\\n\";\n            }\
\n            my $result = rest_get_result($jobid,\
 $resultType->{'identifier'});\n            if (de\
fined($params{'outfile'}) && $params{'outfile'} eq\
 '-') {\n                write_file($params{'outfi\
le'}, $result);\n            }\n            else {\
\n                write_file(\n                   \
 $output_basename . '.'\n                        .\
 $resultType->{'identifier'} . '.'\n              \
          . $resultType->{'fileSuffix'},\n        \
            $result\n                );\n         \
   }\n        }\n    }\n    print_debug_message('g\
et_results', 'End', 1);\n}\n\n=head2 read_file()\n\
\nRead a file into a scalar. The special filename \
'-' can be used to read from\nstandard input (STDI\
N).\n\n  my $data = &read_file($filename);\n\n=cut\
\n\nsub read_file {\n    print_debug_message('read\
_file', 'Begin', 1);\n    my $filename = shift;\n \
   print_debug_message('read_file', 'filename: ' .\
 $filename, 2);\n    my ($content, $buffer);\n    \
if ($filename eq '-') {\n        while (sysread(ST\
DIN, $buffer, 1024)) {\n            $content .= $b\
uffer;\n        }\n    }\n    else {\n        # Fi\
le\n        open(my $FILE, '<', $filename)\n      \
      or die \"Error: unable to open input file $f\
ilename ($!)\";\n        while (sysread($FILE, $bu\
ffer, 1024)) {\n            $content .= $buffer;\n\
        }\n        close($FILE);\n    }\n    print\
_debug_message('read_file', 'End', 1);\n    return\
 $content;\n}\n\n=head2 write_file()\n\nWrite data\
 to a file. The special filename '-' can be used t\
o write to\nstandard output (STDOUT).\n\n  &write_\
file($filename, $data);\n\n=cut\n\nsub write_file \
{\n    print_debug_message('write_file', 'Begin', \
1);\n    my ($filename, $data) = @_;\n    print_de\
bug_message('write_file', 'filename: ' . $filename\
, 2);\n    if ($outputLevel > 0) {\n        print \
STDERR 'Creating result file: ' . $filename . \"\\\
n\";\n    }\n    if ($filename eq '-') {\n        \
print STDOUT $data;\n    }\n    else {\n        op\
en(my $FILE, '>', $filename)\n            or die \\
"Error: unable to open output file $filename ($!)\\
";\n        syswrite($FILE, $data);\n        close\
($FILE);\n    }\n    print_debug_message('write_fi\
le', 'End', 1);\n}\n\n=head2 usage()\n\nPrint prog\
ram usage message.\n\n  &usage();\n\n=cut\n\nsub u\
sage {\n    print STDERR <<EOF\nEMBL-EBI NCBI Blas\
t Perl Client:\n\nSequence similarity search with \
NCBI Blast.\n\n[Required (for job submission)]\n  \
--email               E-mail address.\n  --program\
             The BLAST program to be used for the \
Sequence Similarity\n                        Searc\
h.\n  --stype               Indicates if the seque\
nce is protein or DNA/RNA.\n  --sequence          \
  The query sequence can be entered directly into \
this form.\n                        The sequence c\
an be in GCG, FASTA, EMBL (Nucleotide only),\n    \
                    GenBank, PIR, NBRF, PHYLIP or \
UniProtKB/Swiss-Prot (Protein\n                   \
     only) format. A partially formatted sequence \
is not\n                        accepted. Adding a\
 return to the end of the sequence may\n          \
              help certain applications understand\
 the input. Note that\n                        dir\
ectly using data from word processors may yield\n \
                       unpredictable results as hi\
dden/control characters may be\n                  \
      present.\n  --database            Database.\\
n\n[Optional]\n  --task                Task option\
 (only selectable for blastn).\n  --matrix        \
      (Protein searches) The substitution matrix u\
sed for scoring\n                        alignment\
s when searching the database.\n  --alignments    \
      Maximum number of match alignments reported \
in the result\n                        output.\n  \
--scores              Maximum number of match scor\
e summaries reported in the\n                     \
   result output.\n  --exp                 Limits \
the number of scores and alignments reported based\
 on\n                        the expectation value\
. This is the maximum number of times\n           \
             the match is expected to occur by cha\
nce.\n  --dropoff             The amount a score c\
an drop before gapped extension of word\n         \
               hits is halted.\n  --match_scores  \
      (Nucleotide searches) The match score is the\
 bonus to the\n                        alignment s\
core when matching the same base. The mismatch is\\
n                        the penalty when failing \
to match.\n  --gapopen             Penalty taken a\
way from the score when a gap is created in\n     \
                   sequence. Increasing the gap op\
enning penalty will decrease\n                    \
    the number of gaps in the final alignment.\n  \
--gapext              Penalty taken away from the \
score for each base or residue\n                  \
      in the gap. Increasing the gap extension pen\
alty favors\n                        short gaps in\
 the final alignment, conversly decreasing the\n  \
                      gap extension penalty favors\
 long gaps in the final\n                        a\
lignment.\n  --filter              Filter regions \
of low sequence complexity. This can avoid\n      \
                  issues with low complexity seque\
nces where matches are found\n                    \
    due to composition rather than meaningful sequ\
ence\n                        similarity. However \
in some cases filtering also masks\n              \
          regions of interest and so should be use\
d with caution.\n  --seqrange            Specify a\
 range or section of the input sequence to use in\\
n                        the search. Example: Spec\
ifying '34-89' in an input sequence\n             \
           of total length 100, will tell BLAST to\
 only use residues 34\n                        to \
89, inclusive.\n  --gapalign            This is a \
true/false setting that tells the program the\n   \
                     perform optimised alignments \
within regions involving gaps.\n                  \
      If set to true, the program will perform an \
alignment using\n                        gaps. Oth\
erwise, if it is set to false, it will report only\
\n                        individual HSP where two\
 sequence match each other, and thus\n            \
            will not produce alignments with gaps.\
\n  --wordsize            Word size for wordfinder\
 algorithm.\n  --taxids              Specify one o\
r more TaxIDs so that the BLAST search becomes\n  \
                      taxonomically aware.\n  --co\
mpstats           Use composition-based statistics\
.\n  --align               Formating for the align\
ments.\n  --transltable         Query Genetic code\
 to use in translation.\n\n[General]\n  -h, --help\
            Show this help message and exit.\n  --\
asyncjob            Forces to make an asynchronous\
 query.\n  --title               Title for job.\n \
 --status              Get job status.\n  --result\
Types         Get available result types for job.\\
n  --polljob             Poll for the status of a \
job.\n  --pollFreq            Poll frequency in se\
conds (default 3s).\n  --jobid               JobId\
 that was returned when an asynchronous job was su\
bmitted.\n  --outfile             File name for re\
sults (default is JobId; for STDOUT).\n  --multifa\
sta          Treat input as a set of fasta formatt\
ed sequences.\n  --useSeqId            Use sequenc\
e identifiers for output filenames.\n             \
           Only available in multi-fasta and multi\
-identifier modes.\n  --maxJobs             Maximu\
m number of concurrent jobs. Only\n               \
         available in multifasta or list file mode\
s.\n  --outformat           Result format(s) to re\
trieve. It accepts comma-separated values.\n  --pa\
rams              List input parameters.\n  --para\
mDetail         Display details for input paramete\
r.\n  --quiet               Decrease output.\n  --\
verbose             Increase output.\n  --version \
            Prints out the version of the Client a\
nd exit.\n  --baseUrl             Base URL. Defaul\
ts to:\n                        https://www.ebi.ac\
.uk/Tools/services/rest/ncbiblast\n\nSynchronous j\
ob:\n  The results/errors are returned as soon as \
the job is finished.\n  Usage: perl $scriptName --\
email <your\\@email.com> [options...] <SeqFile|Seq\
ID(s)>\n  Returns: results as an attachment\n\nAsy\
nchronous job:\n  Use this if you want to retrieve\
 the results at a later time. The results\n  are s\
tored for up to 24 hours.\n  Usage: perl $scriptNa\
me --asyncjob --email <your\\@email.com> [options.\
..] <SeqFile|SeqID(s)>\n  Returns: jobid\n\nCheck \
status of Asynchronous job:\n  Usage: perl $script\
Name --status --jobid <jobId>\n\nRetrieve job data\
:\n  Use the jobid to query for the status of the \
job. If the job is finished,\n  it also returns th\
e results/errors.\n  Usage: perl $scriptName --pol\
ljob --jobid <jobId> [--outfile string]\n  Returns\
: string indicating the status of the job and if a\
pplicable, results\n  as an attachment.\n\nFurther\
 information:\n  https://www.ebi.ac.uk/Tools/webse\
rvices and\n    https://github.com/ebi-wp/webservi\
ce-clients\n\nSupport/Feedback:\n  https://www.ebi\
.ac.uk/support/\nEOF\n}\n\n=head1 FEEDBACK/SUPPORT\
\n\nPlease contact us at L<https://www.ebi.ac.uk/s\
upport/> if you have any\nfeedback, suggestions or\
 issues with the service or this client.\n\n=cut\n\
","\n=head1 NAME\n\nwublast_lwp.pl\n\n=head1 DESCR\
IPTION\n\nWU-BLAST (REST) web service Perl client \
using L<LWP>.\n\nTested with:\n\n=over\n\n=item *\\
nL<LWP> 5.79, L<XML::Simple> 2.12 and Perl 5.8.3\n\
\n=item *\nL<LWP> 5.808, L<XML::Simple> 2.18 and P\
erl 5.8.8 (Ubuntu 8.04 LTS)\n\n=item *\nL<LWP> 5.8\
34, L<XML::Simple> 2.18 and Perl 5.10.1 (Ubuntu 10\
.04 LTS)\n\n=item *\nL<LWP> 6.03, L<XML::Simple> 2\
.18 and Perl 5.14.2 (Ubuntu 12.04 LTS)\n\n=back\n\\
nFor further information see:\n\n=over\n\n=item *\\
nL<http://www.ebi.ac.uk/Tools/webservices/services\
/sss/wu_blast_rest>\n\n=item *\nL<http://www.ebi.a\
c.uk/Tools/webservices/tutorials/perl>\n\n=back\n\\
n=head1 LICENSE\n\nCopyright 2012-2013 EMBL - Euro\
pean Bioinformatics Institute\n\nLicensed under th\
e Apache License, Version 2.0 (the \"License\");\n\
you may not use this file except in compliance wit\
h the License.\nYou may obtain a copy of the Licen\
se at\n\n    http://www.apache.org/licenses/LICENS\
E-2.0\n\nUnless required by applicable law or agre\
ed to in writing, software\ndistributed under the \
License is distributed on an \"AS IS\" BASIS,\nWIT\
HOUT WARRANTIES OR CONDITIONS OF ANY KIND, either \
express or implied.\nSee the License for the speci\
fic language governing permissions and\nlimitation\
s under the License.\n\n=head1 VERSION\n\n$Id: wub\
last_lwp.pl 2560 2013-03-20 12:56:31Z hpm $\n\n=cu\
t\n\nuse strict;\nuse warnings;\n\nuse English;\nu\
se LWP;\nuse XML::Simple;\nuse Getopt::Long qw(:co\
nfig no_ignore_case bundling);\nuse File::Basename\
;\nuse Data::Dumper;\n\nmy $baseUrl = 'http://www.\
ebi.ac.uk/Tools/services/rest/wublast';\n\nmy $che\
ckInterval = 3;\n\nmy $outputLevel = 1;\n\nmy $num\
Opts = scalar(@ARGV);\nmy %params = ( 'debugLevel'\
 => 0 );\n\nmy %tool_params = ();\nGetOptions(\n\n\
	# Tool specific options\n	'program|p=s'     => \\\
$tool_params{'program'},      # BLAST program\n	'd\
atabase|D=s'    => \\$params{'database'},     # Se\
arch database\n	'matrix|m=s'      => \\$tool_param\
s{'matrix'},       # Scoring matrix\n	'exp|E=f'   \
      => \\$tool_params{'exp'},          # E-value\
 threshold\n	'viewfilter|e'    => \\$tool_params{'\
viewfilter'},   # Display filtered sequence\n	'fil\
ter|f=s'      => \\$tool_params{'filter'},       #\
 Low complexity filter name\n	'alignments|n=i'  =>\
 \\$tool_params{'alignments'},   # Number of align\
ments\n	'scores|s=i'      => \\$tool_params{'score\
s'},       # Number of scores\n	'sensitivity|S=s' \
=> \\$tool_params{'sensitivity'},  # Search sensit\
ivity\n	'sort|t=s'        => \\$tool_params{'sort'\
},         # Sort hits by...\n	'stats|T=s'       =\
> \\$tool_params{'stats'},        # Scoring statis\
tic to use\n	'strand|d=s'      => \\$tool_params{'\
strand'},       # Strand to use\n	'topcombon|c=i' \
  => \\$tool_params{'topcombon'},    # Consistent \
sets of HSPs\n	'align|A=i'       => \\$tool_params\
{'align'},   # Pairwise alignment format\n	'stype=\
s' => \\$tool_params{'stype'},    # Sequence type \
'protein' or 'dna'\n	'sequence=s' => \\$params{'se\
quence'},         # Query sequence file or DB:ID\n\
	'multifasta' => \\$params{'multifasta'},       # \
Multiple fasta input\n\n	# Compatability options, \
old command-line.\n	'echofilter|e'    => \\$params\
{'echofilter'},   # Display filtered sequence\n	'b\
=i'  => \\$params{'numal'},        # Number of ali\
gnments\n	'appxml=s'        => \\$params{'appxml'}\
,       # Application XML\n\n	# Generic options\n	\
'email=s'       => \\$params{'email'},          # \
User e-mail address\n	'title=s'       => \\$params\
{'title'},          # Job title\n	'outfile=s'     \
=> \\$params{'outfile'},        # Output file name\
\n	'outformat=s'   => \\$params{'outformat'},     \
 # Output file type\n	'jobid=s'       => \\$params\
{'jobid'},          # JobId\n	'help|h'        => \\
\$params{'help'},           # Usage help\n	'async'\
         => \\$params{'async'},          # Asynchr\
onous submission\n	'polljob'       => \\$params{'p\
olljob'},        # Get results\n	'resultTypes'   =\
> \\$params{'resultTypes'},    # Get result types\\
n	'status'        => \\$params{'status'},         \
# Get status\n	'params'        => \\$params{'param\
s'},         # List input parameters\n	'paramDetai\
l=s' => \\$params{'paramDetail'},    # Get details\
 for parameter\n	'quiet'         => \\$params{'qui\
et'},          # Decrease output level\n	'verbose'\
       => \\$params{'verbose'},        # Increase \
output level\n	'debugLevel=i'  => \\$params{'debug\
Level'},     # Debug output level\n	'baseUrl=s'   \
  => \\$baseUrl,                  # Base URL for s\
ervice.\n);\nif ( $params{'verbose'} ) { $outputLe\
vel++ }\nif ( $params{'quiet'} )  { $outputLevel--\
 }\n\n&print_debug_message( 'MAIN', 'LWP::VERSION:\
 ' . $LWP::VERSION,\n	1 );\n\n&print_debug_message\
( 'MAIN', \"params:\\n\" . Dumper( \\%params ),   \
        11 );\n&print_debug_message( 'MAIN', \"too\
l_params:\\n\" . Dumper( \\%tool_params ), 11 );\n\
\nmy $ua;\n\nmy $scriptName = basename( $0, () );\\
n\nif ( $params{'help'} || $numOpts == 0 ) {\n	&us\
age();\n	exit(0);\n}\n\n&print_debug_message( 'MAI\
N', 'baseUrl: ' . $baseUrl, 1 );\n\nif (\n	!(\n		 \
  $params{'polljob'}\n		|| $params{'resultTypes'}\\
n		|| $params{'status'}\n		|| $params{'params'}\n	\
	|| $params{'paramDetail'}\n	)\n	&& !( defined( $A\
RGV[0] ) || defined( $params{'sequence'} ) )\n  )\\
n{\n\n	# Bad argument combination, so print error \
message and usage\n	print STDERR 'Error: bad optio\
n combination', \"\\n\";\n	&usage();\n	exit(1);\n}\
\n\nelsif ( $params{'params'} ) {\n	&print_tool_pa\
rams();\n}\n\nelsif ( $params{'paramDetail'} ) {\n\
	&print_param_details( $params{'paramDetail'} );\n\
}\n\nelsif ( $params{'status'} && defined( $params\
{'jobid'} ) ) {\n	&print_job_status( $params{'jobi\
d'} );\n}\n\nelsif ( $params{'resultTypes'} && def\
ined( $params{'jobid'} ) ) {\n	&print_result_types\
( $params{'jobid'} );\n}\n\nelsif ( $params{'pollj\
ob'} && defined( $params{'jobid'} ) ) {\n	&get_res\
ults( $params{'jobid'} );\n}\n\nelse {\n\n	# Multi\
ple input sequence mode, assume fasta format.\n	if\
 ( $params{'multifasta'} ) {\n		&multi_submit_job(\
);\n	}\n\n	# Entry identifier list file.\n	elsif (\
( defined( $params{'sequence'} ) && $params{'seque\
nce'} =~ m/^\\@/ )\n		|| ( defined( $ARGV[0] ) && \
$ARGV[0] =~ m/^\\@/ ) )\n	{\n		my $list_filename =\
 $params{'sequence'} || $ARGV[0];\n		$list_filenam\
e =~ s/^\\@//;\n		&list_file_submit_job($list_file\
name);\n	}\n\n	# Default: single sequence/identifi\
er.\n	else {\n\n		# Load the sequence data and sub\
mit.\n		&submit_job( &load_data() );\n	}\n}\n\n=he\
ad1 FUNCTIONS\n\n=cut\n\n\n=head2 rest_user_agent(\
)\n\nGet a LWP UserAgent to use to perform REST re\
quests.\n\n  my $ua = &rest_user_agent();\n\n=cut\\
n\nsub rest_user_agent() {\n	print_debug_message( \
'rest_user_agent', 'Begin', 21 );\n	# Create an LW\
P UserAgent for making HTTP calls.\n	my $ua = LWP:\
:UserAgent->new();\n	# Set 'User-Agent' HTTP heade\
r to identifiy the client.\n	'$Revision: 2560 $' =\
~ m/(\\d+)/;\n	$ua->agent(\"EBI-Sample-Client/$1 (\
$scriptName; $OSNAME) \" . $ua->agent());\n	# Conf\
igure HTTP proxy support from environment.\n	$ua->\
env_proxy;\n	print_debug_message( 'rest_user_agent\
', 'End', 21 );\n	return $ua;\n}\n\n=head2 rest_er\
ror()\n\nCheck a REST response for an error condit\
ion. An error is mapped to a die.\n\n  &rest_error\
($response, $content_data);\n\n=cut\n\nsub rest_er\
ror() {\n	print_debug_message( 'rest_error', 'Begi\
n', 21 );\n	my $response = shift;\n	my $contentdat\
a;\n	if(scalar(@_) > 0) {\n		$contentdata = shift;\
\n	}\n	if(!defined($contentdata) || $contentdata e\
q '') {\n		$contentdata = $response->content();\n	\
}\n	# Check for HTTP error codes\n	if ( $response-\
>is_error ) {\n		my $error_message = '';\n		# HTML\
 response.\n		if(	$contentdata =~ m/<h1>([^<]+)<\\\
/h1>/ ) {\n			$error_message = $1;\n		}\n		#  XML \
response.\n		elsif($contentdata =~ m/<description>\
([^<]+)<\\/description>/) {\n			$error_message = $\
1;\n		}\n		die 'http status: ' . $response->code .\
 ' ' . $response->message . '  ' . $error_message;\
\n	}\n	print_debug_message( 'rest_error', 'End', 2\
1 );\n}\n\n=head2 rest_request()\n\nPerform a REST\
 request (HTTP GET).\n\n  my $response_str = &rest\
_request($url);\n\n=cut\n\nsub rest_request {\n	pr\
int_debug_message( 'rest_request', 'Begin', 11 );\\
n	my $requestUrl = shift;\n	print_debug_message( '\
rest_request', 'URL: ' . $requestUrl, 11 );\n\n	# \
Get an LWP UserAgent.\n	$ua = &rest_user_agent() u\
nless defined($ua);\n	# Available HTTP compression\
 methods.\n	my $can_accept;\n	eval {\n	    $can_ac\
cept = HTTP::Message::decodable();\n	};\n	$can_acc\
ept = '' unless defined($can_accept);\n	# Perform \
the request\n	my $response = $ua->get($requestUrl,\
\n		'Accept-Encoding' => $can_accept, # HTTP compr\
ession.\n	);\n	print_debug_message( 'rest_request'\
, 'HTTP status: ' . $response->code,\n		11 );\n	pr\
int_debug_message( 'rest_request',\n		'response le\
ngth: ' . length($response->content()), 11 );\n	pr\
int_debug_message( 'rest_request',\n		'request:' .\
\"\\n\" . $response->request()->as_string(), 32 );\
\n	print_debug_message( 'rest_request',\n		'respon\
se: ' . \"\\n\" . $response->as_string(), 32 );\n	\
# Unpack possibly compressed response.\n	my $retVa\
l;\n	if ( defined($can_accept) && $can_accept ne '\
') {\n	    $retVal = $response->decoded_content();\
\n	}\n	# If unable to decode use orginal content.\\
n	$retVal = $response->content() unless defined($r\
etVal);\n	# Check for an error.\n	&rest_error($res\
ponse, $retVal);\n	print_debug_message( 'rest_requ\
est', 'retVal: ' . $retVal, 12 );\n	print_debug_me\
ssage( 'rest_request', 'End', 11 );\n\n	# Return t\
he response data\n	return $retVal;\n}\n\n=head2 re\
st_get_parameters()\n\nGet list of tool parameter \
names.\n\n  my (@param_list) = &rest_get_parameter\
s();\n\n=cut\n\nsub rest_get_parameters {\n	print_\
debug_message( 'rest_get_parameters', 'Begin', 1 )\
;\n	my $url                = $baseUrl . '/paramete\
rs/';\n	my $param_list_xml_str = rest_request($url\
);\n	my $param_list_xml     = XMLin($param_list_xm\
l_str);\n	my (@param_list)       = @{ $param_list_\
xml->{'id'} };\n	print_debug_message( 'rest_get_pa\
rameters', 'End', 1 );\n	return (@param_list);\n}\\
n\n=head2 rest_get_parameter_details()\n\nGet deta\
ils of a tool parameter.\n\n  my $paramDetail = &r\
est_get_parameter_details($param_name);\n\n=cut\n\\
nsub rest_get_parameter_details {\n	print_debug_me\
ssage( 'rest_get_parameter_details', 'Begin', 1 );\
\n	my $parameterId = shift;\n	print_debug_message(\
 'rest_get_parameter_details',\n		'parameterId: ' \
. $parameterId, 1 );\n	my $url                  = \
$baseUrl . '/parameterdetails/' . $parameterId;\n	\
my $param_detail_xml_str = rest_request($url);\n	m\
y $param_detail_xml     = XMLin($param_detail_xml_\
str);\n	print_debug_message( 'rest_get_parameter_d\
etails', 'End', 1 );\n	return ($param_detail_xml);\
\n}\n\n=head2 rest_run()\n\nSubmit a job.\n\n  my \
$job_id = &rest_run($email, $title, \\%params );\n\
\n=cut\n\nsub rest_run {\n	print_debug_message( 'r\
est_run', 'Begin', 1 );\n	my $email  = shift;\n	my\
 $title  = shift;\n	my $params = shift;\n	print_de\
bug_message( 'rest_run', 'email: ' . $email, 1 );\\
n	if ( defined($title) ) {\n		print_debug_message(\
 'rest_run', 'title: ' . $title, 1 );\n	}\n	print_\
debug_message( 'rest_run', 'params: ' . Dumper($pa\
rams), 1 );\n\n	# Get an LWP UserAgent.\n	$ua = &r\
est_user_agent() unless defined($ua);\n\n	# Clean \
up parameters\n	my (%tmp_params) = %{$params};\n	$\
tmp_params{'email'} = $email;\n	$tmp_params{'title\
'} = $title;\n	foreach my $param_name ( keys(%tmp_\
params) ) {\n		if ( !defined( $tmp_params{$param_n\
ame} ) ) {\n			delete $tmp_params{$param_name};\n	\
	}\n	}\n\n	# Submit the job as a POST\n	my $url = \
$baseUrl . '/run';\n	my $response = $ua->post( $ur\
l, \\%tmp_params );\n	print_debug_message( 'rest_r\
un', 'HTTP status: ' . $response->code, 11 );\n	pr\
int_debug_message( 'rest_run',\n		'request:' .\"\\\
n\" . $response->request()->as_string(), 11 );\n	p\
rint_debug_message( 'rest_run',\n		'response: ' . \
length($response->as_string()) . \"\\n\" . $respon\
se->as_string(), 11 );\n\n	# Check for an error.\n\
	&rest_error($response);\n\n	# The job id is retur\
ned\n	my $job_id = $response->content();\n	print_d\
ebug_message( 'rest_run', 'End', 1 );\n	return $jo\
b_id;\n}\n\n=head2 rest_get_status()\n\nCheck the \
status of a job.\n\n  my $status = &rest_get_statu\
s($job_id);\n\n=cut\n\nsub rest_get_status {\n	pri\
nt_debug_message( 'rest_get_status', 'Begin', 1 );\
\n	my $job_id = shift;\n	print_debug_message( 'res\
t_get_status', 'jobid: ' . $job_id, 2 );\n	my $sta\
tus_str = 'UNKNOWN';\n	my $url        = $baseUrl .\
 '/status/' . $job_id;\n	$status_str = &rest_reque\
st($url);\n	print_debug_message( 'rest_get_status'\
, 'status_str: ' . $status_str, 2 );\n	print_debug\
_message( 'rest_get_status', 'End', 1 );\n	return \
$status_str;\n}\n\n=head2 rest_get_result_types()\\
n\nGet list of result types for finished job.\n\n \
 my (@result_types) = &rest_get_result_types($job_\
id);\n\n=cut\n\nsub rest_get_result_types {\n	prin\
t_debug_message( 'rest_get_result_types', 'Begin',\
 1 );\n	my $job_id = shift;\n	print_debug_message(\
 'rest_get_result_types', 'jobid: ' . $job_id, 2 )\
;\n	my (@resultTypes);\n	my $url                  \
    = $baseUrl . '/resulttypes/' . $job_id;\n	my $\
result_type_list_xml_str = &rest_request($url);\n	\
my $result_type_list_xml     = XMLin($result_type_\
list_xml_str);\n	(@resultTypes) = @{ $result_type_\
list_xml->{'type'} };\n	print_debug_message( 'rest\
_get_result_types',\n		scalar(@resultTypes) . ' re\
sult types', 2 );\n	print_debug_message( 'rest_get\
_result_types', 'End', 1 );\n	return (@resultTypes\
);\n}\n\n=head2 rest_get_result()\n\nGet result da\
ta of a specified type for a finished job.\n\n  my\
 $result = rest_get_result($job_id, $result_type);\
\n\n=cut\n\nsub rest_get_result {\n	print_debug_me\
ssage( 'rest_get_result', 'Begin', 1 );\n	my $job_\
id = shift;\n	my $type   = shift;\n	print_debug_me\
ssage( 'rest_get_result', 'jobid: ' . $job_id, 1 )\
;\n	print_debug_message( 'rest_get_result', 'type:\
 ' . $type,    1 );\n	my $url    = $baseUrl . '/re\
sult/' . $job_id . '/' . $type;\n	my $result = &re\
st_request($url);\n	print_debug_message( 'rest_get\
_result', length($result) . ' characters',\n		1 );\
\n	print_debug_message( 'rest_get_result', 'End', \
1 );\n	return $result;\n}\n\n\n=head2 print_debug_\
message()\n\nPrint debug message at specified debu\
g level.\n\n  &print_debug_message($method_name, $\
message, $level);\n\n=cut\n\nsub print_debug_messa\
ge {\n	my $function_name = shift;\n	my $message   \
    = shift;\n	my $level         = shift;\n	if ( $\
level <= $params{'debugLevel'} ) {\n		print STDERR\
 '[', $function_name, '()] ', $message, \"\\n\";\n\
	}\n}\n\n=head2 print_tool_params()\n\nPrint list \
of tool parameters.\n\n  &print_tool_params();\n\n\
=cut\n\nsub print_tool_params {\n	print_debug_mess\
age( 'print_tool_params', 'Begin', 1 );\n	my (@par\
am_list) = &rest_get_parameters();\n	foreach my $p\
aram ( sort(@param_list) ) {\n		print $param, \"\\\
n\";\n	}\n	print_debug_message( 'print_tool_params\
', 'End', 1 );\n}\n\n=head2 print_param_details()\\
n\nPrint details of a tool parameter.\n\n  &print_\
param_details($param_name);\n\n=cut\n\nsub print_p\
aram_details {\n	print_debug_message( 'print_param\
_details', 'Begin', 1 );\n	my $paramName = shift;\\
n	print_debug_message( 'print_param_details', 'par\
amName: ' . $paramName, 2 );\n	my $paramDetail = &\
rest_get_parameter_details($paramName);\n	print $p\
aramDetail->{'name'}, \"\\t\", $paramDetail->{'typ\
e'}, \"\\n\";\n	print $paramDetail->{'description'\
}, \"\\n\";\n	if(defined($paramDetail->{'values'}-\
>{'value'})) {\n		if(ref($paramDetail->{'values'}-\
>{'value'}) eq 'ARRAY') {\n			foreach my $value ( \
@{ $paramDetail->{'values'}->{'value'} } ) {\n				\
&print_param_value($value);\n			}\n		}\n		else {\n\
				&print_param_value($paramDetail->{'values'}->{\
'value'});\n		}\n	}\n	print_debug_message( 'print_\
param_details', 'End', 1 );\n}\n\n=head2 print_par\
am_value()\n\nPrint details of a tool parameter va\
lue.\n\n  &print_param_details($param_value);\n\nU\
sed by print_param_details() to handle both singlu\
ar and array values.\n\n=cut\n\nsub print_param_va\
lue {\n	my $value = shift;\n	print $value->{'value\
'};\n	if ( $value->{'defaultValue'} eq 'true' ) {\\
n		print \"\\t\", 'default';\n	}\n	print \"\\n\";\\
n	print \"\\t\", $value->{'label'}, \"\\n\";\n	if \
( defined( $value->{'properties'} ) ) {\n		foreach\
\n		  my $key ( sort( keys( %{ $value->{'propertie\
s'}{'property'} } ) ) )\n		{\n			if ( ref( $value-\
>{'properties'}{'property'}{$key} ) eq 'HASH'\n			\
	&& defined( $value->{'properties'}{'property'}{$k\
ey}{'value'} )\n			  )\n			{\n				print \"\\t\", $\
key, \"\\t\",\n				  $value->{'properties'}{'prope\
rty'}{$key}{'value'}, \"\\n\";\n			}\n			else {\n	\
			print \"\\t\", $value->{'properties'}{'property\
'}{'key'},\n				  \"\\t\", $value->{'properties'}{\
'property'}{'value'}, \"\\n\";\n				last;\n			}\n	\
	}\n	}\n}\n\n=head2 print_job_status()\n\nPrint st\
atus of a job.\n\n  &print_job_status($job_id);\n\\
n=cut\n\nsub print_job_status {\n	print_debug_mess\
age( 'print_job_status', 'Begin', 1 );\n	my $jobid\
 = shift;\n	print_debug_message( 'print_job_status\
', 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel > \
0 ) {\n		print STDERR 'Getting status for job ', $\
jobid, \"\\n\";\n	}\n	my $result = &rest_get_statu\
s($jobid);\n	print \"$result\\n\";\n	if ( $result \
eq 'FINISHED' && $outputLevel > 0 ) {\n		print STD\
ERR \"To get results: $scriptName --polljob --jobi\
d \" . $jobid\n		  . \"\\n\";\n	}\n	print_debug_me\
ssage( 'print_job_status', 'End', 1 );\n}\n\n=head\
2 print_result_types()\n\nPrint available result t\
ypes for a job.\n\n  &print_result_types($job_id);\
\n\n=cut\n\nsub print_result_types {\n	print_debug\
_message( 'result_types', 'Begin', 1 );\n	my $jobi\
d = shift;\n	print_debug_message( 'result_types', \
'jobid: ' . $jobid, 1 );\n	if ( $outputLevel > 0 )\
 {\n		print STDERR 'Getting result types for job '\
, $jobid, \"\\n\";\n	}\n	my $status = &rest_get_st\
atus($jobid);\n	if ( $status eq 'PENDING' || $stat\
us eq 'RUNNING' ) {\n		print STDERR 'Error: Job st\
atus is ', $status,\n		  '. To get result types th\
e job must be finished.', \"\\n\";\n	}\n	else {\n	\
	my (@resultTypes) = &rest_get_result_types($jobid\
);\n		if ( $outputLevel > 0 ) {\n			print STDOUT '\
Available result types:', \"\\n\";\n		}\n		foreach\
 my $resultType (@resultTypes) {\n			print STDOUT \
$resultType->{'identifier'}, \"\\n\";\n			if ( def\
ined( $resultType->{'label'} ) ) {\n				print STDO\
UT \"\\t\", $resultType->{'label'}, \"\\n\";\n			}\
\n			if ( defined( $resultType->{'description'} ) \
) {\n				print STDOUT \"\\t\", $resultType->{'desc\
ription'}, \"\\n\";\n			}\n			if ( defined( $resul\
tType->{'mediaType'} ) ) {\n				print STDOUT \"\\t\
\", $resultType->{'mediaType'}, \"\\n\";\n			}\n		\
	if ( defined( $resultType->{'fileSuffix'} ) ) {\n\
				print STDOUT \"\\t\", $resultType->{'fileSuffi\
x'}, \"\\n\";\n			}\n		}\n		if ( $status eq 'FINIS\
HED' && $outputLevel > 0 ) {\n			print STDERR \"\\\
n\", 'To get results:', \"\\n\",\n			  \"  $script\
Name --polljob --jobid \" . $params{'jobid'} . \"\\
\n\",\n			  \"  $scriptName --polljob --outformat \
<type> --jobid \"\n			  . $params{'jobid'} . \"\\n\
\";\n		}\n	}\n	print_debug_message( 'result_types'\
, 'End', 1 );\n}\n\n=head2 submit_job()\n\nSubmit \
a job to the service.\n\n  &submit_job($seq);\n\n=\
cut\n\nsub submit_job {\n	print_debug_message( 'su\
bmit_job', 'Begin', 1 );\n\n	# Set input sequence\\
n	$tool_params{'sequence'} = shift;\n\n	# Load par\
ameters\n	&load_params();\n\n	# Submit the job\n	m\
y $jobid = &rest_run( $params{'email'}, $params{'t\
itle'}, \\%tool_params );\n\n	# Simulate sync/asyn\
c mode\n	if ( defined( $params{'async'} ) ) {\n		p\
rint STDOUT $jobid, \"\\n\";\n		if ( $outputLevel \
> 0 ) {\n			print STDERR\n			  \"To check status: \
$scriptName --status --jobid $jobid\\n\";\n		}\n	}\
\n	else {\n		if ( $outputLevel > 0 ) {\n			print S\
TDERR \"JobId: $jobid\\n\";\n		}\n		sleep 1;\n		&g\
et_results($jobid);\n	}\n	print_debug_message( 'su\
bmit_job', 'End', 1 );\n}\n\n=head2 multi_submit_j\
ob()\n\nSubmit multiple jobs assuming input is a c\
ollection of fasta formatted sequences.\n\n  &mult\
i_submit_job();\n\n=cut\n\nsub multi_submit_job {\\
n	print_debug_message( 'multi_submit_job', 'Begin'\
, 1 );\n	my $jobIdForFilename = 1;\n	$jobIdForFile\
name = 0 if ( defined( $params{'outfile'} ) );\n	m\
y (@filename_list) = ();\n\n	# Query sequence\n	if\
 ( defined( $ARGV[0] ) ) {    # Bare option\n		if \
( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # File\n		\
	push( @filename_list, $ARGV[0] );\n		}\n		else {\\
n			warn 'Warning: Input file \"' . $ARGV[0] . '\"\
 does not exist'\n		}\n	}\n	if ( $params{'sequence\
'} ) {                   # Via --sequence\n		if ( \
-f $params{'sequence'} || $params{'sequence'} eq '\
-' ) {    # File\n			push( @filename_list, $params\
{'sequence'} );\n		}\n		else {\n			warn 'Warning: \
Input file \"' . $params{'sequence'} . '\" does no\
t exist'\n		}\n	}\n\n	$/ = '>';\n	foreach my $file\
name (@filename_list) {\n		my $INFILE;\n		if($file\
name eq '-') { # STDIN.\n			open( $INFILE, '<-' )\\
n			  or die 'Error: unable to STDIN (' . $! . ')'\
;\n		} else { # File.\n			open( $INFILE, '<', $fil\
ename )\n			  or die 'Error: unable to open file '\
 . $filename . ' (' . $! . ')';\n		}\n		while (<$I\
NFILE>) {\n			my $seq = $_;\n			$seq =~ s/>$//;\n	\
		if ( $seq =~ m/(\\S+)/ ) {\n				print STDERR \"S\
ubmitting job for: $1\\n\"\n				  if ( $outputLeve\
l > 0 );\n				$seq = '>' . $seq;\n				&print_debug\
_message( 'multi_submit_job', $seq, 11 );\n				&su\
bmit_job($seq);\n				$params{'outfile'} = undef if\
 ( $jobIdForFilename == 1 );\n			}\n		}\n		close $\
INFILE;\n	}\n	print_debug_message( 'multi_submit_j\
ob', 'End', 1 );\n}\n\n=head2 list_file_submit_job\
()\n\nSubmit multiple jobs using a file containing\
 a list of entry identifiers as \ninput.\n\n  &lis\
t_file_submit_job($list_filename)\n\n=cut\n\nsub l\
ist_file_submit_job {\n	print_debug_message( 'list\
_file_submit_job', 'Begin', 11 );\n	my $filename  \
       = shift;\n	my $jobIdForFilename = 1;\n	$job\
IdForFilename = 0 if ( defined( $params{'outfile'}\
 ) );\n\n	# Iterate over identifiers, submitting e\
ach job\n	my $LISTFILE;\n	if($filename eq '-') { #\
 STDIN.\n		open( $LISTFILE, '<-' )\n		  or die 'Er\
ror: unable to STDIN (' . $! . ')';\n	} else { # F\
ile.\n		open( $LISTFILE, '<', $filename )\n		  or \
die 'Error: unable to open file ' . $filename . ' \
(' . $! . ')';\n	}\n	while (<$LISTFILE>) {\n		my $\
line = $_;\n		chomp($line);\n		if ( $line ne '' ) \
{\n			&print_debug_message( 'list_file_submit_job'\
, 'line: ' . $line, 2 );\n			if ( $line =~ m/\\w:\\
\w/ ) {    # Check this is an identifier\n				prin\
t STDERR \"Submitting job for: $line\\n\"\n				  i\
f ( $outputLevel > 0 );\n				&submit_job($line);\n\
			}\n			else {\n				print STDERR\n\"Warning: line\
 \\\"$line\\\" is not recognised as an identifier\\
\n\";\n			}\n		}\n		$params{'outfile'} = undef if \
( $jobIdForFilename == 1 );\n	}\n	close $LISTFILE;\
\n	print_debug_message( 'list_file_submit_job', 'E\
nd', 11 );\n}\n\n=head2 load_data()\n\nLoad sequen\
ce data from file or option specified on the comma\
nd-line.\n\n  &load_data();\n\n=cut\n\nsub load_da\
ta {\n	print_debug_message( 'load_data', 'Begin', \
1 );\n	my $retSeq;\n\n	# Query sequence\n	if ( def\
ined( $ARGV[0] ) ) {    # Bare option\n		if ( -f $\
ARGV[0] || $ARGV[0] eq '-' ) {    # File\n			$retS\
eq = &read_file( $ARGV[0] );\n		}\n		else {       \
                              # DB:ID or sequence\\
n			$retSeq = $ARGV[0];\n		}\n	}\n	if ( $params{'s\
equence'} ) {                   # Via --sequence\n\
		if ( -f $params{'sequence'} || $params{'sequence\
'} eq '-' ) {    # File\n			$retSeq = &read_file( \
$params{'sequence'} );\n		}\n		else {    # DB:ID o\
r sequence\n			$retSeq = $params{'sequence'};\n		}\
\n	}\n	print_debug_message( 'load_data', 'End', 1 \
);\n	return $retSeq;\n}\n\n=head2 load_params()\n\\
nLoad job parameters from command-line options.\n\\
n  &load_params();\n\n=cut\n\nsub load_params {\n	\
print_debug_message( 'load_params', 'Begin', 1 );\\
n\n	# Database(s) to search\n	my (@dbList) = split\
 /[ ,]/, $params{'database'};\n	$tool_params{'data\
base'} = \\@dbList;\n\n	# Compatability options, o\
ld command-line.\n	if(!$tool_params{'viewfilter'} \
&& $params{'echofilter'}) {\n		$tool_params{'viewf\
ilter'} = 'true';\n	}\n	if(!$tool_params{'alignmen\
ts'} && $params{'numal'}) {\n		$tool_params{'align\
ments'} = $params{'numal'};\n	}\n	# TODO: set alig\
nment format option to get NCBI BLAST XML.\n	if($p\
arams{'appxml'}) {\n		$tool_params{'align'} = '';\\
n	}\n\n	print_debug_message( 'load_params', 'End',\
 1 );\n}\n\n=head2 client_poll()\n\nClient-side jo\
b polling.\n\n  &client_poll($job_id);\n\n=cut\n\n\
sub client_poll {\n	print_debug_message( 'client_p\
oll', 'Begin', 1 );\n	my $jobid  = shift;\n	my $st\
atus = 'PENDING';\n\n	my $errorCount = 0;\n	while \
($status eq 'RUNNING'\n		|| $status eq 'PENDING'\n\
		|| ( $status eq 'ERROR' && $errorCount < 2 ) )\n\
	{\n		$status = rest_get_status($jobid);\n		print \
STDERR \"$status\\n\" if ( $outputLevel > 0 );\n		\
if ( $status eq 'ERROR' ) {\n			$errorCount++;\n		\
}\n		elsif ( $errorCount > 0 ) {\n			$errorCount--\
;\n		}\n		if (   $status eq 'RUNNING'\n			|| $stat\
us eq 'PENDING'\n			|| $status eq 'ERROR' )\n		{\n\
\n			# Wait before polling again.\n			sleep $check\
Interval;\n		}\n	}\n	print_debug_message( 'client_\
poll', 'End', 1 );\n	return $status;\n}\n\n=head2 \
get_results()\n\nGet the results for a job identif\
ier.\n\n  &get_results($job_id);\n\n=cut\n\nsub ge\
t_results {\n	print_debug_message( 'get_results', \
'Begin', 1 );\n	my $jobid = shift;\n	print_debug_m\
essage( 'get_results', 'jobid: ' . $jobid, 1 );\n\\
n	# Verbose\n	if ( $outputLevel > 1 ) {\n		print '\
Getting results for job ', $jobid, \"\\n\";\n	}\n\\
n	# Check status, and wait if not finished\n	clien\
t_poll($jobid);\n\n	# Use JobId if output file nam\
e is not defined\n	unless ( defined( $params{'outf\
ile'} ) ) {\n		$params{'outfile'} = $jobid;\n	}\n\\
n	# Get list of data types\n	my (@resultTypes) = r\
est_get_result_types($jobid);\n\n	# Get the data a\
nd write it to a file\n	if ( defined( $params{'out\
format'} ) ) {    # Specified data type\n		my $sel\
ResultType;\n		foreach my $resultType (@resultType\
s) {\n			if ( $resultType->{'identifier'} eq $para\
ms{'outformat'} ) {\n				$selResultType = $resultT\
ype;\n			}\n		}\n		if ( defined($selResultType) ) \
{\n			my $result =\n			  rest_get_result( $jobid, \
$selResultType->{'identifier'} );\n			if ( $params\
{'outfile'} eq '-' ) {\n				write_file( $params{'o\
utfile'}, $result );\n			}\n			else {\n				write_f\
ile(\n					$params{'outfile'} . '.'\n					  . $sel\
ResultType->{'identifier'} . '.'\n					  . $selRes\
ultType->{'fileSuffix'},\n					$result\n				);\n		\
	}\n		}\n		else {\n			die 'Error: unknown result f\
ormat \"' . $params{'outformat'} . '\"';\n		}\n	}\\
n	else {    # Data types available\n		      # Writ\
e a file for each output type\n		for my $resultTyp\
e (@resultTypes) {\n			if ( $outputLevel > 1 ) {\n\
				print STDERR 'Getting ', $resultType->{'identi\
fier'}, \"\\n\";\n			}\n			my $result = rest_get_r\
esult( $jobid, $resultType->{'identifier'} );\n			\
if ( $params{'outfile'} eq '-' ) {\n				write_file\
( $params{'outfile'}, $result );\n			}\n			else {\\
n				write_file(\n					$params{'outfile'} . '.'\n	\
				  . $resultType->{'identifier'} . '.'\n					  \
. $resultType->{'fileSuffix'},\n					$result\n				\
);\n			}\n		}\n	}\n	print_debug_message( 'get_resu\
lts', 'End', 1 );\n}\n\n=head2 read_file()\n\nRead\
 a file into a scalar. The special filename '-' ca\
n be used to read from \nstandard input (STDIN).\n\
\n  my $data = &read_file($filename);\n\n=cut\n\ns\
ub read_file {\n	print_debug_message( 'read_file',\
 'Begin', 1 );\n	my $filename = shift;\n	print_deb\
ug_message( 'read_file', 'filename: ' . $filename,\
 2 );\n	my ( $content, $buffer );\n	if ( $filename\
 eq '-' ) {\n		while ( sysread( STDIN, $buffer, 10\
24 ) ) {\n			$content .= $buffer;\n		}\n	}\n	else \
{    # File\n		open( my $FILE, '<', $filename )\n	\
	  or die \"Error: unable to open input file $file\
name ($!)\";\n		while ( sysread( $FILE, $buffer, 1\
024 ) ) {\n			$content .= $buffer;\n		}\n		close($\
FILE);\n	}\n	print_debug_message( 'read_file', 'En\
d', 1 );\n	return $content;\n}\n\n=head2 write_fil\
e()\n\nWrite data to a file. The special filename \
'-' can be used to write to \nstandard output (STD\
OUT).\n\n  &write_file($filename, $data);\n\n=cut\\
n\nsub write_file {\n	print_debug_message( 'write_\
file', 'Begin', 1 );\n	my ( $filename, $data ) = @\
_;\n	print_debug_message( 'write_file', 'filename:\
 ' . $filename, 2 );\n	if ( $outputLevel > 0 ) {\n\
		print STDERR 'Creating result file: ' . $filenam\
e . \"\\n\";\n	}\n	if ( $filename eq '-' ) {\n		pr\
int STDOUT $data;\n	}\n	else {\n		open( my $FILE, \
'>', $filename )\n		  or die \"Error: unable to op\
en output file $filename ($!)\";\n		syswrite( $FIL\
E, $data );\n		close($FILE);\n	}\n	print_debug_mes\
sage( 'write_file', 'End', 1 );\n}\n\n=head2 usage\
()\n\nPrint program usage message.\n\n  &usage();\\
n\n=cut\n\nsub usage {\n	print STDERR <<EOF\nWU-BL\
AST\n========\n   \nRapid sequence database search\
 programs utilizing the BLAST algorithm\n    \n[Re\
quired]\n\n  -p, --program      : str  : BLAST pro\
gram to use, see --paramDetail program\n  -D, --da\
tabase     : str  : database(s) to search, space s\
eparated. See\n                              --par\
amDetail database\n      --stype        : str  : q\
uery sequence type, see --paramDetail stype\n  seq\
File            : file : query sequence (\"-\" for\
 STDIN, \\@filename for\n                         \
     identifier list file)\n\n[Optional]\n\n  -m, \
--matrix       : str  : scoring matrix, see --para\
mDetail matrix\n  -e, --exp          : real : 0<E<\
= 1000. Statistical significance threshold \n     \
                         for reporting database se\
quence matches.\n  -e, --viewfilter   :      : dis\
play the filtered query sequence\n  -f, --filter  \
     : str  : filter the query sequence for low co\
mplexity \n                              regions, \
see --paramDetail filter\n  -A, --align        : i\
nt  : pairwise alignment format, see --paramDetail\
 align\n  -s, --scores       : int  : number of sc\
ores to be reported\n  -b, --alignments   : int  :\
 number of alignments to report\n  -S, --sensitivi\
ty  : str  : sensitivity of the search, \n        \
                      see --paramDetail sensitivit\
y\n  -t, --sort	     : str  : sort order for hits,\
 see --paramDetail sort\n  -T, --stats        : st\
r  : statistical model, see --paramDetail stats\n \
 -d, --strand       : str  : DNA strand to search \
with,\n                              see --paramDe\
tail strand\n  -c, --topcombon    : str  : consist\
ent sets of HSPs\n      --multifasta   :      : tr\
eat input as a set of fasta formatted sequences\n\\
n[General]\n\n  -h, --help         :      : prints\
 this help text\n      --async        :      : for\
ces to make an asynchronous query\n      --email  \
      : str  : e-mail address\n      --title      \
  : str  : title for job\n      --status       :  \
    : get job status\n      --resultTypes  :      \
: get available result types for job\n      --poll\
job      :      : poll for the status of a job\n  \
    --jobid        : str  : jobid that was returne\
d when an asynchronous job \n                     \
         was submitted.\n      --outfile      : st\
r  : file name for results (default is jobid;\n   \
                           \"-\" for STDOUT)\n    \
  --outformat    : str  : result format to retriev\
e\n      --params       :      : list input parame\
ters\n      --paramDetail  : str  : display detail\
s for input parameter\n      --quiet        :     \
 : decrease output\n      --verbose      :      : \
increase output\n   \nSynchronous job:\n\n  The re\
sults/errors are returned as soon as the job is fi\
nished.\n  Usage: $scriptName --email <your\\@emai\
l> [options...] seqFile\n  Returns: results as an \
attachment\n\nAsynchronous job:\n\n  Use this if y\
ou want to retrieve the results at a later time. T\
he results \n  are stored for up to 24 hours. 	\n \
 Usage: $scriptName --async --email <your\\@email>\
 [options...] seqFile\n  Returns: jobid\n\n  Use t\
he jobid to query for the status of the job. If th\
e job is finished, \n  it also returns the results\
/errors.\n  Usage: $scriptName --polljob --jobid <\
jobId> [--outfile string]\n  Returns: string indic\
ating the status of the job and if applicable, res\
ults \n  as an attachment.\n\nFurther information:\
\n\n  http://www.ebi.ac.uk/Tools/webservices/servi\
ces/sss/wu_blast_rest\n  http://www.ebi.ac.uk/Tool\
s/webservices/tutorials/perl\n\nSupport/Feedback:\\
n\n  http://www.ebi.ac.uk/support/\nEOF\n}\n\n=hea\
d1 FEEDBACK/SUPPORT\n\nPlease contact us at L<http\
://www.ebi.ac.uk/support/> if you have any \nfeedb\
ack, suggestions or issues with the service or thi\
s client.\n\n=cut\n","\n\n\nmy $PROBTRESH = 0.3;# \
base pairs below this prob threshold will be ignor\
ed\nmy $WEIGHT = 100.0; # float!!\nmy $NUCALPH = \\
"ACGTUNRYMKSWHBVD\";\nuse vars qw($NUCALPH $WEIGHT\
);\n\nmy $myname = basename($0);\n\nuse strict;\nu\
se warnings;\n\nuse File::Basename;\nuse Getopt::L\
ong;\nuse File::Glob ':glob';\nuse File::Spec;\nus\
e File::Temp qw/ tempfile tempdir /;\n\n\n\n\nsub \
tcoffeelib_header($;$)\n{\n    my ($nseq, $fd) = @\
_;\n    if (! defined($fd)) {\n        $fd = *STDO\
UT;\n    }\n    printf $fd \"! TC_LIB_FORMAT_01\\n\
\";\n    printf $fd \"%d\\n\", $nseq;\n}\n\n\nsub \
tcoffeelib_header_addseq($$;$)\n{\n    my ($id, $s\
eq, $fd) = @_;\n    if (! defined($fd)) {\n       \
 $fd = *STDOUT;\n    }\n    printf $fd \"%s %d %s\\
\n\", $id, length($seq), $seq;\n}\n\n\nsub tcoffee\
lib_comment($;$)\n{\n    my ($comment, $fd) = @_;\\
n    if (! defined($fd)) {\n        $fd = *STDOUT;\
\n    }\n    printf $fd \"!\" . $comment . \"\\n\"\
;\n}\n\n\nsub tcoffeelib_struct($$$;$)\n{\n    my \
($nseq, $len, $bpm, $fd) = @_;\n\n    if (! define\
d($fd)) {\n        $fd = *STDOUT;\n    }\n\n    # \
output basepair indices with fixed weight\n    pri\
ntf $fd \"#%d %d\\n\", $nseq, $nseq;\n    # output\
 basepairs (only once) and with unit-offset\n    f\
or (my $i=0; $i<$len; $i++) {\n        for (my $j=\
$i+1; $j<$len; $j++) {\n            if (! defined(\
$bpm->[$i][$j])) {\n                print STDERR \\
"ERROR: \\$bpm->[$i][$j] undefined\\n\";\n        \
    }\n            if ($bpm->[$i][$j]>0) {\n      \
          print $fd $i+1;\n                print $\
fd \" \";\n                print $fd $j+1;\n      \
          print $fd \" \" . $bpm->[$i][$j] . \"\\n\
\";\n            }\n        }\n    }\n}\n\n\nsub t\
coffeelib_footer(;$)\n{\n    my ($fd) = @_;\n    i\
f (! defined($fd)) {\n        $fd = *STDOUT;\n    \
}\n    print $fd \"! SEQ_1_TO_N\\n\";\n}\n\n\n    \
\nsub plfold($$$)\n{    \n    my ($id, $seq, $prob\
tresh) = @_;\n    my (@struct);# return\n    my ($\
templ, $fhtmp, $fnametmp, $cmd, $ctr, $window_size\
);\n    our $ntemp++;\n    \n    $templ = $myname \
. \".pid-\" . $$ .$ntemp .\".XXXXXX\";\n    ($fhtm\
p, $fnametmp) = tempfile($templ, UNLINK => 1); \n \
   print $fhtmp \">$id\\n$seq\\n\";\n\n    # --- i\
nit basepair array\n    #\n    for (my $i=0; $i<le\
ngth($seq); $i++) {\n        for (my $j=$i+1; $j<l\
ength($seq); $j++) {\n            $struct[$i][$j]=\
0;\n        }\n    }\n\n\n    # --- call rnaplfold\
 and drop a readme\n    #\n    $window_size=(lengt\
h($seq)<70)?length($seq):70;\n    $cmd = \"RNAplfo\
ld -W $window_size < $fnametmp >/dev/null\";\n    \
system($cmd);\n    \n    if ($? != 0) {\n        p\
rintf STDERR \"ERROR: RNAplfold ($cmd) exited with\
 error status %d\\n\", $? >> 8;\n        return;\n\
    }\n    #unlink($fnametmp);\n    my $fps = spri\
ntf(\"%s_dp.ps\", $id); # check long name\n    \n \
   if (! -s $fps) {\n      {\n\n	$fps = sprintf(\"\
%s_dp.ps\", substr($id,0,12)); # check short name\\
n 	if (! -s $fps)\n	  {\n	    die(\"couldn't find \
expected file $fps\\n\");\n	    return;\n	  }\n   \
   }\n    }\n\n    \n    # --- read base pairs fro\
m created postscript\n    #\n    open(FH, $fps);\n\
    while (my $line = <FH>) {\n        my ($nti, $\
ntj, $prob);\n        chomp($line);        \n     \
   # line: bp bp sqrt-prob ubox\n        my @match\
 = ($line =~ m/^([0-9]+) +([0-9]+) +([0-9\\.]+) +u\
box$/);\n        if (scalar(@match)) {\n          \
  $nti=$1;\n            $ntj=$2;\n            $pro\
b=$3*$3;# prob stored as square root\n\n          \
  if ($prob>$probtresh) {\n                #printf\
 STDERR \"\\$struct[$nti][$ntj] sqrtprob=$3 prob=$\
prob > $probtresh\\n\";\n                $struct[$\
nti-1][$ntj-1] = $WEIGHT\n            }\n         \
   # store with zero-offset\n        }\n    }\n   \
 close(FH);\n\n    # remove or gzi postscript\n   \
 #\n    unlink($fps);\n    #\n    # or gzip\n    #\
$cmd = \"gzip -qf $fps\";\n    #system($cmd);\n   \
 #if ($? != 0) {\n    #    printf STDERR \"ERROR: \
gzip ($cmd) exited with error status %d\\n\", $? >\
> 8;\n    #}\n\n    return \\@struct;\n}\n\n\n\n\n\
\nsub rnaseqfmt($)\n{\n    my ($seq) = @_;\n    # \
remove gaps\n    $seq =~ s/-//g;\n    # uppercase \
RNA\n    $seq = uc($seq);\n    # T -> U\n    $seq \
=~ s/T/U/g;\n    # check for invalid charaters\n  \
  $_ = $seq;\n    s/[^$NUCALPH]//g;\n    return $_\
;\n}\n\n\n\n\nsub usage(;$)\n{    \n    my ($errms\
g) = @_;\n    if ($errmsg) {\n        print STDERR\
 \"ERROR: $errmsg\\n\";\n    }\n    print STDERR <\
< \"EOF\";\n$myname:\n Creates a T-Coffee RNA stru\
cture library from RNAplfold prediction.\n See FIX\
ME:citation\nUsage:\n $myname -in seq_file -out tc\
offee_lib\nEOF\n    exit(1);\n}\n\nsub read_fasta_\
seq \n  {\n    my $f=$_[0];\n    my %hseq;\n    my\
 (@seq, @com, @name);\n    my ($a, $s,$nseq);\n\n \
   open (F, $f);\n    while (<F>)\n      {\n	$s.=$\
_;\n      }\n    close (F);\n\n    \n    @name=($s\
=~/>(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.\
*.*\\n([^>]*)/g);\n    @com =($s=~/>(\\S*)(.*)\\n(\
[^>]*)/g);\n\n\n    $nseq=$#name+1;\n  \n    for (\
$a=0; $a<$nseq; $a++)\n      {\n	my $n=$name[$a];\\
n	my $s;\n	$hseq{$n}{name}=$n;\n	$s=$seq[$a];$s=~s\
/\\s//g;\n	\n	$hseq{$n}{seq}=$s;\n	$hseq{$n}{com}=\
$com[$a];\n      }\n    return %hseq;\n  }\n\n\n\n\
\n\n\n\nmy $fmsq = \"\";\nmy $flib = \"\";\nmy %OP\
TS;\nmy %seq;\nmy ($id, $nseq, $i);\nmy @nl;\n\nGe\
tOptions(\"in=s\" => \\$fmsq, \"out=s\" => \\$flib\
);\n\nif (! -s $fmsq) {\n    usage(\"empty or non-\
existant file \\\"$fmsq\\\"\")\n}\nif (length($fli\
b)==0) {\n    usage(\"empty out-filename\")\n}\n\n\
\n\n\n\n\n%seq=read_fasta_seq($fmsq);\n\n\n@nl=key\
s(%seq);\n\n$nseq=$#nl+1;\nopen FD_LIB, \">$flib\"\
 or die \"can't open $flib!\";\ntcoffeelib_header(\
$nseq, *FD_LIB);\nforeach $id (keys (%seq))\n  {\n\
    my ($seq, $fmtseq);\n    \n    $seq = $seq{$id\
}{seq};\n    \n    $fmtseq = rnaseqfmt($seq);# che\
ck here, formatting for folding important later\n \
   if (length($seq)!=length($fmtseq)) {\n        p\
rint STDERR \"ERROR: invalid sequence $id is not a\
n RNA sequence. read seq is: $seq\\n\";\n        e\
xit\n      }\n   \n    tcoffeelib_header_addseq($i\
d, uc($seq), *FD_LIB);\n  }\ntcoffeelib_comment(\"\
generated by $myname on \" . localtime(), *FD_LIB)\
;\n\n\n\n$i=0;\nforeach $id (keys (%seq))\n  {\n  \
  my ($cleanid, $seq, $bpm);\n    $seq=$seq{$id}{s\
eq};\n    $cleanid = $id;\n    $cleanid =~ s,[/ ],\
_,g;# needed for rnaplfold\n    $seq = rnaseqfmt($\
seq);\n    \n    $bpm = plfold($cleanid, rnaseqfmt\
($seq), $PROBTRESH);       \n    \n    tcoffeelib_\
struct($i+1, length($seq), $bpm, *FD_LIB);\n    $i\
++;\n}\n\n\ntcoffeelib_footer(*FD_LIB);\nclose FD_\
LIB;\nexit (0);\n\n","\n\n\n\n\n$cmd=join ' ', @AR\
GV;\nif ($cmd=~/-infile=(\\S+)/){ $seqfile=$1;}\ni\
f ($cmd=~/-outfile=(\\S+)/){ $libfile=$1;}\n\n\n\n\
%s=read_fasta_seq ($seqfile);\n\nopen (F, \">$libf\
ile\");\nforeach $name (keys (%s))\n  {\n    my $t\
clib=\"$name.RNAplfold_tclib\";\n    print (F \">$\
name _F_ $tclib\\n\");\n    seq2RNAplfold2tclib ($\
name, $s{$name}{seq}, $tclib);\n  }\nclose (F);\ne\
xit (EXIT_SUCCESS);\n\nsub seq2RNAplfold2tclib\n  \
{\n    my ($name, $seq, $tclib)=@_;\n    my ($tmp)\
;\n    $n++;\n    $tmp=\"tmp4seq2RNAplfold_tclib.$\
$.$n.pep\";\n    open (RF, \">$tmp\");\n    print \
(RF \">$name\\n$seq\\n\");\n    close (RF);\n    \\
n    system \"t_coffee -other_pg RNAplfold2tclib.p\
l -in=$tmp -out=$tclib\";\n    \n    unlink ($tmp)\
;\n    return $tclib;\n  }\n    \n    \nsub read_f\
asta_seq \n  {\n    my $f=@_[0];\n    my %hseq;\n \
   my (@seq, @com, @name);\n    my ($a, $s,$nseq);\
\n\n    open (F, $f);\n    while (<F>)\n      {\n	\
$s.=$_;\n      }\n    close (F);\n\n    \n    @nam\
e=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s\
=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\
\n([^>]*)/g);\n\n    \n    $nseq=$#name+1;\n    \n\
    for ($a=0; $a<$nseq; $a++)\n      {\n	my $n=$n\
ame[$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}=$s\
eq[$a];\n	$hseq{$n}{com}=$com[$a];\n      }\n    r\
eturn %hseq;\n  }\n","use Getopt::Long;\nuse File:\
:Path;\nuse Env;\nuse FileHandle;\nuse Cwd;\nuse S\
ys::Hostname;\nour $PIDCHILD;\nour $ERROR_DONE;\no\
ur @TMPFILE_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT\
_SUCCESS=0;\n\nour $REFDIR=getcwd;\nour $EXIT_SUCC\
ESS=0;\nour $EXIT_FAILURE=1;\n\nour $PROGRAM=\"tc_\
generic_method.pl\";\nour $CL=$PROGRAM;\n\nour $CL\
EAN_EXIT_STARTED;\nour $debug_lock=$ENV{\"DEBUG_LO\
CK\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\\
nif (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR\
=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$EN\
V{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (\
isshellpid(getppid())){lock4tc(getppid(), \"LLOCK\\
", \"LSET\", \"$$\\n\");}\n      \nour $print;\nmy\
 ($fmsq1, $fmsq2, $output, $outfile, $arch, $psv, \
$hmmtop_home, $trim, $cov, $sample, $mode, $gor_ho\
me, $gor_seq, $gor_obs);\n\nGetOptions(\"-in=s\" =\
> \\$fmsq1,\"-output=s\" =>\\$output ,\"-out=s\" =\
> \\$outfile, \"-arch=s\" => \\$arch,\"-psv=s\" =>\
 \\$psv, \"-hmmtop_home=s\", \\$hmmtop_home,\"-tri\
m=s\" =>\\$trim ,\"-print=s\" =>\\$print,\"-cov=s\\
" =>\\$cov , \"-sample=s\" =>\\$sample, \"-mode=s\\
" =>\\$mode, \"-gor_home=s\"=>\\$gor_home, \"-gor_\
seq=s\"=>\\$gor_seq,\"-gor_obs=s\"=>\\$gor_obs);\n\
\n\nif (!$mode){$mode = \"hmmtop\"}\nelsif ($mode \
eq \"hmmtop\"){;}\nelsif ($mode eq \"gor\"){;}\nel\
se {myexit(flush_error (\"-mode=$mode is unknown\"\
));}\n\n\nour $HOME=$ENV{\"HOME\"};\nour $MCOFFEE=\
($ENV{\"MCOFFEE_4_TCOFFEE\"})?$ENV{\"MCOFFEE_4_TCO\
FFEE\"}:\"$HOME/.t_coffee/mcoffee\";\n\nif ($mode \
eq \"hmmtop\")\n  {\n    \n    check_configuration\
 (\"hmmtop\");\n    if (-e $arch){$ENV{'HMMTOP_ARC\
H'}=$arch;}\n    elsif (-e $ENV{HMMTOP_ARCH}){$arc\
h=$ENV{HMMTOP_ARCH};}\n    elsif (-e \"$MCOFFEE/hm\
mtop.arch\"){$arch=$ENV{'HMMTOP_ARCH'}=\"$MCOFFEE/\
hmmtop.arch\";}\n    elsif (-e \"$hmmtop_home/hmmt\
op.arch\"){$arch=$ENV{'HMMTOP_ARCH'}=\"$hmmtop_hom\
e/hmmtop.arch\";}\n    else {myexit(flush_error ( \
\"Could not find ARCH file for hmmtop\"));}\n    \\
n    \n    if (-e $psv){$ENV{'HMMTOP_PSV'}=$psv;}\\
n    elsif (-e $ENV{HMMTOP_PSV}){$psv=$ENV{HMMTOP_\
PSV};}\n    elsif (-e \"$MCOFFEE/hmmtop.psv\"){$ps\
v=$ENV{'HMMTOP_PSV'}=\"$MCOFFEE/hmmtop.psv\";}\n  \
  elsif (-e \"$hmmtop_home/hmmtop.psv\"){$psv=$ENV\
{'HMMTOP_PSV'}=\"$hmmtop_home/hmmtop.psv\";}\n    \
else {myexit(flush_error ( \"Could not find PSV fi\
le for hmmtop\"));}\n\n  }\nelsif ($mode eq \"gor\\
")\n  {\n    our $GOR_SEQ;\n    our $GOR_OBS;\n   \
 \n    check_configuration (\"gorIV\");\n    if (-\
e $gor_seq){$GOR_SEQ=$gor_seq;}\n    elsif (-e $EN\
V{GOR_SEQ}){$GOR_SEQ=$ENV{GOR_SEQ};}\n    elsif (-\
e \"$MCOFFEE/New_KS.267.seq\"){$GOR_SEQ=\"$MCOFFEE\
/New_KS.267.seq\";}\n    elsif (-e \"$gor_home/New\
_KS.267.seq\"){$GOR_SEQ=\"$gor_home/New_KS.267.seq\
\";}\n    else {myexit(flush_error ( \"Could not f\
ind SEQ file for gor\"));}\n\n    if (-e $gor_obs)\
{$GOR_OBS=$gor_obs;}\n    elsif (-e $ENV{GOR_OBS})\
{$GOR_OBS=$ENV{GOR_OBS};}\n    elsif (-e \"$MCOFFE\
E/New_KS.267.obs\"){$GOR_OBS=\"$MCOFFEE/New_KS.267\
.obs\";}\n    elsif (-e \"$gor_home/New_KS.267.obs\
\"){$GOR_OBS=\"$gor_home/New_KS.267.obs\";}\n    e\
lse {myexit(flush_error ( \"Could not find OBS fil\
e for gor\"));}\n  }\n\n\nif ( ! -e $fmsq1){myexit\
(flush_error (\"Could Not Read Input file $fmsq1\"\
));}\n\n\nmy $fmsq2=vtmpnam();\nmy $fmsq3=vtmpnam(\
);\nmy $tmpfile=vtmpnam();\nmy $predfile=vtmpnam()\
;\n\nif ($trim){$trim_action=\" +trim _aln_%%$trim\
\\_K1 \";}\nif ($cov) {$cov_action= \" +sim_filter\
 _aln_c$cov \";}\n&safe_system(\"t_coffee -other_p\
g seq_reformat -in $fmsq1 -action +convert 'BOUJXZ\
-' $cov_action $trim_action -output fasta_aln -out\
 $fmsq2\");\nmy (%pred, %seq, %predA);\n\n\n%seq=r\
ead_fasta_seq($fmsq2);\n%seq=fasta2sample(\\%seq, \
$sample);\n\nif (1==2 &&$mode eq \"hmmtop\" && $ou\
tput eq \"cons\")\n  {\n    fasta2hmmtop_cons($out\
file,\\%seq);\n  }\nelse\n  {\n   \n    %pred=fast\
a2pred(\\%seq, $mode);\n    %predA=pred2aln (\\%pr\
ed, \\%seq);\n    \n    \n    if (!$output || $out\
put eq \"prediction\"){output_fasta_seq (\\%predA,\
 $outfile);}\n    elsif ($output eq \"color_html\"\
){pred2color (\\%pred,\\%seq, $outfile);}\n    els\
if ($output eq \"cons\"){pred2cons($outfile,\\%pre\
dA);}\n    else {flush_error (\"$output is an unkn\
own output mode\");}\n  }\n\nsub fasta2sample\n  {\
\n    my $SR=shift;\n    my $it=shift;\n    my %S=\
%$SR;\n    \n    my $seq=index2seq_name (\\%S, 1);\
\n    my $l=length($S{$seq}{seq});\n    my @sl=key\
s(%S);\n    my $nseq=$#sl+1;\n    my $index=$nseq;\
\n  \n    if (!$sample) {return %S;}\n    for (my \
$a=0; $a<$it; $a++)\n      {\n	my $newseq=\"\";\n	\
my $nname=\"$seq\\_sampled_$index\";\n	for (my $p=\
0; $p<$l; $p++)\n	  {\n	    my $i=int(rand($nseq))\
;\n	    \n	    my $name = $sl[$i];\n	    my $seq=$\
S{$name}{seq};\n	    my $r=substr ($seq, $p, 1);\n\
	    $newseq.=$r;\n	  }\n	$S{$nname}{name}=$nname;\
\n	$S{$nname}{seq}=$newseq;\n	$S{$nname}{com}=\"sa\
mpled\";\n	$S{$nname}{index}=++$index;\n      }\n \
   return %S;\n  }\n	      \nsub fasta2pred\n  {\n\
    my $s=shift;\n    my $mode=shift;\n\n    if ( \
$mode eq \"hmmtop\"){return fasta2hmmtop_pred($s);\
}\n    elsif ($mode eq \"gor\"){return fasta2gor_p\
red ($s);}\n  }\nsub fasta2hmmtop_cons\n  {\n    m\
y $outfile=shift;\n    my $SR=shift;\n    \n    my\
 $o = new FileHandle;\n    my $i = new FileHandle;\
\n    my $tmp_in =vtmpnam();\n    my $tmp_out=vtmp\
nam();\n    my %seq=%$SR;\n    my %pred;\n    my $\
N=keys(%seq);\n    \n    output_fasta_seq (\\%seq,\
$tmp_in, \"seq\");\n    `hmmtop -pi=mpred -if=$tmp\
_in -sf=FAS -pl 2>/dev/null >$tmp_out`;\n    open \
($o, \">$outfile\");\n    open ($i, \"$tmp_out\");\
\n    while (<$i>)\n      {\n	my $l=$_;\n	if (($l=\
~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\">$2 NSEQ:\
 $N\\n\";print $o \"$line\";}\n	elsif ( ($l=~/.*pr\
ed(.*)/))  {my $line=\"$1\\n\";print $o \"$line\";\
}\n      }\n    close ($o);\n    close ($i);\n    \
return read_fasta_seq($tmp);\n  }\nsub fasta2hmmto\
p_pred\n  {\n    my $SR=shift;\n    my $o = new Fi\
leHandle;\n    my $i = new FileHandle;\n    my $tm\
p    =vtmpnam();\n    my $tmp_in =vtmpnam();\n    \
my $tmp_out=vtmpnam();\n    my %seq=%$SR;\n    my \
%pred;\n    \n\n    output_fasta_seq (\\%seq,$tmp_\
in, \"seq\");\n\n    \n    `hmmtop -if=$tmp_in -sf\
=FAS -pl 2>/dev/null >$tmp_out`;\n    \n\n    \n  \
  \n    open ($o, \">$tmp\");\n    open ($i, \"$tm\
p_out\");\n    while (<$i>)\n      {\n	my $l=$_;\n\
	if (($l=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\"\
>$2\\n\";print $o \"$line\";}\n	elsif ( ($l=~/.*pr\
ed(.*)/))  {my $line=\"$1\\n\";print $o \"$line\";\
}\n      }\n    close ($o);\n    close ($i);\n    \
return read_fasta_seq($tmp);\n  }\n    \n	\n	\n	  \
  \n	\n	\n\n	\nsub fasta2gor_pred\n  {\n    my $SR\
=shift;\n    my $o = new FileHandle;\n    my $i = \
new FileHandle;\n    my $tmp    =vtmpnam();\n    m\
y $tmp_in =vtmpnam();\n    my $tmp_out=vtmpnam();\\
n    my %seq=%$SR;\n    my %pred;\n    \n\n    out\
put_fasta_seq (\\%seq,$tmp_in, \"seq\");\n    `gor\
IV -prd $tmp_in -seq $GOR_SEQ -obs $GOR_OBS >$tmp_\
out`;\n    open ($o, \">$tmp\");\n    open ($i, \"\
$tmp_out\");\n    while (<$i>)\n      {\n	my $l=$_\
;\n\n	\n	if ( $l=~/>/){print $o \"$l\";}\n	elsif (\
 $l=~/Predicted Sec. Struct./){$l=~s/Predicted Sec\
. Struct\\.//;print $o \"$l\";}\n      }\n    clos\
e ($o);\n    close ($i);\n    return read_fasta_se\
q($tmp);\n  }\n			\n			     \nsub index2seq_name\n\
  {\n    \n    my $SR=shift;\n    my $index=shift;\
\n    \n    \n    my %S=%$SR;\n    \n    foreach m\
y $s (%S)\n      {\n	if ( $S{$s}{index}==$index){r\
eturn $s;}\n      }\n    return \"\";\n  }\n\nsub \
pred2cons\n  {\n    my $outfile=shift;\n    my $pr\
edR=shift;\n    my $seq=shift;\n    my %P=%$predR;\
\n    my %C;\n    my ($s,@r,$nseq);\n    my $f= ne\
w FileHandle;\n\n    open ($f, \">$outfile\");\n\n\
    if (!$seq){$seq=index2seq_name(\\%P,1);}\n    \
foreach my $s (keys(%P))\n      {\n	$nseq++;\n	$st\
ring= $P{$s}{seq};\n	$string = uc $string;\n	my @r\
=split (//,$string);\n	for (my $a=0; $a<=$#r; $a++\
)\n	  {\n	    if (($r[$a]=~/[OHICE]/)){$C{$a}{$r[$\
a]}++;}\n	  }\n      }\n    @l=keys(%C);\n    \n  \
  \n    $s=$P{$seq}{seq};\n    print $f \">$seq pr\
ed based on $nseq\\n\";\n    @r=split (//,$s);\n  \
  \n    for (my $x=0; $x<=$#r; $x++)\n      {\n	if\
 ($r[$x] ne \"-\")\n	  {\n	    my $h=$C{$x}{H};\n	\
    my $i=$C{$x}{I};\n	    my $o=$C{$x}{O};\n	    \
my $c=$C{$x}{C};\n	    my $e=$C{$x}{E};\n	    my $\
l=$i+$o;\n	    \n	    if ($h>=$i && $h>=$o && $h>=\
$c && $h>=$e){$r[$x]='H';}\n	    elsif ($i>=$o && \
$i>=$c && $i>=$e){$r[$x]='I';}\n	    elsif ($o>=$c\
 && $o>=$e){$r[$x]='O';}\n	    elsif ($c>=$e){$r[$\
x]='C';}\n	    else {$r[$x]='E';}\n	  }\n      }\n\
    $j=join ('', @r);\n    print $f \"$j\\n\";\n  \
  close ($f);\n    return $j;\n  }\n\nsub pred2aln\
\n  {\n    my $PR=shift;\n    my $AR=shift;\n    \\
n    my $f=new FileHandle;\n    my %P=%$PR;\n    m\
y %A=%$AR;\n    my %PA;\n    my $tmp=vtmpnam();\n \
   my $f= new FileHandle;\n    \n    open ($f, \">\
$tmp\");\n    foreach my $s (sort{$A{$a}{index}<=>\
$A{$b}{index}}(keys (%A)))\n      {\n	my (@list, $\
seq, @plist, @pseq, $L, $PL, $c, $w);\n	my $seq;\n\
	my $seq=$A{$s}{seq};\n	my $pred=$P{$s}{seq};\n	$s\
eq=pred2alnS($P{$s}{seq},$A{$s}{seq});\n	print $f \
\">$s\\n$seq\\n\";\n      }\n    close ($f);\n    \
return read_fasta_seq ($tmp);\n  }\nsub pred2alnS\\
n  {\n    my $pred=shift;\n    my $aln= shift;\n  \
  my ($j,$a,$b);\n    my @P=split (//, $pred);\n  \
  my @A=split (//, $aln);\n    for ($a=$b=0;$a<=$#\
A; $a++)\n      {\n	if ($A[$a] ne \"-\"){$A[$a]=$P\
[$b++];}\n      }\n    if ($b!= ($#P+1)){add_warni\
ng (\"Could not thread sequence: $b $#P\");}\n    \
\n    $j= join ('', @A);\n    return $j;\n  }\nsub\
 pred2color\n  {\n    my $predP=shift;\n    my $al\
nP=shift;\n    my $out=shift;\n    my $F=new FileH\
andle;\n    my $struc=vtmpnam();\n    my $aln=vtmp\
nam();\n    \n\n    output_fasta_seq ($alnP, $aln)\
;\n    my %p=%$predP;\n    \n    open ($F, \">$str\
uc\");\n    \n    \n    foreach my $s (keys(%p))\n\
      {\n	\n	print $F \">$s\\n\";\n	my $s=uc($p{$s\
}{seq});\n	\n	$s=~s/[Oo]/0/g;\n	$s=~s/[Ee]/0/g;\n	\
\n	$s=~s/[Ii]/5/g;\n	$s=~s/[Cc]/5/g;\n	\n	$s=~s/[H\
h]/9/g;\n	\n	print $F \"$s\\n\";\n      }\n    clo\
se ($F);\n    \n    \n    \n    safe_system ( \"t_\
coffee -other_pg seq_reformat -in $aln -struc_in $\
struc -struc_in_f number_fasta -output color_html \
-out $out\");\n    return;\n  }\n	  \n    \nsub di\
splay_fasta_seq\n  {\n    my $SR=shift;\n    my %S\
=%$SR;\n    \n    foreach my $s (sort{$S{$a}{index\
}<=>$S{$b}{index}}(keys (%S)))\n      {\n	print ST\
DERR \">$s\\n$S{$s}{seq}\\n\";\n      }\n    close\
 ($f);\n  }\nsub output_fasta_seq\n  {\n    my $SR\
=shift;\n    my $outfile=shift;\n    my $mode =shi\
ft;\n    my $f= new FileHandle;\n    my %S=%$SR;\n\
    \n    \n    open ($f, \">$outfile\");\n    for\
each my $s (sort{$S{$a}{index}<=>$S{$b}{index}}(ke\
ys (%S)))\n      {\n	my $seq=$S{$s}{seq};\n	if ( $\
mode eq \"seq\"){$seq=~s/\\-//g;}\n	print $f \">$s\
\\n$seq\\n\";\n      }\n    close ($f);\n  }\n    \
  \nsub read_fasta_seq \n  {\n    my $f=$_[0];\n  \
  my %hseq;\n    my (@seq, @com, @name);\n    my (\
$a, $s,$nseq);\n    my $index;\n    open (F, $f);\\
n    while (<F>)\n      {\n	$s.=$_;\n      }\n    \
close (F);\n\n    \n    @name=($s=~/>(\\S*).*\\n[^\
>]*/g);\n    \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\
\n    @com =($s=~/>.*(.*)\\n([^>]*)/g);\n\n\n    $\
nseq=$#name+1;\n    \n  \n    for ($a=0; $a<$nseq;\
 $a++)\n      {\n	my $n=$name[$a];\n	my $s;\n	$hse\
q{$n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\n	$hse\
q{$n}{index}=++$index;\n	$hseq{$n}{seq}=$s;\n	$hse\
q{$n}{com}=$com[$a];\n      }\n    return %hseq;\n\
  }\n\n\nsub file2head\n      {\n	my $file = shift\
;\n	my $size = shift;\n	my $f= new FileHandle;\n	m\
y $line;\n	open ($f,$file);\n	read ($f,$line, $siz\
e);\n	close ($f);\n	return $line;\n      }\nsub fi\
le2tail\n      {\n	my $file = shift;\n	my $size = \
shift;\n	my $f= new FileHandle;\n	my $line;\n	\n	o\
pen ($f,$file);\n	seek ($f,$size*-1, 2);\n	read ($\
f,$line, $size);\n	close ($f);\n	return $line;\n  \
    }\n\n\nsub vtmpnam\n      {\n	my $r=rand(10000\
0);\n	my $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\\
n	    $f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $\
f);\n	return $f;\n      }\n\nsub myexit\n  {\n    \
my $code=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){r\
eturn;}\n    else {$CLEAN_EXIT_STARTED=1;}\n    ##\
# ONLY BARE EXIT\n    exit ($code);\n  }\nsub set_\
error_lock\n    {\n      my $name = shift;\n      \
my $pid=$$;\n\n      \n      &lock4tc ($$,\"LERROR\
\", \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\
\n      return;\n    }\nsub set_lock\n  {\n    my \
$pid=shift;\n    my $msg= shift;\n    my $p=getppi\
d();\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p\
$msg\\n\");\n  }\nsub unset_lock\n   {\n     \n   \
 my $pid=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LR\
ELEASE\",\"\");\n  }\nsub shift_lock\n  {\n    my \
$from=shift;\n    my $to=shift;\n    my $from_type\
=shift;\n    my $to_type=shift;\n    my $action=sh\
ift;\n    my $msg;\n    \n    if (!&lock4tc($from,\
 $from_type, \"LCHECK\", \"\")){return 0;}\n    $m\
sg=&lock4tc ($from, $from_type, \"LREAD\", \"\");\\
n    &lock4tc ($from, $from_type,\"LRELEASE\", $ms\
g);\n    &lock4tc ($to, $to_type, $action, $msg);\\
n    return;\n  }\nsub isshellpid\n  {\n    my $p=\
shift;\n    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\\
")){return 0;}\n    else\n      {\n	my $c=lock4tc(\
$p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){r\
eturn 1;}\n      }\n    return 0;\n  }\nsub isroot\
pid\n  {\n    if(lock4tc (getppid(), \"LLOCK\", \"\
LCHECK\")){return 0;}\n    else {return 1;}\n  }\n\
sub lock4tc\n	{\n	  my ($pid,$type,$action,$value)\
=@_;\n	  my $fname;\n	  my $host=hostname;\n	  \n	\
  if ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.\
$host.lock4tcoffee\";}\n	  elsif ( $type eq \"LERR\
OR\"){ $fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\
\";}\n	  elsif ( $type eq \"LWARNING\"){ $fname=\"\
$LOCKDIR/.$pid.$host.warning4tcoffee\";}\n	  \n	  \
if ($debug_lock)\n	    {\n	      print STDERR \"\\\
n\\t---lock4tc(tcg): $action => $fname =>$value (R\
D: $LOCKDIR)\\n\";\n	    }\n\n	  if    ($action eq\
 \"LCHECK\") {return -e $fname;}\n	  elsif ($actio\
n eq \"LREAD\"){return file2string($fname);}\n	  e\
lsif ($action eq \"LSET\") {return string2file ($v\
alue, $fname, \">>\");}\n	  elsif ($action eq \"LR\
ESET\") {return string2file ($value, $fname, \">\"\
);}\n	  elsif ($action eq \"LRELEASE\") \n	    {\n\
	      if ( $debug_lock)\n		{\n		  my $g=new FileH\
andle;\n		  open ($g, \">>$fname\");\n		  print $g\
 \"\\nDestroyed by $$\\n\";\n		  close ($g);\n		  \
safe_system (\"mv $fname $fname.old\");\n		}\n	   \
   else\n		{\n		  unlink ($fname);\n		}\n	    }\n	\
  return \"\";\n	}\n	\nsub file2string\n	{\n	  my \
$file=@_[0];\n	  my $f=new FileHandle;\n	  my $r;\\
n	  open ($f, \"$file\");\n	  while (<$f>){$r.=$_;\
}\n	  close ($f);\n	  return $r;\n	}\nsub string2f\
ile \n    {\n    my ($s,$file,$mode)=@_;\n    my $\
f=new FileHandle;\n    \n    open ($f, \"$mode$fil\
e\");\n    print $f  \"$s\";\n    close ($f);\n  }\
\n\nBEGIN\n    {\n      srand;\n    \n      $SIG{'\
SIGUP'}='signal_cleanup';\n      $SIG{'SIGINT'}='s\
ignal_cleanup';\n      $SIG{'SIGQUIT'}='signal_cle\
anup';\n      $SIG{'SIGILL'}='signal_cleanup';\n  \
    $SIG{'SIGTRAP'}='signal_cleanup';\n      $SIG{\
'SIGABRT'}='signal_cleanup';\n      $SIG{'SIGEMT'}\
='signal_cleanup';\n      $SIG{'SIGFPE'}='signal_c\
leanup';\n      \n      $SIG{'SIGKILL'}='signal_cl\
eanup';\n      $SIG{'SIGPIPE'}='signal_cleanup';\n\
      $SIG{'SIGSTOP'}='signal_cleanup';\n      $SI\
G{'SIGTTIN'}='signal_cleanup';\n      $SIG{'SIGXFS\
Z'}='signal_cleanup';\n      $SIG{'SIGINFO'}='sign\
al_cleanup';\n      \n      $SIG{'SIGBUS'}='signal\
_cleanup';\n      $SIG{'SIGALRM'}='signal_cleanup'\
;\n      $SIG{'SIGTSTP'}='signal_cleanup';\n      \
$SIG{'SIGTTOU'}='signal_cleanup';\n      $SIG{'SIG\
VTALRM'}='signal_cleanup';\n      $SIG{'SIGUSR1'}=\
'signal_cleanup';\n\n\n      $SIG{'SIGSEGV'}='sign\
al_cleanup';\n      $SIG{'SIGTERM'}='signal_cleanu\
p';\n      $SIG{'SIGCONT'}='signal_cleanup';\n    \
  $SIG{'SIGIO'}='signal_cleanup';\n      $SIG{'SIG\
PROF'}='signal_cleanup';\n      $SIG{'SIGUSR2'}='s\
ignal_cleanup';\n\n      $SIG{'SIGSYS'}='signal_cl\
eanup';\n      $SIG{'SIGURG'}='signal_cleanup';\n \
     $SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG\
{'SIGXCPU'}='signal_cleanup';\n      $SIG{'SIGWINC\
H'}='signal_cleanup';\n      \n      $SIG{'INT'}='\
signal_cleanup';\n      $SIG{'TERM'}='signal_clean\
up';\n      $SIG{'KILL'}='signal_cleanup';\n      \
$SIG{'QUIT'}='signal_cleanup';\n      \n      our \
$debug_lock=$ENV{\"DEBUG_LOCK\"};\n      \n      \\
n      \n      \n      foreach my $a (@ARGV){$CL.=\
\" $a\";}\n      if ( $debug_lock ){print STDERR \\
"\\n\\n\\n********** START PG: $PROGRAM **********\
***\\n\";}\n      if ( $debug_lock ){print STDERR \
\"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ **\
***********\\n\";}\n      if ( $debug_lock ){print\
 STDERR \"\\n --- $$ -- $CL\\n\";}\n      \n	     \
\n      \n      \n    }\nsub flush_error\n  {\n   \
 my $msg=shift;\n    return add_error ($EXIT_FAILU\
RE,$$, $$,getppid(), $msg, $CL);\n  }\nsub add_err\
or \n  {\n    my $code=shift;\n    my $rpid=shift;\
\n    my $pid=shift;\n    my $ppid=shift;\n    my \
$type=shift;\n    my $com=shift;\n    \n    $ERROR\
_DONE=1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\
\"$pid -- ERROR: $type\\n\");\n    lock4tc ($$, \"\
LERROR\",\"LSET\", \"$pid -- COM: $com\\n\");\n   \
 lock4tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK\
: $ppid -> $pid\\n\");\n   \n    return $code;\n  \
}\nsub add_warning \n  {\n    my $rpid=shift;\n   \
 my $pid =shift;\n    my $command=shift;\n    my $\
msg=\"$$ -- WARNING: $command\\n\";\n    print STD\
ERR \"$msg\";\n    lock4tc ($$, \"LWARNING\", \"LS\
ET\", $msg);\n  }\n\nsub signal_cleanup\n  {\n    \
print dtderr \"\\n**** $$ (tcg) was killed\\n\";\n\
    &cleanup;\n    exit ($EXIT_FAILURE);\n  }\nsub\
 clean_dir\n  {\n    my $dir=@_[0];\n    if ( !-d \
$dir){return ;}\n    elsif (!($dir=~/tmp/)){return\
 ;}#safety check 1\n    elsif (($dir=~/\\*/)){retu\
rn ;}#safety check 2\n    else\n      {\n	`rm -rf \
$dir`;\n      }\n    return;\n  }\nsub cleanup\n  \
{\n    #print stderr \"\\n----tc: $$ Kills $PIDCHI\
LD\\n\";\n    #kill (SIGTERM,$PIDCHILD);\n    my $\
p=getppid();\n    $CLEAN_EXIT_STARTED=1;\n    \n  \
  \n    \n    if (&lock4tc($$,\"LERROR\", \"LCHECK\
\", \"\"))\n      {\n	my $ppid=getppid();\n	if (!$\
ERROR_DONE) \n	  {\n	    &lock4tc($$,\"LERROR\", \\
"LSET\", \"$$ -- STACK: $p -> $$\\n\");\n	    &loc\
k4tc($$,\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\\
");\n	  }\n      }\n    my $warning=&lock4tc($$, \\
"LWARNING\", \"LREAD\", \"\");\n    my $error=&loc\
k4tc($$,  \"LERROR\", \"LREAD\", \"\");\n    #rele\
ase error and warning lock if root\n    \n    if (\
isrootpid() && ($warning || $error) )\n      {\n	\\
n	print STDERR \"**************** Summary ********\
*****\\n$error\\n$warning\\n\";\n\n	&lock4tc($$,\"\
LERROR\",\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNI\
NG\",\"RELEASE\",\"\");\n      } \n    \n    \n   \
 foreach my $f (@TMPFILE_LIST)\n      {\n	if (-e $\
f){unlink ($f);} \n      }\n    foreach my $d (@TM\
PDIR_LIST)\n      {\n	clean_dir ($d);\n      }\n  \
  #No More Lock Release\n    #&lock4tc($$,\"LLOCK\\
",\"LRELEASE\",\"\"); #release lock \n\n    if ( $\
debug_lock ){print STDERR \"\\n\\n\\n********** EN\
D PG: $PROGRAM ($$) *************\\n\";}\n    if (\
 $debug_lock ){print STDERR \"\\n\\n\\n**********(\
tcg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n  \
}\nEND \n  {\n    \n    &cleanup();\n  }\n   \n\ns\
ub safe_system \n{\n  my $com=shift;\n  my $ntry=s\
hift;\n  my $ctry=shift;\n  my $pid;\n  my $status\
;\n  my $ppid=getppid();\n  if ($com eq \"\"){retu\
rn 1;}\n  \n  \n\n  if (($pid = fork ()) < 0){retu\
rn (-1);}\n  if ($pid == 0)\n    {\n      set_lock\
($$, \" -SHELL- $com (tcg)\");\n      exec ($com);\
\n    }\n  else\n    {\n      lock4tc ($$, \"LLOCK\
\", \"LSET\", \"$pid\\n\");#update parent\n      $\
PIDCHILD=$pid;\n    }\n  if ($debug_lock){printf S\
TDERR \"\\n\\t .... safe_system (fasta_seq2hmm)  p\
: $$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,W\
TERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\\
"LWARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAILU\
RE || lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\
\n    {\n      if ($ntry && $ctry <$ntry)\n	{\n	  \
add_warning ($$,$$,\"$com failed [retry: $ctry]\")\
;\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\\
");\n	  return safe_system ($com, $ntry, ++$ctry);\
\n	}\n      elsif ($ntry == -1)\n	{\n	  if (!shift\
_lock ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\\
"))\n	    {\n	      add_warning ($$,$$,\"$com fail\
ed\");\n	    }\n	  else\n	    {\n	      lock4tc ($\
pid, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  \
return $?;}\n      else\n	{\n	  if (!shift_lock ($\
pid,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\\
n	      myexit(add_error ($EXIT_FAILURE,$$,$pid,ge\
tppid(), \"UNSPECIFIED system\", $com));\n	    }\n\
	}\n    }\n  return $?;\n}\n\nsub check_configurat\
ion \n    {\n      my @l=@_;\n      my $v;\n      \
foreach my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMA\
IL\")\n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n\
		add_warning($$,$$,\"Could Not Use EMAIL\");\n		m\
yexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"E\
MAIL\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p \
eq \"INTERNET\")\n	    {\n	      if ( !&check_inte\
rnet_connection())\n		{\n		  myexit(add_error ($EX\
IT_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\
\n		}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\n\
	      if (!&pg_is_installed (\"wget\") && !&pg_is\
_installed (\"curl\"))\n		{\n		  myexit(add_error \
($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:\
wget\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is\
_installed ($p)))\n	    {\n	      myexit(add_error\
 ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED\
:$p\",\"$CL\"));\n	    }\n	}\n      return 1;\n   \
 }\nsub pg_is_installed\n  {\n    my @ml=@_;\n    \
my $r, $p, $m;\n    my $supported=0;\n    \n    my\
 $p=shift (@ml);\n    if ($p=~/::/)\n      {\n	if \
(safe_system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){\
return 1;}\n	else {return 0;}\n      }\n    else\n\
      {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \\
"\"){return 0;}\n	else {return 1;}\n      }\n  }\n\
\n\n\nsub check_internet_connection\n  {\n    my $\
internet;\n    my $tmp;\n    &check_configuration \
( \"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \n\
    if     (&pg_is_installed    (\"wget\")){`wget \
www.google.com -O$tmp >/dev/null 2>/dev/null`;}\n \
   elsif  (&pg_is_installed    (\"curl\")){`curl w\
ww.google.com -o$tmp >/dev/null 2>/dev/null`;}\n  \
  \n    if ( !-e $tmp || -s $tmp < 10){$internet=0\
;}\n    else {$internet=1;}\n    if (-e $tmp){unli\
nk $tmp;}\n\n    return $internet;\n  }\nsub check\
_pg_is_installed\n  {\n    my @ml=@_;\n    my $r=&\
pg_is_installed (@ml);\n    if (!$r && $p=~/::/)\n\
      {\n	print STDERR \"\\nYou Must Install the p\
erl package $p on your system.\\nRUN:\\n\\tsudo pe\
rl -MCPAN -e 'install $pg'\\n\";\n      }\n    els\
if (!$r)\n      {\n	myexit(flush_error(\"\\nProgra\
m $p Supported but Not Installed on your system\")\
);\n      }\n    else\n      {\n	return 1;\n      \
}\n  }\n\n\n\n","\n\n\n\n\nmy $FMODEL =\"\"; \nmy \
$TMPDIR = \"/tmp\";\n\n\n\n\nmy $NUCALPH = \"ACGTU\
NRYMKSWHBVD\";\nmy $PRIMNUCALPH = \"ACGTUN\";\nuse\
 vars qw($NUCALPH $PRIMNUCALPH $TMPDIR);\n\n\nmy $\
errmsg;\nuse vars qw($errmsg);\n\n\n\nuse Getopt::\
Long;\nuse Cwd;\nuse File::Basename;\nuse File::Te\
mp qw/ tempfile tempdir /;\nuse File::Copy;\nuse F\
ile::Path;\n\n\n\nsub usage(;$)\n{\n    my ($errms\
g) = @_;\n    my $myname = basename($0);\n\n    if\
 ($errmsg) {\n        print STDERR \"ERROR: $errms\
g\\n\";\n    }\n\n    print STDERR << \"EOF\";\n  \
  \n$myname: align two sequences by means of consa\
n\\'s sfold\nUsage:\n $myname -i file -o file -d p\
ath\nOptions:\n -i|--in : pairwise input sequence \
file\n -o|--out: output alignment\n -d|--directory\
 containing data\n\nEOF\n}\n\nsub read_stk_aln \n \
 {\n    my $f=$_[0];\n    my ($seq, $id);\n    \n \
   my %hseq;\n\n    open (STK, \"$f\");\n    while\
 (<STK>)\n      {\n	if ( /^#/ || /^\\/\\// || /^\\\
s*$/){;}\n	else\n	  {\n	    ($id,$seq)=/(\\S+)\\s+\
(\\S+)/;\n	    $hseq{$id}{'seq'}.=$seq;\n	  }\n   \
   }\n    close (STK);\n    return %hseq;\n  }\nsu\
b read_fasta_seq \n  {\n    my $f=$_[0];\n    my %\
hseq;\n    my (@seq, @com, @name);\n    my ($a, $s\
,$nseq);\n\n    open (F, $f);\n    while (<F>)\n  \
    {\n	$s.=$_;\n      }\n    close (F);\n\n    \n\
    @name=($s=~/>(.*).*\\n[^>]*/g);\n    \n    @se\
q =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>.*(\
.*)\\n([^>]*)/g);\n\n    \n    $nseq=$#name+1;\n  \
  \n    for ($a=0; $a<$nseq; $a++)\n      {\n	my $\
n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq\
}=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n      }\n \
   return %hseq;\n  }\n\n\n\nsub sfold_parseoutput\
($$)\n{\n    my ($frawout, $foutfa) = @_;\n    my \
%haln;\n    my ($fstk, $cmd, $id);\n    open FOUTF\
A, \">$foutfa\";\n    \n    $fstk = $frawout . \".\
stk\";\n    \n    # first line of raw out contains\
 info\n    # remaining stuff is stockholm formatte\
d\n    $cmd = \"sed -e '1d' $frawout\";\n    syste\
m(\"$cmd > $fstk\");\n    if ($? != 0) {\n        \
$errmsg = \"command failed with exit status $?.\";\
\n        $errmsg .=  \"Command was \\\"$cmd\\\"\"\
;\n        return -1;\n    }\n\n    # this gives a\
n error message. just ignore it...\n    %haln=read\
_stk_aln ( $fstk);\n    foreach $i (keys (%haln))\\
n      {\n	my $s;\n	$s=$haln{$i}{'seq'};\n	$s =~ s\
/\\./-/g;\n	print FOUTFA \">$i\\n$s\\n\";\n      }\
\n    close FOUTFA;\n    return 0;\n}\n\n\n\n\nsub\
 sfold_wrapper($$$$)\n{\n    \n    my ($fs1, $fs2,\
 $fmodel, $foutfa) = @_;\n    \n\n    my ($cmd, $f\
rawout, $ferrlog, $freadme, $ftimelog, $fstk);\n\n\
    # add  basename($fmsqin) (unknown here!)\n    \
$frawout = \"sfold.log\";\n    $ferrlog = \"sfold.\
err\";\n    $ftimelog = \"sfold.time\";\n    $frea\
dme =  \"sfold.README\";\n    $fstk = \"sfold.stk\\
";\n    \n    # prepare execution...\n    #\n    #\
 ./tmp is essential for dswpalign\n    # otherwise\
 you'll get a segfault\n    mkdir \"./tmp\";\n    \
\n    $cmd = \"sfold -m $fmodel $fs1 $fs2\";\n    \
open(FREADME,\">$freadme\");\n    print FREADME \"\
$cmd\\n\"; \n    close(FREADME);\n\n    # and go\n\
    #\n    system(\"/usr/bin/time -p -o $ftimelog \
$cmd >$frawout 2>$ferrlog\");\n    if ($? != 0) {\\
n        $errmsg = \"command failed with exit stat\
us $?\";\n        $errmsg .= \"command was \\\"$cm\
d\\\". See \" . getcwd . \"\\n\";\n        return \
-1;\n    }\n\n    return sfold_parseoutput($frawou\
t, $foutfa);\n}\n\n\n\n\n\n\n\nmy ($help, $fmsqin,\
 $fmsaout);\nGetOptions(\"help\"  => \\$help,\n   \
        \"in=s\" => \\$fmsqin,\n           \"out=s\
\" => \\$fmsaout,\n	   \"data=s\" => \\$ref_dir);\\
n\n\n\nif ($help) {\n    usage();\n    exit(0);\n}\
\nif (! defined($fmsqin)) {\n    usage('missing in\
put filename');\n    exit(1);\n}\nif (! defined($f\
msaout)) {\n    usage('missing output filename');\\
n    exit(1);\n\n}\nif (scalar(@ARGV)) {\n    usag\
e('Unknown remaining args');\n    exit(1);\n}\n\n$\
FMODEL = \"$ref_dir/mix80.mod\";\nif (! -e \"$FMOD\
EL\") {\n    die(\"couldn't find sfold grammar mod\
el file. Expected $FMODEL\\n\");\n}\n\n\nmy %hseq=\
read_fasta_seq ($fmsqin);\nmy $id;\n\nforeach $id \
(keys(%hseq))\n  {\n    push(@seq_array, $hseq{$id\
});\n  }\n\nif ( scalar(@seq_array) != 2 ) {\n    \
die(\"Need *exactly* two sequences as input (pairw\
ise alignment!).\")\n}\n\n\n\nmy ($sec, $min, $hou\
r, $mday, $mon, $year, $wday, $yday, $isdst) = loc\
altime(time);\nmy $datei = sprintf(\"%4d-%02d-%02d\
\", $year+1900, $mon+1, $mday);\nmy $templ = basen\
ame($0) . \".\" . $datei . \".pid-\" . $$ . \".XXX\
XXX\";\nmy $wd = tempdir ( $templ, DIR => $TMPDIR)\
;\n\ncopy($fmsqin, \"$wd/\" . basename($fmsqin) . \
\".org\"); # for reproduction\ncopy($FMODEL, \"$wd\
\");\nmy $fmodel = basename($FMODEL);\nmy $orgwd =\
 getcwd;\nchdir $wd;\n\n\n\nmy @sepseqfiles;\nfore\
ach $id (keys(%hseq)) {\n    my ($seq, $orgseq, $f\
name, $sout);\n    $seq=$hseq{$id}{'seq'};\n    \n\
    $fname = basename($fmsqin) . \"_$id.fa\";\n   \
 # replace funnies in file/id name (e.g. \"/\" \" \
\" etc)\n    $fname =~ s,[/ ],_,g;\n    open (PF, \
\">$fname\");\n    print (PF \">$id\\n$seq\\n\");\\
n    close (PF);\n\n    push(@sepseqfiles, $fname)\
;\n}\n\nmy ($f1, $f2, $fout);\n$f1 = $sepseqfiles[\
0];\n$f2 = $sepseqfiles[1];\n$fout = $wd . basenam\
e($fmsqin) . \".out.fa\";\nif (sfold_wrapper($f1, \
$f2, $fmodel, \"$fout\") != 0) {\n    printf STDER\
R \"ERROR: See logs in $wd\\n\";\n    exit(1);\n} \
else {\n    chdir $orgwd;\n    copy($fout, $fmsaou\
t);\n    rmtree($wd);\n   exit(0);\n}\n","\nuse En\
v qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\\
n\n\n$tmp=clean_cr ($ARGV[0]);\nopen (F, $tmp);\n\\
nwhile ( <F>)\n  {\n    my $l=$_;\n    if ( $l=~/^\
# STOCKHOLM/){$stockholm=1;}\n    elsif ( $stockho\
lm && $l=~/^#/)\n      {\n	$l=~/^#(\\S+)\\s+(\\S+)\
\\s+(\\S*)/g;\n	$l=\"_stockholmhasch_$1\\_stockhol\
mspace_$2 $3\\n\";\n      }\n    $file.=$l;\n  }\n\
close (F);\nunlink($tmp);\n$file1=$file;\n\n$file=\
~s/\\#/_hash_symbol_/g;\n$file=~s/\\@/_arobase_sym\
bol_/g;\n\n\n$file=~s/\\n[\\.:*\\s]+\\n/\\n\\n/g;\\
n\n$file=~s/\\n[ \\t\\r\\f]+(\\b)/\\n\\1/g;\n\n\n$\
file=~s/(\\n\\S+)(\\s+)(\\S)/\\1_blank_\\3/g;\n\n$\
file=~s/[ ]//g;\n$file=~s/_blank_/ /g;\n\n\n\n$fil\
e =~s/\\n\\s*\\n/#/g;\n\n$file.=\"#\";\n$file =~s/\
\\n/@/g;\n\n\n\n\n@blocks=split /\\#/, $file;\nshi\
ft (@blocks);\n@s=split /\\@/, $blocks[0];\n$nseq=\
$#s+1;\n\n\n\n$file=join '@', @blocks;\n@lines=spl\
it /\\@/,$file;\n\n$c=0;\n\nforeach $l (@lines)\n \
 {\n    if (!($l=~/\\S/)){next;}\n    elsif ($stoc\
kholm && ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){next;\
}#get read of STOCHOLM Terminator\n   \n    $l=~/(\
\\S+)\\s+(\\S*)/g;\n    $n=$1; $s=$2;\n    \n    $\
seq[$c].=$s;\n    $name[$c]=$n;\n    $c++;\n    \n\
    if ( $c==$nseq){$c=0;}\n    \n  } \n\nif ( $c!\
=0)\n      {\n	print STDERR \"ERROR: $ARGV[0] is N\
OT an MSA in Clustalw format: make sure there is n\
o blank line within a block [ERROR]\\n\";\n	exit (\
EXIT_FAILURE);\n      }\n\nfor ($a=0; $a< $nseq; $\
a++)\n  {\n    $name[$a]=cleanstring ($name[$a]);\\
n    $seq[$a]=cleanstring ($seq[$a]);\n    $seq[$a\
]=breakstring($seq[$a], 60);\n    \n    $line=\">$\
name[$a]\\n$seq[$a]\\n\";\n    \n    print \"$line\
\";\n  }\nexit (EXIT_SUCCESS);\n\nsub cleanstring\\
n  {\n    my $s=@_[0];\n    $s=~s/_hash_symbol_/\\\
#/g;\n    $s=~s/_arobase_symbol_/\\@/g;\n    $s=~s\
/[ \\t]//g;\n    return $s;\n  }\nsub breakstring\\
n  {\n    my $s=@_[0];\n    my $size=@_[1];\n    m\
y @list;\n    my $n,$ns, $symbol;\n    \n    @list\
=split //,$s;\n    $n=0;$ns=\"\";\n    foreach $sy\
mbol (@list)\n      {\n	if ( $n==$size)\n	  {\n	  \
  $ns.=\"\\n\";\n	    $n=0;\n	  }\n	$ns.=$symbol;\\
n	$n++;\n      }\n    return $ns;\n    }\n\nsub cl\
ean_cr\n  {\n    my $f=@_[0];\n    my $file;\n    \
\n    $tmp=\"f$.$$\";\n    \n    \n    open (IN, $\
f);\n    open (OUT, \">$tmp\");\n    \n    while (\
 <IN>)\n      {\n	$file=$_;\n	$file=~s/\\r\\n/\\n/\
g;\n	$file=~s/\\n\\r/\\n/g;\n	$file=~s/\\r\\r/\\n/\
g;\n	$file=~s/\\r/\\n/g;\n	print OUT \"$file\";\n \
     }\n    \n    close (IN);\n    close (OUT);\n \
   return $tmp;\n  }\n","use strict;\nuse FileHand\
le;\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env\
 qw(USER);\n\nmy $format=file2format ($ARGV[0]);\n\
\nif    ($format eq \"clustalw\"){clustalw2fasta($\
ARGV[0]);}\nelsif ($format eq \"fasta\")   {fasta2\
fasta($ARGV[0]);}\nelsif ($format eq \"msf\")   {m\
sf2fasta($ARGV[0]);}\nelsif ($format eq \"phylip\"\
)   {phylip2fasta($ARGV[0]);}\nelsif ($format eq \\
"nameseq\") {display_file ($ARGV[0]);}\n \nexit (0\
);\n\nsub file2format\n  {\n    my $f=shift;\n    \
\n    my $l=file2n_lines($f,2);\n    \n    if ( $l\
=~/^CLUSTAL/){return \"clustalw\";}\n    elsif ($l\
=~/^SAGA/){return \"clustalw\";}\n    elsif ($l=~/\
^>/){return \"fasta\";}\n    elsif ($l=~/^PileUp/)\
{return \"msf\";}\n    elsif ($l=~/\\s+\\d+\\s+\\d\
+\\s/){return \"phylip\";}\n    elsif ($l=~/\\#NAM\
ESEQ_01/){return \"nameseq\";}\n    else \n      {\
\n	print STDERR \"ERROR: $f FILE is NOT a supporte\
d format [ERROR]\\n\";\n	system (\"cp $f /Users/cn\
otredame/cedric1.txt\");\n	exit (1);\n      }\n  }\
\nsub display_file\n    {\n       my $file=shift;\\
n       my $F= new FileHandle;\n       open ($F, $\
file);\n       while (<$F>){print \"$_\";}\n      \
 close ($F);\n     }\nsub phylip2fasta\n    {\n   \
   my $file=shift;\n      my $F= new FileHandle;\n\
      my ($seq, $name,$seq);\n      my $query_star\
t=-1;\n      my $query_end=-1;\n      my $in_aln=0\
;\n      my %list;\n      my ($first,$seq,$name, $\
cn, $nseq, $l,%len);\n      \n      open ($F, $fil\
e);\n      <$F>;\n      $l=$_;\n      $l=~/\\s*(\\\
d+)\\s*(\\d+)/;\n      $first=1;\n      $cn=0;\n  \
    while (<$F>)\n	{\n	  my $l=$_;\n	  if (!($l=~/\
\\S/))\n	    {\n	      $cn=0;\n	      $first=0;\n	\
    }\n	  elsif ($first==1)\n	    {\n	      $l=~/\\
\s*(\\S+)(.*)/;\n	      my $name=$1;\n	      my $s\
eq=$2;\n	      chomp ($seq);\n	      $seq=~s/\\s//\
g;\n	      $list{$cn}{'name'}=$name;\n	      $list\
{$cn}{'seq'}.=$seq;\n	      $cn++;\n	      $nseq++\
;\n	    }\n	  else\n	    {\n	      chomp ($l);\n	 \
     $l=~s/\\s//g;\n	      $list{$cn}{'seq'}.=$l;\\
n	      $cn++;\n	    }\n	}\n      close ($F);\n   \
   \n      for (my $a=0; $a<$nseq; $a++)\n	{\n	  p\
rint \">$list{$a}{'name'}\\n$list{$a}{'seq'}\\n\";\
\n	}\n    }\n      \nsub msf2fasta\n    {\n      m\
y $file=shift;\n      my $F= new FileHandle;\n    \
  my ($seq, $name,$seq);\n      my $query_start=-1\
;\n      my $query_end=-1;\n      my $in_aln=0;\n \
     my %list;\n      my ($seq,$name, $n, $nseq, $\
l,%len);\n      \n      open ($F, $file);\n      w\
hile (<$F>)\n	{\n	  if ( /\\/\\//){$in_aln=1;}\n	 \
 elsif ( $in_aln && /(\\S+)\\s+(.*)/)\n	    {\n	  \
    $name=$1;\n	      $seq=$2;\n	      $seq=~s/\\s\
//g;\n	      $seq=~s/\\~/\\-/g;\n	      $seq=~s/\\\
./\\-/g;\n	      if ( $list{$n}{'name'} && $list{$\
n}{'name'} ne $name)\n		{\n		  print \"$list{$n}{'\
name'} Vs $name\";\n		  \n		  exit (1);\n		}\n	   \
   else\n		{\n		  $list{$n}{'name'}= $name;\n		}\n\
	      \n	      $list{$n}{'seq'}=$list{$n}{'seq'}.\
$seq;\n	      \n	      $nseq=++$n;\n	      \n	    \
}\n	  else\n	    {$n=0;}\n	}\n      close ($F);\n \
     \n      for (my $a=0; $a<$nseq; $a++)\n	{\n	 \
 my $nl=length ($list{$a}{'name'});\n	  my $sl=len\
gth ($list{$a}{'seq'});\n	  print \">$list{$a}{'na\
me'}\\n$list{$a}{'seq'}\\n\";\n	}\n    }\n    \nsu\
b fasta2fasta\n    {\n      my $file=shift;\n     \
 my $F= new FileHandle;\n      my ($seq, $name,$n,\
$l,%len);\n      my $started=0;\n      open ($F, $\
file);\n      while (<$F>)\n	{\n	  if ( /^>(\\S+)/\
){$n++;$seq=\"\";$name=$1;}\n	  else\n	    {\n	   \
   $l=$_;\n	      chomp ($l);\n	      \n	      $se\
q.=$l;\n	      $len{$name}=length($seq);\n	    }\n\
	}\n      close ($F);\n      \n      open ($F, $fi\
le);\n      while (<$F>)\n	{\n	  my $l=$_;\n	  $l=\
~s/\\r[\\n]*/\\n/gm;\n	  if ( ($l=~/^>(\\S+)(.*)\\\
n/))\n	    {\n	      my $name=$1;\n	      my $comm\
ent=$2;\n	      my $nl=length ($name);\n	      my \
$sl=$len{$name};\n	      if ($comment)\n		{\n		  $\
comment=~s/^\\s+//g;\n		  my $cl=length ($comment)\
;\n		}\n	      if (!$started){$started=1;print \">\
$name\\n\";}\n	      else {print \"\\n>$name\\n\"}\
\n	    }\n	  else\n	    {\n	      $l=$_;\n	      c\
homp ($l);\n	      $l=~s/\\W//g;\n	      print \"$\
l\";\n	    }\n	}\n      print \"\\n\";\n      clos\
e ($F);\n    }\nsub clustalw2fasta\n  {\n    my $f\
name=shift;\n    my ($file1, $file);\n    my (@blo\
cks, @lines,@s, $n,$nseq, $c);\n    my (@name, @se\
q);\n    my $F= new FileHandle;\n    my $stockholm\
;\n   \n    \n    open ($F, $fname);\n    \n    wh\
ile ( <$F>)\n      {\n	my $l=$_;\n	$l=clean_cr($l)\
;\n	if ( $l=~/^# STOCKHOLM/){$stockholm=1;}\n	elsi\
f ( $stockholm && $l=~/^#/)\n	  {\n	    $l=~/^#(\\\
S+)\\s+(\\S+)\\s+(\\S*)/g;\n	    $l=\"_stockholmha\
sch_$1\\_stockholmspace_$2 $3\\n\";\n	  }\n	$file.\
=$l;\n      }\n    close ($F);\n        \n    #Pro\
tect # and @\n    $file=~s/\\#/_hash_symbol_/g;\n \
   $file=~s/\\@/_arobase_symbol_/g;\n    \n    \n \
   #Remove annotation\n    $file=~s/\\n[\\.:*\\s]+\
\\n/\\n\\n/g;\n    \n    #Remove White spaces befo\
re the sequence name\n    $file=~s/\\n[ \\t\\r\\f]\
+(\\b)/\\n\\1/g;\n    \n    \n    #Remove Internal\
 Blanks\n    $file=~s/(\\n\\S+)(\\s+)(\\S)/\\1_bla\
nk_\\3/g;\n    \n    $file=~s/[ ]//g;\n    $file=~\
s/_blank_/ /g;\n    \n    \n    #Identify Double B\
lank lines\n    \n    $file =~s/\\n\\s*\\n/#/g;\n \
   \n    $file.=\"#\";\n    $file =~s/\\n/@/g;\n  \
  \n    \n    \n    \n    #count nseq\n    @blocks\
=split /\\#/, $file;\n    shift (@blocks);\n    @s\
=split /\\@/, $blocks[0];\n    $nseq=$#s+1;\n    \\
n    #Merge all the sequences and split every Nseq\
\n    \n    \n    $file=join '@', @blocks;\n    @l\
ines=split /\\@/,$file;\n    \n    $c=0;\n    \n  \
  foreach my $l (@lines)\n      {\n	my ($n, $s);\n\
	\n	if (!($l=~/\\S/)){next;}\n	elsif ($stockholm &\
& ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){next;}#get r\
ead of STOCHOLM Terminator\n	\n	$l=~/(\\S+)\\s+(\\\
S*)/g;\n	$n=$1; $s=$2;\n	\n	$seq[$c].=$s;\n	$name[\
$c]=$n;\n	$c++;\n	\n	if ( $c==$nseq){$c=0;}\n	\n  \
    } \n    \n    if ( $c!=0)\n      {\n	print STD\
ERR \"ERROR: $fname is NOT an MSA in Clustalw form\
at: make sure there is no blank line within a bloc\
k [ERROR]\\n\";\n	exit (1);\n      }\n    \n    \n\
    for (my $a=0; $a< $nseq; $a++)\n      {\n	$nam\
e[$a]=cleanstring ($name[$a]);\n	$seq[$a]=cleanstr\
ing ($seq[$a]);\n	print \">$name[$a]\\n$seq[$a]\\n\
\";\n      }\n  }\nsub cleanstring\n    {\n      m\
y $s=@_[0];\n      $s=~s/_hash_symbol_/\\#/g;\n   \
   $s=~s/_arobase_symbol_/\\@/g;\n      $s=~s/[ \\\
t]//g;\n      return $s;\n    }\n\nsub clean_cr\n \
 {\n    my $f=shift;\n    $f=~s/\\r\\n/\\n/g;\n   \
 $f=~s/\\n\\r/\\n/g;\n    $f=~s/\\r\\r/\\n/g;\n   \
 $f=~s/\\r/\\n/g;\n    return $f;\n  }\n\nsub file\
2n_lines\n    {\n      my $file=shift;\n      my $\
nl=shift;\n      my $ret;\n      my $F=new FileHan\
dle;\n      my $n=0;\n      open ($F, $file);\n\n \
     while (<$F>)\n	{\n	  $ret.=$_;\n	  $n++;\n	  \
\n	  if ($n>=$n){close ($F); return $ret;}\n	}\n  \
    close ($F);\n      return $ret;\n    }\n","use\
 strict;\nuse FileHandle;\nuse Env qw(HOST);\nuse \
Env qw(HOME);\nuse Env qw(USER);\nmy %name;\nmy $n\
seq;\nmy $fasta;\nif ($ARGV[2] eq \"-fasta\"){$fas\
ta=1;}\nmy $F= new FileHandle;\n\nopen ($F, $ARGV[\
1]);\nwhile(<$F>)\n  {\n    my $l=$_;\n    if ($l=\
~/^#/){;}\n    elsif (($l=~/\\d+\\s+\\d+\\s+(\\S+)\
\\s+(\\S+)/))\n      {\n	my $n=$1;\n	$name{$1}++;\\
n      }\n  }\nclose ($F);\n\nopen ($F, $ARGV[0]);\
\nwhile(<$F>)\n  {\n    my $l=$_;\n    if ($l=~/^#\
/){;}\n    elsif ($l=~/\\d+\\s+\\d+\\s+(\\S+)\\s+(\
\\S+)/)\n      {\n	my $n=$1;\n	$name{$n}++;\n	if (\
$name{$n}==2){$nseq++;}\n      }\n  }\nclose ($F);\
\n\nif (!$fasta && $nseq>0)\n  {\n    print \"#NAM\
ESEQ_01\\n\";\n    print \"# $nseq\\n\";\n  }\nope\
n ($F, $ARGV[0]);\nwhile(<$F>)\n  {\n    my $l=$_;\
\n    if ($l=~/^#/){;}\n    elsif ($l=~/.\\d+\\s+\\
\d+\\s+(\\S+)\\s+(\\S+)/)\n      {\n	my $n=$1;\n	m\
y $s=$2;\n	if ($name{$n}==2)\n	  {\n	    if ($fast\
a)\n	      {\n		print \">$n\\n$s\\n\";\n	      }\n\
	    else\n	      {\n		print \"$l\";\n	      }\n	 \
 }\n      }\n  }\nclose ($F);\nexit (0);\n\n\n","u\
se Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(US\
ER);\n\n\n$query_start=-1;\n$query_end=-1;\n\nwhil\
e (<>)\n  {\n    if ( /\\/\\//){$in_aln=1;}\n    e\
lsif ( $in_aln && /(\\S+)\\s+(.*)/)\n      {\n\n\n\
	$name=$1;\n	\n\n	$seq=$2;\n	$seq=~s/\\s//g;\n    \
    $seq=~s/\\~/\\-/g;\n	$seq=~s/\\./\\-/g;\n	if (\
 $list{$n}{'name'} && $list{$n}{'name'} ne $name)\\
n	  {\n	    print \"$list{$n}{'name'} Vs $name\";\\
n	    \n	    exit (EXIT_FAILURE);\n	  }\n	else\n	 \
 {\n	    $list{$n}{'name'}= $name;\n	  }\n\n	$list\
{$n}{'seq'}=$list{$n}{'seq'}.$seq;\n	\n	$nseq=++$n\
;\n	\n      }\n    else\n      {$n=0;}\n  }\n\n\nf\
or ($a=0; $a<$nseq; $a++)\n  {\n    print \">$list\
{$a}{'name'}\\n$list{$a}{'seq'}\\n\";\n  }\n      \
\n","$run_anyway=2;\nmy $msaf=\"msa.in.tmp.$$\";\n\
my $msaoutf=\"msa.out.tmp.$$\";\nmy $err=\"msa.out\
.err.$$\";\nopen  (F, $ARGV[0]);\nopen  (OUT, \">$\
msaf\");\n$nseq=0;\nwhile (<F>)\n  {\n    $l=$_;\n\
    if ( $l=~/^>(\\S+)/)\n      {\n	$s=$seqname{$n\
seq++}=$1;\n	print OUT \"$l\";\n	\n      }\n    el\
se \n      {\n	$l=uc($l);\n	print OUT \"$l\";\n   \
   }\n  }\n\nclose (F);\nclose(OUT);\n\nsystem (\"\
msa $msaf > $msaoutf 2>$err\");\nopen (F, \"$msaou\
tf\");\n$read=0;\n$cn=0;\nwhile (<F>)\n  {\n    $l\
=$_;\n    if ($read)\n      {\n	if ($l=~/End gaps \
not penalized/){$read=0;}\n	elsif (!($l=~/\\S/))\n\
	  {\n	    $cn=0;\n	  }\n	else\n	  {\n	    \n	    \
chomp ($l);\n	    $seqal{$cn++}.=$l;\n	    $tot++;\
\n	  }\n      }\n    elsif ($l=~/Optimal Multiple \
Alignment/)\n      {\n	$read=1;\n      }\n  }\nclo\
se (F);\n\nif ($tot<1 && $run_anyway==1)\n  {\n   \
 print STDERR \"\\nWarning: MSA returned a NULL fi\
le -- Use T-Coffee instead\\n\";\n    open (F,$err\
);\n    while (<F>){print \"$_\";}\n      \n    sy\
stem (\"t_coffee -seq $msaf -outfile $ARGV[1]  -qu\
iet\");\n  }\nelsif ($tot<1 && $run_anyway==2)\n  \
{\n    \n    \n    $nseq/=2;\n    $nseq=int ($nseq\
);\n    if ($nseq<2){$nseq=2;}\n    print \"RUN MS\
A with NSeq=$nseq\\n\";\n    #print (\"t_coffee -d\
pa -dpa_nseq $nseq -seq $ARGV[0] -dpa_tree codnd -\
outfile $ARGV[1] -dpa_method msa_msa\");\n    syst\
em (\"t_coffee -dpa -dpa_nseq $nseq -seq $ARGV[0] \
-dpa_tree codnd -outfile $ARGV[1] -dpa_method msa_\
msa>/dev/null\");\n\n  }\nelsif ($tot<1)\n  {\n   \
 exit (EXIT_FAILURE);\n  }\nelse\n  {\n    open (O\
UT, \">$ARGV[1]\");\n    for ($a=0; $a<$nseq;$a++)\
\n      {\n	print OUT \">$seqname{$a}\\n$seqal{$a}\
\\n\";\n      }\n    close (OUT);\n  }\n\n\n\nunli\
nk ($msaf);\nunlink ($msaoutf);\nunlink ($err);\n"\
,"use strict;\nuse Cwd;\nuse File::Basename;\nmy $\
test=0;\n\nmy $tmpdir=\"/tmp/tco/aligners/upp/\";\\
nmymkdir ($tmpdir);\n\n\n\nif ($ARGV[0] eq \"one\"\
)\n  {\n    seq2msa ($ARGV[1], $ARGV[2]);\n  }\nel\
sif ($ARGV[0] eq \"all\")\n  {\n    listseq2listms\
a ($ARGV[1]);\n  }\n\nsub listseq2listmsa\n  {\n  \
  my $list=shift;\n    my $cdir = getcwd;\n    my \
$dir=random_string(10);\n    $dir=\"$tmpdir/$dir/\\
";\n    my %h;\n    my $n;\n    mkdir  ($dir);\n\n\
    open (F, $list);\n    while (<F>)\n      {\n  \
      my $l=$_;\n\n        chomp($l);\n        my \
@f=split (/\\s+/, $l);\n	if ( -e $f[0])\n         \
 {\n            $h{$n}{in}=$f[0];\n            ($h\
{$n}{name},$h{$n}{path})=fileparse ($f[0]);\n     \
       $h{$n}{NFin}= \"$dir/$h{$n}{name}.seq\";\n	\
    \n            $h{$n}{NFout}=\"$dir/$h{$n}{name\
}.aln\";\n\n            $h{$n}{out}=$f[1];\n\n    \
        fasta2fastaupp ($h{$n}{in}, $h{$n}{NFin});\
\n            $n++;\n          }\n      }\n    clo\
se (F);\n    chdir ($dir);\n    \n    if (!$test)\\
n      {\n	system (\"fbname=\\$(basename `ls *.seq\
` .seq); \\\n             run_upp.py -s \\${fbname\
}.seq -m amino --cpu 1 -d outdir -o \\${fbname}.al\
n; \\\n             mv outdir/\\${fbname}.aln_alig\
nment.fasta \\${fbname}.aln;\");\n      }\n    \n \
   foreach my $n (keys (%h))\n      {\n	if ($test)\
\n	  {\n	    system (\"cp $h{$n}{NFin} $h{$n}{NFou\
t}\");\n	    print \"$h{$n}{NFin} $h{$n}{NFout} $h\
{$n}{out}\\n\";\n	  }\n        fastaupp2fasta ($h{\
$n}{NFout},$h{$n}{out});\n      }\n    chdir ($cdi\
r);\n  }\n\nsub seq2msa\n    {\n      my ($in, $ou\
t)=@_;\n      my $cdir=getcwd;\n      \n      \n  \
    if (!($in=~/\\//)){$in=$cdir.\"/\".$in;}\n    \
  if (!($out=~/\\//)){$out=$cdir.\"/\".$out;}\n   \
   \n      my $file=random_string(10);\n      $fil\
e=\"$tmpdir/$file\";\n      open (F, \">$file\");\\
n      print F \"$in $out\\n\";\n      close (F);\\
n      \n      return listseq2listmsa ($file);\n  \
  }\n	\nsub fasta2fastaupp\n  {\n    my ($in, $out\
)=@_;\n    my ($name, $seq, $n);\n    \n    if (!-\
e $in){return;}\n    \n    open (IN, \"$in\");\n  \
  open (OUT, \">$out\");\n    local $/ = \"\\n>\";\
  # read by FASTA record\n    \n    while (<IN>)\n\
      {\n	my $l=$_;\n	$l=~s/>//g;\n	$l=\">\".$l;\n\
	\n	$l=~/^>(.*)/;\n	$name=$1;\n	\n	$l=~s/^>*.+\\n/\
/;\n	$l=~s/\\n//g;\n	$seq=$l;\n	\n	$seq=~s/u/x/g;\\
n	$seq=~s/U/X/g;\n	print OUT \">$name\\n$seq\\n\";\
\n	$n++;\n      }\n    if ($n==2)\n      {\n	print\
 OUT \">fake_seq4upp\\n$seq\\n\";\n      }\n    cl\
ose (IN);\n    close (OUT);\n    local $/=\"\\n\";\
\n  }\n\nsub fastaupp2fasta\n  {\n    my ($in, $ou\
t)=@_;\n    my ($name, $seq, $n);\n    \n    if (!\
-e $in){return;}\n    \n    open (IN, \"$in\");\n \
   open (OUT, \">$out\");\n    local $/ = \"\\n>\"\
;  # read by FASTA record\n    \n    while (<IN>)\\
n      {\n	my $l=$_;\n	$l=~s/>//g;\n	$l=\">\".$l;\\
n	\n	$l=~/^>(.*)/;\n	$name=$1;\n	\n	$l=~s/^>*.+\\n\
//;\n	$l=~s/\\n//g;\n	$seq=$l;\n	\n	$seq=~s/x/u/g;\
\n	$seq=~s/X/U/g;\n	\n	if (!($name=~/fake_seq4upp/\
))\n	  {\n	    print OUT \">$name\\n$seq\\n\";\n	 \
 }\n      }\n    close (IN);\n    close (OUT);\n  \
  local $/=\"\\n\";\n  }\n\nsub random_string\n   \
 {\n      my $len=shift;\n      my @chars = (\"A\"\
..\"Z\", \"a\"..\"z\");\n      my $string;\n      \
$string .= $chars[rand @chars] for 1..$len;\n     \
 return $string;\n    }\n\nsub mymkdir\n      {\n	\
my $d=shift;\n	my $cd='/';\n	\n	foreach my $e (spl\
it (/\\//, $d))\n	  {\n	    $cd.=\"$e/\";\n	    if\
 ( !-d $cd){mkdir ($cd);}\n	  }\n	return;\n      }\
\n      \n			  \n      \n","use strict;\nuse Cwd;\\
nuse File::Basename;\n\n\nmy $tmpdir=\"/tmp/tco/al\
igners/clustalo/\";\nmymkdir ($tmpdir);\n\n\n\nif \
($ARGV[0] eq \"one\")\n  {\n    seq2msa ($ARGV[1],\
 $ARGV[2]);\n  }\nelsif ($ARGV[0] eq \"all\")\n  {\
\n    listseq2listmsa ($ARGV[1]);\n  }\n\n\n\nsub \
listseq2listmsa\n  {\n    my $list=shift;\n    my \
$cdir = getcwd;\n    my $dir=random_string(10);\n \
   $dir=\"$tmpdir/$dir/\";\n    my %h;\n    my $n;\
\n    mkdir  ($dir);\n    \n    open (F, $list);\n\
    while (<F>)\n      {\n	my $l=$_;\n\n	chomp($l)\
;\n	my @f=split (/\\s+/, $l);\n	#print \"$l: 0:$f[\
0], 1:$f[1]\\n\";\n	if ( -e $f[0])\n	  {\n	    $h{\
$n}{in}=$f[0];\n	    ($h{$n}{name},$h{$n}{path})=f\
ileparse ($f[0]);\n	    $h{$n}{NFin}= \"$dir/$h{$n\
}{name}.seq4nf\";\n	    $h{$n}{NFout}=\"$dir/$h{$n\
}{name}.aln\";\n	    \n	    $h{$n}{out}=$f[1];\n	 \
   \n	    translate_fasta_seq (\"uU\", \"X\",$h{$n\
}{in}, $h{$n}{NFin});\n	    $n++;\n	  }\n      }\n\
    close (F);\n    \n    \n    chdir ($dir);\n   \
 dump_nf (\"nf\");\n    dump_config ();\n   \n    \
#system (\"nextflow run nf  --name \\'*.seq4nf\\' \
>/dev/null 2>/dev/null\");\n    system (\"nextflow\
 run nf  --name \\'*.seq4nf\\'\");\n    foreach my\
 $n (keys (%h))\n      {\n	translate_fasta_seq (\"\
uU\", \"X\",$h{$n}{NFout},$h{$n}{out});\n      }\n\
    chdir ($cdir);\n  }\nsub seq2msa\n    {\n     \
 my ($in, $out)=@_;\n      my $cdir=getcwd;\n     \
 \n      \n      if (!($in=~/\\//)){$in=$cdir.\"/\\
".$in;}\n      if (!($out=~/\\//)){$out=$cdir.\"/\\
".$out;}\n      \n      my $file=random_string(10)\
;\n      $file=\"$tmpdir/$file\";\n      open (F, \
\">$file\");\n      print F \"$in $out\\n\";\n    \
  close (F);\n      \n      return listseq2listmsa\
 ($file);\n    }\n	\nsub seq2msa_old\n  {\n    my \
($in, $out)=@_;\n    my $cdir = getcwd;\n    my $d\
ir=random_string(10);\n    $dir=\"/tmp/upp.nf4tcof\
fee/$dir\";\n    my $seq=random_string(10);\n    $\
seq.=\".fa\";\n    my $aln=$seq;\n    $aln.=\".aln\
\";\n    \n    mkdir ($dir);\n    translate_fasta_\
seq (\"uU\", \"X\",$in, \"$dir/$seq\");\n    chdir\
 ($dir);\n    \n    dump_nf (\"nf\");\n    dump_co\
nfig ();\n    print \"IN: $in OUT: $cdir/$out\\nDI\
R: $dir\\nnextflow run nf  --name \\'*.fa\\' \\n\"\
;\n    system (\"nextflow run nf  --name \\'*.fa\\\
' \");\n    print \"$dir/$aln $cdir/$out\\n\";\n  \
  translate_fasta_seq (\"xX\", \"U\",$aln, \"$cdir\
/$out\");\n    chdir ($cdir);\n   } \nsub translat\
e_fasta_seq\n  {\n    my ($from, $to, $in, $out)=@\
_;\n    my $n;\n    my $skip;\n    my $l;\n    my \
$cseq;\n    if (!-e $in){return;}\n    \n    open \
(IN, \"$in\");\n    open (OUT, \">$out\");\n   \n \
   while (<IN>)\n      {\n	$l=$_;\n	if ($l=~\">\")\
{$n++;$cseq=\"\";}\n	else { $l=~s/[$from]/$to/;$cs\
eq.=$l;}\n\n	if ($skip){$skip=0;}\n	elsif ($l=~/>f\
ake_seq/){$skip=1;}\n	else\n	  {\n	    print OUT \\
"$l\";\n	  }\n      }\n    if ($n==2 && $from eq \\
"uU\")\n      {\n	print OUT \">fake_seq\\n$cseq\\n\
\";\n      }\n    close (IN);\n    close (OUT);\n \
 }\n\nsub dump_config\n    {\n      open (F, \">ne\
xtflow.config\");\n\n      print F \"docker.enable\
d = true\\n\";\n      print F \"process.container \
= \\'cbcrg/benchfam_large_scale\\'\\n\";\n      cl\
ose (F);\n    }\n\nsub dump_nf\n  {\n    my $nff=s\
hift;\n    open (F,\">$nff\");\n    print F \"#!/u\
sr/bin/env nextflow\\n\";\n    print F \"params.ba\
se_dir=\\\"./\\\"\\n\";\n    print F \"params.out_\
dir=\\\"./\\\"\\n\";\n    print F \"Channel.fromPa\
th(params.name)\\n\";\n    print F \"\\t.map{ tupl\
e(it.baseName, it) }\\n\";\n    \n    print F \"\\\
t.set{ file_names_1 }\\n\";\n    print F \"process\
 clustalo_align{\\n\";\n    print F \"\\tpublishDi\
r params.out_dir, mode: \\\"copy\\\"\\n\";\n    pr\
int F \"tag \\\"\\${name}\\\"\";\n    print F \"\\\
n\";\n    print F \"\\tinput:\\n\";\n    print F \\
"\\tset name, file(seq_file) from file_names_1\\n\\
";\n    print F \"\\toutput:\\n\";\n    print F \"\
\\tfile \\\"\\${name}.aln\\\"\\n\";\n    print F \\
"\\n\";\n    print F \" \\\"\\\"\\\"\\n\";\n    pr\
int F \" clustalo -i \\$seq_file -o \\${name}.aln\\
\n\";\n    print F \"\\\"\\\"\\\"\\n\\n\";\n    pr\
int F \"}\\n\";\n    close (F);\n  }\n\nsub random\
_string\n    {\n      my $len=shift;\n      my @ch\
ars = (\"A\"..\"Z\", \"a\"..\"z\");\n      my $str\
ing;\n      $string .= $chars[rand @chars] for 1..\
$len;\n      return $string;\n    }\n\nsub mymkdir\
\n      {\n	my $d=shift;\n	my $cd='/';\n	\n	foreac\
h my $e (split (/\\//, $d))\n	  {\n	    $cd.=\"$e/\
\";\n	    if ( !-d $cd){mkdir ($cd);}\n	  }\n	retu\
rn;\n      }\n      \n			  \n      \n","\nmy $msaf\
=\"msa.in.tmp.$$\";\nmy $msaoutf=\"msa.out.tmp.$$\\
";\nmy $cost=\"blosum62.tmp.$$\";\n\nopen  (F, $AR\
GV[0]);\nopen  (OUT, \">$msaf\");\n$nseq=0;\nwhile\
 (<F>)\n  {\n    $l=$_;\n    if ( $l=~/^>(\\S+)/)\\
n      {\n	my $simple=\"Seq$nseq\";\n	$s=$seqname{\
$nseq++}=$1;\n	$translate{$simple}=$s;\n	\n	print \
OUT \">$simple\\n\";\n	\n      }\n    else\n      \
{\n	$l=uc($l);\n	print OUT \"$l\";\n      }\n  }\n\
close (F);\nclose(OUT);\n\ndump_blosum ($cost);\ns\
ystem (\"dca -c $cost -q $msaf> $msaoutf 2>/dev/nu\
ll\");\nopen (F, \"$msaoutf\");\nopen (OUT, \">$AR\
GV[1]\");\n\n$read=0;\nwhile (<F>)\n  {\n    $l=$_\
;\n    if ($l=~/^>(\\S+)/)\n      {\n	$read=1;\n	$\
name=$translate{$1};\n	print OUT \">$name\\n\";\n \
     }\n    elsif ($read && ($l=~/\\S/))\n      {\\
n	print OUT \"$l\";\n      }\n    else\n      {\n	\
$read=0;\n      }\n  }\nclose (F);\n\nunlink ($cos\
t);\nunlink ($msaf);\nunlink ($msaoutf);\n\nsub du\
mp_blosum\n  {\n    my $f=shift;\n    open (F, \">\
$f\");\n\n    print F \"6\\n\";\n    print F \"- -\
   0\\n\";\n    print F \"W W   0\\n\";\n    print\
 F \"Y Y   4\\n\";\n    print F \"F F   5\\n\";\n \
   print F \"V V   7\\n\";\n    print F \"L L   7\\
\n\";\n    print F \"I I   7\\n\";\n    print F \"\
M M   6\\n\";\n    print F \"K K   6\\n\";\nprint \
F \"R R   6\\n\";\n    print F \"H H   3\\n\";\n  \
  print F \"Q Q   6\\n\";\n    print F \"E E   6\\\
n\";\n    print F \"D D   5\\n\";\n    print F \"N\
 N   5\\n\";\n    print F \"G G   5\\n\";\n    pri\
nt F \"A A   7\\n\";\n    print F \"P P   4\\n\";\\
n    print F \"T T   6\\n\";\n    print F \"S S   \
7\\n\";\n    print F \"C C   2\\n\";\n    print F \
\"- C  10 \\n\";\n    print F \"- S  10\\n\";\n   \
 print F \"- T  10 \\n\";\n    print F \"- P  10\\\
n\";\n    print F \"- A  10 \\n\";\n    print F \"\
- G  10\\n\";\n    print F \"- N  10 \\n\";\n    p\
rint F \"- D  10\\n\";\n    print F \"- E  10 \\n\\
";\n    print F \"- Q  10\\n\";\nprint F \"- H  10\
 \\n\";\n    print F \"- R  10\\n\";\n    print F \
\"- K  10 \\n\";\n    print F \"- M  10\\n\";\n   \
 print F \"- I  10 \\n\";\n    print F \"- L  10\\\
n\";\n    print F \"- V  10 \\n\";\n    print F \"\
- F  10\\n\";\n    print F \"- Y  10 \\n\";\n    p\
rint F \"- W  10\\n\";\n    print F \"W C  13 \\n\\
";\n    print F \"W S  14\\n\";\n    print F \"W T\
  13 \\n\";\n    print F \"W P  15\\n\";\n    prin\
t F \"W A  14 \\n\";\n    print F \"W G  13\\n\";\\
n    print F \"W N  15 \\n\";\n    print F \"W D  \
15\\n\";\n    print F \"W E  14 \\n\";\n    print \
F \"W Q  13\\n\";\n    print F \"W H  13 \\n\";\n \
   print F \"W R  14\\n\";\n    print F \"W K  14 \
\\n\";\n    print F \"W M  12\\n\";\n    print F \\
"W I  14 \\n\";\n    print F \"W L  13\\n\";\n    \
print F \"W V  14 \\n\";\n    print F \"W F  10\\n\
\";\n    print F \"W Y   9 \\n\";\n    print F \"Y\
 C  13\\n\";\n    print F \"Y S  13 \\n\";\n    pr\
int F \"Y T  13\\n\";\n    print F \"Y P  14 \\n\"\
;\n    print F \"Y A  13\\n\";\n    print F \"Y G \
 14 \\n\";\n    print F \"Y N  13\\n\";\n    print\
 F \"Y D  14 \\n\";\n    print F \"Y E  13\\n\";\n\
    print F \"Y Q  12 \\n\";\n    print F \"Y H   \
9\\n\";\n    print F \"Y R  13 \\n\";\n    print F\
 \"Y K  13\\n\";\n    print F \"Y M  12 \\n\";\n  \
  print F \"Y I  12\\n\";\n    print F \"Y L  12 \\
\n\";\n    print F \"Y V  12\\n\";\n    print F \"\
Y F   8 \\n\";\n    print F \"F C  13\\n\";\nprint\
 F \"F S  13 \\n\";\n    print F \"F T  13\\n\";\n\
    print F \"F P  15 \\n\";\n    print F \"F A  1\
3\\n\";\n    print F \"F G  14 \\n\";\n    print F\
 \"F N  14\\n\";\n    print F \"F D  14 \\n\";\n  \
  print F \"F E  14\\n\";\n    print F \"F Q  14 \\
\n\";\n    print F \"F H  12\\n\";\n    print F \"\
F R  14 \\n\";\n    print F \"F K  14\\n\";\n    p\
rint F \"F M  11 \\n\";\n    print F \"F I  11\\n\\
";\n    print F \"F L  11 \\n\";\n    print F \"F \
V  12\\n\";\n    print F \"V C  12 \\n\";\n    pri\
nt F \"V S  13\\n\";\n    print F \"V T  11 \\n\";\
\n    print F \"V P  13\\n\";\n    print F \"V A  \
11 \\n\";\n    print F \"V G  14\\n\";\n    print \
F \"V N  14 \\n\";\n    print F \"V D  14\\n\";\np\
rint F \"V E  13 \\n\";\nprint F \"V Q  13\\n\";\n\
print F \"V H  14 \\n\";\nprint F \"V R  14\\n\";\\
nprint F \"V K  13 \\n\";\nprint F \"V M  10\\n\";\
\nprint F \"V I   8 \\n\";\nprint F \"V L  10\\n\"\
;\nprint F \"L C  12 \\n\";\nprint F \"L S  13\\n\\
";\nprint F \"L T  12 \\n\";\nprint F \"L P  14\\n\
\";\nprint F \"L A  12 \\n\";\nprint F \"L G  15\\\
n\";\nprint F \"L N  14 \\n\";\nprint F \"L D  15\\
\n\";\nprint F \"L E  14 \\n\";\nprint F \"L Q  13\
\\n\";\nprint F \"L H  14 \\n\";\nprint F \"L R  1\
3\\n\";\nprint F \"L K  13 \\n\";\nprint F \"L M  \
 9\\n\";\nprint F \"L I   9 \\n\";\nprint F \"I C \
 12\\n\";\nprint F \"I S  13 \\n\";\nprint F \"I T\
  12\\n\";\nprint F \"I P  14 \\n\";\nprint F \"I \
A  12\\n\";\nprint F \"I G  15 \\n\";\nprint F \"I\
 N  14\\n\";\nprint F \"I D  14 \\n\";\nprint F \"\
I E  14\\n\";\nprint F \"I Q  14 \\n\";\nprint F \\
"I H  14\\n\";\nprint F \"I R  14 \\n\";\nprint F \
\"I K  14\\n\";\nprint F \"I M  10 \\n\";\nprint F\
 \"M C  12\\n\";\nprint F \"M S  12 \\n\";\nprint \
F \"M T  12\\n\";\nprint F \"M P  13 \\n\";\nprint\
 F \"M A  12\\n\";\nprint F \"M G  14 \\n\";\nprin\
t F \"M N  13\\n\";\nprint F \"M D  14 \\n\";\npri\
nt F \"M E  13\\n\";\nprint F \"M Q  11 \\n\";\npr\
int F \"M H  13\\n\";\nprint F \"M R  12 \\n\";\np\
rint F \"M K  12\\n\";\nprint F \"K C  14 \\n\";\n\
print F \"K S  11\\n\";\nprint F \"K T  12 \\n\";\\
nprint F \"K P  12\\n\";\nprint F \"K A  12 \\n\";\
\nprint F \"K G  13\\n\";\nprint F \"K N  11 \\n\"\
;\nprint F \"K D  12\\n\";\nprint F \"K E  10 \\n\\
";\nprint F \"K Q  10\\n\";\nprint F \"K H  12 \\n\
\";\nprint F \"K R   9\\n\";\nprint F \"R C  14 \\\
n\";\nprint F \"R S  12\\n\";\nprint F \"R T  12 \\
\n\";\nprint F \"R P  13\\n\";\nprint F \"R A  12 \
\\n\";\nprint F \"R G  13\\n\";\nprint F \"R N  11\
 \\n\";\nprint F \"R D  13\\n\";\nprint F \"R E  1\
1 \\n\";\nprint F \"R Q  10\\n\";\nprint F \"R H  \
11 \\n\";\nprint F \"H C  14\\n\";\nprint F \"H S \
 12 \\n\";\nprint F \"H T  13\\n\";\nprint F \"H P\
  13 \\n\";\nprint F \"H A  13\\n\";\nprint F \"H \
G  13 \\n\";\nprint F \"H N  10\\n\";\nprint F \"H\
 D  12 \\n\";\nprint F \"H E  11\\n\";\nprint F \"\
H Q  11 \\n\";\nprint F \"Q C  14\\n\";\nprint F \\
"Q S  11 \\n\";\nprint F \"Q T  12\\n\";\nprint F \
\"Q P  12 \\n\";\nprint F \"Q A  12\\n\";\nprint F\
 \"Q G  13 \\n\";\nprint F \"Q N  11\\n\";\nprint \
F \"Q D  11 \\n\";\nprint F \"Q E   9\\n\";\nprint\
 F \"E C  15 \\n\";\nprint F \"E S  11\\n\";\nprin\
t F \"E T  12 \\n\";\nprint F \"E P  12\\n\";\npri\
nt F \"E A  12 \\n\";\nprint F \"E G  13\\n\";\npr\
int F \"E N  11 \\n\";\nprint F \"E D   9\\n\";\np\
rint F \"D C  14 \\n\";\nprint F \"D S  11\\n\";\n\
print F \"D T  12 \\n\";\nprint F \"D P  12\\n\";\\
nprint F \"D A  13 \\n\";\nprint F \"D G  12\\n\";\
\nprint F \"D N  10 \\n\";\nprint F \"N C  14\\n\"\
;\nprint F \"N S  10 \\n\";\nprint F \"N T  11\\n\\
";\nprint F \"N P  13 \\n\";\nprint F \"N A  13\\n\
\";\nprint F \"N G  11 \\n\";\nprint F \"G C  14\\\
n\";\nprint F \"G S  11 \\n\";\nprint F \"G T  13\\
\n\";\nprint F \"G P  13 \\n\";\nprint F \"G A  11\
\\n\";\nprint F \"A C  11 \\n\";\nprint F \"A S  1\
0\\n\";\nprint F \"A T  11 \\n\";\nprint F \"A P  \
12\\n\";\nprint F \"P C  14 \\n\";\nprint F \"P S \
 12\\n\";\nprint F \"P T  12 \\n\";\nprint F \"T C\
  12\\n\";\nprint F \"T S  10 \\n\";\nprint F \"S \
C  12\\n\";\nclose (F);\n    return;\n  }\n    \n"\
,"\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env \
qw(USER);\n\n                                     \
                   \nuse strict;                  \
                           \nuse warnings;\nuse di\
agnostics;\n\nmy $in_hit_list, my $in_aln=0, my(%n\
ame_list)=(),my (%list)=(),my $n_seq=0; my $test=0\
;\nmy($j)=0, my $n=0, my $nom, my $lg_query, my %v\
u=();\n\nopen (F, \">tmp\");\n\n$/=\"\\n\";\nwhile\
 (<>)\n{\n    print F $_;\n    if($_ =~ /Query=\\s\
*(.+?)\\s/i) { $nom=$1;}\n\n    if ( /Sequences pr\
oducing significant alignments/){$in_hit_list=1;}\\
n    \n    if ($_=~ /^pdb\\|/i) { $_=~ s/pdb\\|//g\
; }\n    if ($_=~ /^(1_\\d+)\\s+\\d+/) { $_=~ s/$1\
/QUERY/;}\n      \n    if ( /^(\\S+).+?\\s+[\\d.]+\
\\s+([\\de.-]+)\\s+$/ && $in_hit_list)	\n    {\n	m\
y($id)=$1; # \n	$id=~ s/\\|/_/g; #\n	if ($id =~ /.\
+_$/) { chop($id) }; #\n	$name_list{$n_seq++}=$id;\
\n	$name_list{$n_seq-1}=~ s/.*\\|//g;     \n    }\\
n  \n    if (/query/i) {$in_aln=1;}\n    if ( /^(\\
\S+)\\s+(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)/ || /^(\\\
S+)(\\s+)(\\-+)(\\s+)/ && ($in_aln == 1))\n    {\n\
	my $name=$1;\n	my $start=$2;\n	my $seq=$3;\n	my $\
end=$4;\n		\n	if ($name =~ /QUERY/i) { $lg_query=l\
ength($seq); }\n\n	unless ($test > $n) #m\n	{\n	  \
  my(@seqq)= split('',$seq);\n	    my($gap_missing\
)= scalar(@seqq);\n	    \n	    while ($gap_missing\
 != $lg_query)  { unshift (@seqq,\"-\"); $gap_miss\
ing= scalar(@seqq); }\n	    $seq=join('',@seqq);  \
#m\n	}\n	\n	if ($name =~ /QUERY/i)\n	{\n	    $n=0;\
 %vu=(); $j=0;\n	    $list{$n}{'real_name'}=\"$nom\
\";\n	}	\n	else\n	{\n	    unless (exists $vu{$name\
}) { ++$j;}	\n	    $list{$n}{'real_name'}=$name_li\
st{$j-1};\n	}\n		\n	$list{$n}{'name'}=$name;\n\n	$\
seq=~tr/a-z/A-Z/;\n	$list{$n}{'seq'}=$list{$n}{'se\
q'};\n	$list{$n}{'seq'}.=$seq;\n\n	$n++;\n	$vu{$na\
me}++;\n	$test++;\n   } \n    \n}\n\nmy @numero=()\
;\n\nfor (my $a=0; $a<$n; $a++) #m\n{\n    my $lon\
g=length($list{0}{'seq'});  \n    my $long1= lengt\
h($list{$a}{'seq'});\n  \n    while ($long1 ne $lo\
ng)\n    {\n	$list{$a}{'seq'}.=\"-\";\n	$long1= le\
ngth ($list{$a}{'seq'});\n    } \n \n    push (@nu\
mero,\"$list{$a}{'name'} $list{$a}{'real_name'}\\n\
\");\n}\n\nmy %dejavu=();\n\n\nfor (my $i=0; $i<=$\
#numero; $i++)\n{\n    my $s=\">$list{$i}{'real_na\
me'}\\n$list{$i}{'seq'}\\n\";\n    my $k=0;\n    \\
n    if (exists $dejavu{$numero[$i]}) {next;}\n   \
 else\n    {	\n	for ($j=0; $j<$n ; $j++)\n	{\n	   \
 if (\"$numero[$i]\" eq \"$numero[$j]\" && $j != $\
i )\n	    {\n		++$k;\n		$s .=\">$list{$j}{'real_na\
me'}\\n$list{$j}{'seq'}\\n\";\n	    }\n	}	\n    }\\
n    \n    if ($k>0) \n    {\n	my $cons;\n	open (S\
OR,\">tempo_aln2cons\"); print SOR $s;  close SOR \
;\n	open (COM,\"t_coffee -other_pg seq_reformat -i\
n tempo_aln2cons -action +aln2cons +upper |\") ; \\
n     	while (<COM>)\n	{	\n	    if (/^>/) { $cons \
=\">$list{$i}{'real_name'}\\n\"; next;}\n	    $_=~\
 s/\\n//g;\n	    $cons .=$_;\n	}\n	close COM; unli\
nk (\"tempo_aln2cons\");\n	print $cons,\"\\n\"; pr\
int F $cons,\"\\n\";\n    }	\n    else  { print $s\
;  print F $s; }\n    \n    $dejavu{$numero[$i]}++\
;\n} #m\n\nexit;\n\n\n\n\n\n\n\n\n\n\n\n","use Env\
;\n\n\n$tmp_dir=\"\";\n$init_dir=\"\";\n$program=\\
"tc_generic_method.pl\";\n\n$blast=@ARGV[0];\n\n$n\
ame=\"query\";$seq=\"\";\n%p=blast_xml2profile($na\
me,$seq,100, 0, 0, $blast);\n&output_profile (%p);\
\n\n\nsub output_profile\n  {\n    my (%profile)=(\
@_);\n    my ($a);\n    for ($a=0; $a<$profile{n};\
 $a++)\n      {\n	\n	print \">$profile{$a}{name} $\
profile{$a}{comment}\\n$profile{$a}{seq}\\n\";\n  \
    }\n    return;\n  }\nsub file_contains \n  {\n\
    my ($file, $tag, $max)=(@_);\n    my ($n);\n  \
  $n=0;\n    \n    if ( !-e $file && ($file =~/$ta\
g/)) {return 1;}\n    elsif ( !-e $file){return 0;\
}\n    else \n      {\n	open (FC, \"$file\");\n	wh\
ile ( <FC>)\n	  {\n	    if ( ($_=~/$tag/))\n	     \
 {\n		close (FC);\n		return 1;\n	      }\n	    els\
if ($max && $n>$max)\n	      {\n		close (FC);\n		r\
eturn 0;\n	      }\n	    $n++;\n	  }\n      }\n   \
 close (FC);\n    return 0;\n  }\n	    \n	  \nsub \
file2string\n  {\n    my $f=@_[0];\n    my $string\
, $l;\n    open (F,\"$f\");\n    while (<F>)\n    \
  {\n\n	$l=$_;\n	#chomp ($l);\n	$string.=$l;\n    \
  }\n    close (F);\n    $string=~s/\\r\\n//g;\n  \
  $string=~s/\\n//g;\n    return $string;\n  }\n\n\
\n\nsub tag2value \n  {\n    \n    my $tag=(@_[0])\
;\n    my $word=(@_[1]);\n    my $return;\n    \n \
   $tag=~/$word=\"([^\"]+)\"/;\n    $return=$1;\n \
   return $return;\n  }\n      \nsub hit_tag2pdbid\
\n  {\n    my $tag=(@_[0]);\n    my $pdbid;\n     \
  \n    $tag=~/id=\"(\\S+)\"/;\n    $pdbid=$1;\n  \
  $pdbid=~s/_//;\n    return $pdbid;\n  }\nsub id2\
pdbid \n  {\n    my $id=@_[0];\n  \n    if ($id =~\
/pdb/)\n      {\n	$id=~/pdb(.*)/;\n	$id=$1;\n     \
 }\n    $id=~s/[|�_]//g;\n    return $id;\n  }\nsu\
b set_blast_type \n  {\n    my $file =@_[0];\n    \
if (&file_contains ($file,\"EBIApplicationResult\"\
,100)){$BLAST_TYPE=\"EBI\";}\n    elsif (&file_con\
tains ($file,\"NCBI_BlastOutput\",100)) {$BLAST_TY\
PE=\"NCBI\";}\n    else\n      {\n	$BLAST_TYPE=\"\\
";\n      }\n    return $BLAST_TYPE;\n  }\nsub bla\
st_xml2profile \n  {\n    my ($name,$seq,$maxid, $\
minid, $mincov, $file)=(@_);\n    my (%p, $a, $str\
ing, $n);\n    \n\n\n    if ($BLAST_TYPE eq \"EBI\\
" || &file_contains ($file,\"EBIApplicationResult\\
",100)){%p=ebi_blast_xml2profile(@_);}\n    elsif \
($BLAST_TYPE eq \"NCBI\" || &file_contains ($file,\
\"NCBI_BlastOutput\",100)){%p=ncbi_blast_xml2profi\
le(@_);}\n    else \n      {\n	print \"***********\
* ERROR: Blast Returned an unknown XML Format ****\
******************\";\n	die;\n      }\n    for ($a\
=0; $a<$p{n}; $a++)\n      {\n	my $name=$p{$a}{nam\
e};\n	$p{$name}{seq}=$p{$a}{seq};\n      }\n    re\
turn %p;\n  }\nsub ncbi_blast_xml2profile \n  {\n \
   my ($name,$seq,$maxid, $minid, $mincov, $string\
)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhits,@identi\
fyerL);\n    \n    \n    $seq=~s/[^a-zA-Z]//g;\n  \
  $L=length ($seq);\n    \n    %hit=&xml2tag_list \
($string, \"Hit\");\n    \n    \n    for ($nhits=0\
,$a=0; $a<$hit{n}; $a++)\n      {\n	my ($ldb,$id, \
$identity, $expectation, $start, $end, $coverage, \
$r);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\"\";\n\n	%ID=\
&xml2tag_list ($hit{$a}{body}, \"Hit_id\");\n	$ide\
ntifyer=$ID{0}{body};\n	\n	%DE=&xml2tag_list ($hit\
{$a}{body}, \"Hit_def\");\n	$definition=$DE{0}{bod\
y};\n	\n	%HSP=&xml2tag_list ($hit{$a}{body}, \"Hsp\
\");\n	for ($b=0; $b<$HSP{n}; $b++)\n	  {\n	    my\
 (%START,%END,%E,%I,%Q,%M);\n\n	 \n	    %START=&xm\
l2tag_list ($HSP{$b}{body}, \"Hsp_query-from\");\n\
	    %HSTART=&xml2tag_list ($HSP{$b}{body}, \"Hsp_\
hit-from\");\n	    \n	    %LEN=  &xml2tag_list ($H\
SP{$b}{body}, \"Hsp_align-len\");\n	    %END=  &xm\
l2tag_list ($HSP{$b}{body}, \"Hsp_query-to\");\n	 \
   %HEND=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_hi\
t-to\");\n	    %E=&xml2tag_list     ($HSP{$b}{body\
}, \"Hsp_evalue\");\n	    %I=&xml2tag_list     ($H\
SP{$b}{body}, \"Hsp_identity\");\n	    %Q=&xml2tag\
_list     ($HSP{$b}{body}, \"Hsp_qseq\");\n	    %M\
=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_hseq\");\
\n	    \n	    for ($e=0; $e<$Q{n}; $e++)\n\n	     \
 {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{body};\n		if\
 ($seq eq\"\"){$seq=$qs;$L=length($seq);}\n		\n		$\
expectation=$E{$e}{body};\n		$identity=($LEN{$e}{b\
ody}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;\n		$st\
art=$START{$e}{body};\n		$end=$END{$e}{body};\n		$\
Hstart=$HSTART{$e}{body};\n		$Hend=$HEND{$e}{body}\
;\n	\n		$coverage=(($end-$start)*100)/$L;\n\n	\n		\
if ($identity>$maxid || $identity<$minid || $cover\
age<$mincov){next;}\n		@lr1=(split (//,$qs));\n		@\
lr2=(split (//,$ms));\n		$l=$#lr1+1;\n		for ($c=0;\
$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		for ($d=0,$c\
=0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$c];\n		   \
 if ( $r=~/[A-Za-z]/)\n		      {\n			\n			$p[$nhit\
s][$d + $start-1]=$lr2[$c];\n			$d++;\n		      }\n\
		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$ms;\\
n		$QstartL[$nhits]=$start;\n		$HstartL[$nhits]=$H\
start;\n		$identityL[$nhits]=$identity;\n		$endL[$\
nhits]=$end;\n		$definitionL[$nhits]=$definition;\\
n		$identifyerL[$nhits]=$identifyer;\n		$comment[$\
nhits]=\"$ldb|$identifyer [Eval=$expectation][id=$\
identity%][start=$Hstart end=$Hend]\";\n		$nhits++\
;\n	      }\n	  }\n      }\n    \n    $profile{n}=\
0;\n    $profile{$profile{n}}{name}=$name;\n    $p\
rofile{$profile{n}}{seq}=$seq;\n    $profile {n}++\
;\n    \n    for ($a=0; $a<$nhits; $a++)\n      {\\
n	$n=$a+1;\n	\n	$profile{$n}{name}=\"$name\\_$a\";\
\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{Qseq}=$Q\
seq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	$profil\
e{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}{Hstart}\
=$HstartL[$a];\n	$profile{$n}{identity}=$identityL\
[$a];\n	$profile{$n}{definition}=$definitionL[$a];\
\n	$profile{$n}{identifyer}=$identifyerL[$a];\n	$p\
rofile{$n}{comment}=$comment[$a];\n	for ($b=0; $b<\
$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	      {\n	\
	$profile{$n}{seq}.=$p[$a][$b];\n	      }\n	    el\
se\n	      {\n		$profile{$n}{seq}.=\"-\";\n	      \
}\n	  }\n      }\n    \n    $profile{n}=$nhits+1;\\
n    return %profile;\n  }\nsub ebi_blast_xml2prof\
ile \n  {\n    my ($name,$seq,$maxid, $minid, $min\
cov, $string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$n\
hits,@identifyerL,$identifyer);\n    \n\n    \n   \
 $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq);\n   \
 %hit=&xml2tag_list ($string, \"hit\");\n    \n   \
 for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\n	\
my ($ldb,$id, $identity, $expectation, $start, $en\
d, $coverage, $r);\n	my (%Q,%M,%E,%I);\n	\n	$ldb=&\
tag2value ($hit{$a}{open}, \"database\");\n	$ident\
ifyer=&tag2value ($hit{$a}{open}, \"id\");\n\n	$de\
scription=&tag2value ($hit{$a}{open}, \"descriptio\
n\");\n	\n	%Q=&xml2tag_list ($hit{$a}{body}, \"que\
rySeq\");\n	%M=&xml2tag_list ($hit{$a}{body}, \"ma\
tchSeq\");\n	%E=&xml2tag_list ($hit{$a}{body}, \"e\
xpectation\");\n	%I=&xml2tag_list ($hit{$a}{body},\
 \"identity\");\n	\n\n	for ($b=0; $b<$Q{n}; $b++)\\
n	  {\n	    \n	    \n	    $qs=$Q{$b}{body};\n	    \
$ms=$M{$b}{body};\n	    if ($seq eq\"\"){$seq=$qs;\
$L=length($seq);}\n\n	    $expectation=$E{$b}{body\
};\n	    $identity=$I{$b}{body};\n	    \n	    	   \
 \n	    $start=&tag2value ($Q{$b}{open}, \"start\"\
);\n	    $end=&tag2value ($Q{$b}{open}, \"end\");\\
n	    $startM=&tag2value ($M{$b}{open}, \"start\")\
;\n	    $endM=&tag2value ($M{$b}{open}, \"end\");\\
n	    $coverage=(($end-$start)*100)/$L;\n	    \n	 \
  # print \"$id: ID: $identity COV: $coverage [$st\
art $end]\\n\";\n	    \n	    \n	    if ($identity>\
$maxid || $identity<$minid || $coverage<$mincov){n\
ext;}\n	    # print \"KEEP\\n\";\n\n	    \n	    @l\
r1=(split (//,$qs));\n	    @lr2=(split (//,$ms));\\
n	    $l=$#lr1+1;\n	    for ($c=0;$c<$L;$c++){$p[$\
nhits][$c]=\"-\";}\n	    for ($d=0,$c=0; $c<$l; $c\
++)\n	      {\n		$r=$lr1[$c];\n		if ( $r=~/[A-Za-z\
]/)\n		  {\n		    \n		    $p[$nhits][$d + $start-1\
]=$lr2[$c];\n		    $d++;\n		  }\n	      }\n	  \n	 \
   \n	    $identifyerL[$nhits]=$identifyer;\n	    \
$comment[$nhits]=\"$ldb|$identifyer [Eval=$expecta\
tion][id=$identity%][start=$startM end=$endM]\";\n\
	    $nhits++;\n	  }\n      }\n    \n    $profile{\
n}=0;\n    $profile{$profile{n}}{name}=$name;\n   \
 $profile{$profile{n}}{seq}=$seq;\n    $profile {n\
}++;\n    \n    for ($a=0; $a<$nhits; $a++)\n     \
 {\n	$n=$a+1;\n	$profile{$n}{name}=\"$name\\_$a\";\
\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{identify\
er}=$identifyerL[$a];\n	\n	$profile{$n}{comment}=$\
comment[$a];\n	for ($b=0; $b<$L; $b++)\n	  {\n	   \
 if ($p[$a][$b])\n	      {\n		$profile{$n}{seq}.=$\
p[$a][$b];\n	      }\n	    else\n	      {\n		$prof\
ile{$n}{seq}.=\"-\";\n	      }\n	  }\n      }\n   \
 $profile{n}=$nhits+1;\n    \n    return %profile;\
\n  }\n\nsub blast_xml2hit_list\n  {\n    my $stri\
ng=(@_[0]);\n    return &xml2tag_list ($string, \"\
hit\");\n  }\nsub xml2tag_list  \n  {\n    my ($st\
ring_in,$tag)=@_;\n    my $tag_in, $tag_out;\n    \
my %tag;\n    \n    if (-e $string_in)\n      {\n	\
$string=&file2string ($string_in);\n      }\n    e\
lse\n      {\n	$string=$string_in;\n      }\n    $\
tag_in1=\"<$tag \";\n    $tag_in2=\"<$tag>\";\n   \
 $tag_out=\"/$tag>\";\n    $string=~s/>/>##1/g;\n \
   $string=~s/</##2</g;\n    $string=~s/##1/<#/g;\\
n    $string=~s/##2/#>/g;\n    @l=($string=~/(\\<[\
^>]+\\>)/g);\n    $tag{n}=0;\n    $in=0;$n=-1;\n  \
\n \n\n    foreach $t (@l)\n      {\n\n	$t=~s/<#//\
;\n	$t=~s/#>//;\n	\n	if ( $t=~/$tag_in1/ || $t=~/$\
tag_in2/)\n	  {\n	 \n	    $in=1;\n	    $tag{$tag{n\
}}{open}=$t;\n	    $n++;\n	    \n	  }\n	elsif ($t=\
~/$tag_out/)\n	  {\n	    \n\n	    $tag{$tag{n}}{cl\
ose}=$t;\n	    $tag{n}++;\n	    $in=0;\n	  }\n	els\
if ($in)\n	  {\n	   \n	    $tag{$tag{n}}{body}.=$t\
;\n	  }\n      }\n  \n    return %tag;\n  }\n\n\n\\
n\n","use Env qw(HOST);\nuse Env qw(HOME);\nuse En\
v qw(USER);\nwhile (<>)\n  {\n    if ( /^>(\\S+)/)\
\n      {\n	if ($list{$1})\n	  {\n	    print \">$1\
_$list{$1}\\n\";\n	    $list{$1}++;\n	  }\n	else\n\
	  {\n	    print $_;\n	    $list{$1}=1;\n	  }\n   \
   }\n    else\n      {\n	print $_;\n      }\n  }\\
n      \n","\n\n\nuse Env qw(HOST);\nuse Env qw(HO\
ME);\nuse Env qw(USER);\n\n\nopen (F,$ARGV[0]);\nw\
hile ( <>)\n  {\n    @x=/([^:,;\\)\\(\\s]+):[^:,;\\
\)\\(]*/g;\n    @list=(@list,@x);\n  }\n$n=$#list+\
1;\nforeach $n(@list){print \">$n\\nsequence\\n\";\
}\n\n\nclose (F);\n","\nopen (F, $ARGV[0]);\n\nwhi\
le ( <F>)\n  {\n    @l=($_=~/(\\S+)/g);\n    \n   \
 $name=shift @l;\n    \n    print STDOUT \"\\n>$na\
me\\n\";\n    foreach $e (@l){$e=($e eq \"0\")?\"O\
\":\"I\";print \"$e\";}\n  }\nclose (F);\n\n		    \
   \n    \n","use strict;\nuse FileHandle;\nuse En\
v qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\\
nmy %name;\nmy $nseq;\nmy $F= new FileHandle;\nope\
n ($F, $ARGV[0]);\nwhile(<$F>)\n  {\n    \n    my \
$l=$_;\n    if ($l=~/^#/){;}\n    elsif (($l=~/\\d\
+\\s+\\d+\\s+(\\S+)\\s+(\\S+)/))\n      {\n	my $na\
me=$1;\n	my $seq=$2;\n	print \">$name\\n$seq\\n\";\
\n      }\n  }\nclose ($F);\nexit (0);\n\n\n","use\
 Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER\
);\n\n$tmp=\"$ARGV[0].$$\";\nopen (IN, $ARGV[0]);\\
nopen (OUT, \">$tmp\");\n\nwhile ( <IN>)\n  {\n   \
 $file=$_;\n    $file=~s/\\r\\n/\\n/g;\n    $file=\
~s/\\n\\r/\\n/g;\n    $file=~s/\\r\\r/\\n/g;\n    \
$file=~s/\\r/\\n/g;\n    print OUT \"$file\";\n  }\
\nclose (IN);\nclose (OUT);\n\nopen (OUT, \">$ARGV\
[0]\");\nopen (IN, \"$tmp\");\n\nwhile ( <IN>)\n{\\
n  print OUT \"$_\";\n}\nclose (IN);\nclose (OUT);\
\nunlink ($tmp);\n\n"};
