#!/usr/bin/env perl
#This script reads an Interleaved MSA and Outputs a Fasta_msa
#if two sequences have the same name, the second one is renamed name_1 and so on
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
# VersionTag  1.72

$tmp="$ARGV[0].$$";
open (IN, $ARGV[0]);
open (OUT, ">$tmp");

while ( <IN>)
  {
    $file=$_;
    $file=~s/\r\n/\n/g;
    $file=~s/\n\r/\n/g;
    $file=~s/\r\r/\n/g;
    $file=~s/\r/\n/g;
    print OUT "$file";
  }
close (IN);
close (OUT);

open (OUT, ">$ARGV[0]");
open (IN, "$tmp");

while ( <IN>)
{
  print OUT "$_";
}
close (IN);
close (OUT);
unlink ($tmp);

