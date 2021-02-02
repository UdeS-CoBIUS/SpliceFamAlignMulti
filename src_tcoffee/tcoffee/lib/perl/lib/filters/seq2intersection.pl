#!/usr/bin/env perl
use strict;
use FileHandle;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
my %name;
my $nseq;
my $fasta;
if ($ARGV[2] eq "-fasta"){$fasta=1;}
my $F= new FileHandle;

open ($F, $ARGV[1]);
while(<$F>)
  {
    my $l=$_;
    if ($l=~/^#/){;}
    elsif (($l=~/\d+\s+\d+\s+(\S+)\s+(\S+)/))
      {
	my $n=$1;
	$name{$1}++;
      }
  }
close ($F);

open ($F, $ARGV[0]);
while(<$F>)
  {
    my $l=$_;
    if ($l=~/^#/){;}
    elsif ($l=~/\d+\s+\d+\s+(\S+)\s+(\S+)/)
      {
	my $n=$1;
	$name{$n}++;
	if ($name{$n}==2){$nseq++;}
      }
  }
close ($F);

if (!$fasta && $nseq>0)
  {
    print "#NAMESEQ_01\n";
    print "# $nseq\n";
  }
open ($F, $ARGV[0]);
while(<$F>)
  {
    my $l=$_;
    if ($l=~/^#/){;}
    elsif ($l=~/.\d+\s+\d+\s+(\S+)\s+(\S+)/)
      {
	my $n=$1;
	my $s=$2;
	if ($name{$n}==2)
	  {
	    if ($fasta)
	      {
		print ">$n\n$s\n";
	      }
	    else
	      {
		print "$l";
	      }
	  }
      }
  }
close ($F);
exit (0);


