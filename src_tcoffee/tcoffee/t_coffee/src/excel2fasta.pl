#!/usr/bin/env perl

open (F, $ARGV[0]);

while ( <F>)
  {
    @l=($_=~/(\S+)/g);
    
    $name=shift @l;
    
    print STDOUT "\n>$name\n";
    foreach $e (@l){$e=($e eq "0")?"O":"I";print "$e";}
  }
close (F);

		       
    
# Jeu 19 nov 2020 09:17:32 EST
