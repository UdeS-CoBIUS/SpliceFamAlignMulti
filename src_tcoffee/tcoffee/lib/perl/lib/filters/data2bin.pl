#!/usr/bin/env perl
use HTTP::Date;
use strict;
use FileHandle;
my $p;



$p=process_param (@ARGV);
if ($p->{file}=~/\s/)
  {
    my $list=$p->{file};
    my @l=split (/\s+/,$list);
    foreach my $f (@l)
      {
	$p->{file}=$f;
	process_record ($p);
      }
  }
else
  {
    process_record ($p);
  }
die;


sub process_record
  {
    my $p=shift;
    
    my $data;
    my $LOG =new FileHandle;
    my $DATA=new FileHandle;
    
    if ( !-e $p->{file})
      {
	print STDERR "*** ERROR: $p->{file} Does not exist\n";
	return 0;
      }
    
    else
      {
	print STDERR "*** Process Record: $p->{file}\n";
      }

    $p->{out}=$p->{file};
    $p->{out}=~s/\.[^\.]*$//;
    
    $p->{log}="$p->{out}.log";
    $p->{outfile}="$p->{out}.binned";
    
    open ($LOG, ">$p->{log}");
    open ($DATA, ">$p->{outfile}");
    
    
    
    $data=undump_data ($p->{file});
    $data=data2normalize($data, $p);
    $data=data2bin($data, $p);
    display_stat  ($data,$p,$LOG);
    display_param ($p, $LOG);
    
    dump_data($data, $DATA);
    close ($DATA);
    close ($LOG);
    print STDERR "\tLOGFILE:  $p->{log}\n";
    print STDERR "\tDATAFILE: $p->{outfile}\n";
    return 1;
  }
sub file2field_list
  {
    my $file=shift;
    my $F =new FileHandle;
    my $f={};
    my $nf;
    
    open ($F, $file);
    while (<$F>)
      {
	my $v=$_;
	chomp($v);
	$f->{$nf++}=$v;
      }
    return $f;
    close ($F);
  }
sub display_stat 
  {
    my $d=shift;
    my $p=shift;
    my $F=shift;
    
    my $tot;
    my $count={};
     
    foreach my $exp (sort {$a<=>$b}(keys(%$d)))
	  {
	    foreach my $record (sort {$a<=>$b}(keys(%{$d->{$exp}})))
	      {
		$count->{$d->{$exp}{$record}{bin}}{bin}++;
		$count->{$d->{$exp}{$record}{bin}}{score}+=$d->{$exp}{$record}{$p->{bin_field}};
		
		$tot++;
	      }
	  }
    foreach my $b (sort{$a<=>$b}keys (%$count))
      {
	printf $F "BIN: %2d    %10d    %.3f  %10d\n", $b, $count->{$b}{bin}, $count->{$b}{bin}/$tot,$count->{$b}{score}/$count->{$b}{bin};
      }
  }

sub pairplot2field_list
  {
    my $field={};
    my $nf;
    
    $field->{$nf++}="chr";
    $field->{$nf++}="window";
    $field->{$nf++}="Mbp";
    $field->{$nf++}="value::inv";
    $field->{$nf++}="value::wrong_order";
    $field->{$nf++}="value::wrong_distance";
    $field->{$nf++}="value::smaller";
    $field->{$nf++}="value::larger";
    $field->{$nf++}="value::good_pair";
    return $field;
   }
sub dump_data
      {
	my $d=shift;
	my $F=shift;
	
	foreach my $exp (sort {$a<=>$b}(keys(%$d)))
	  {
	    foreach my $record (sort {$a<=>$b}(keys(%{$d->{$exp}})))
	      {
		print $F "#d;$exp;$record;";
		foreach my $k (sort (keys (%{$d->{$exp}{$record}})))
		  {
		    
		    print $F "$k;$d->{$exp}{$record}{$k};";
		  }
		print $F "\n";
	      }
	  }
      }

sub undump_data
  {
    my $file =shift;
    my $field_mode=shift;
    my $F=new FileHandle;
    my $n;
    my $field;
    my $data={};
    
    if ($field_mode eq "pairplot" || !$field_mode)
      {
	$field=pairplot2field_list();
      }
    elsif (-e $field_mode)
      {
	$field=file2field_list($field_mode);
      }
    
    open ($F, $file);
    while (<$F>)
      {
	my $line=$_;
	if (!($line=~/#/))
	  {
	    my @l=($line=~/(\S+)/g);
	    my $exp=$l[0];
	    my $rec=$n++;
	    for (my $a=0;$a<=$#l; $a++)
	      {
		$data->{$exp}{$rec}{$field->{$a}}=$l[$a];
	      }
	    $data->{$exp}{$rec}{'value::nogood_pair'}+=$data->{$exp}{$rec}{'value::inv'};
	    $data->{$exp}{$rec}{'value::nogood_pair'}+=$data->{$exp}{$rec}{'value::wrong_order'};
	    $data->{$exp}{$rec}{'value::nogood_pair'}+=$data->{$exp}{$rec}{'value::wrong_distance'};
	    $data->{$exp}{$rec}{'value::nogood_pair'}+=$data->{$exp}{$rec}{'value::smaller'};
	    $data->{$exp}{$rec}{'value::nogood_pair'}+=$data->{$exp}{$rec}{'value::larger'};
	  }
      }
    
    close ($F);
    return $data;
  }
sub data2bin
    {
      my $d=shift;
      my $p=shift;

      if    ($p->{dim}==1){return data2bin1d($d,$p);}
      elsif ($p->{dim}==2){return data2bin2d($d,$p);}
    }
sub data2bin2d
    {
      my $data=shift;
      my $p=shift;

      my $bin_field1=$p->{bin_field1};
      my $bin_field2=$p->{bin_field2};
      my $nbin1=$p->{nbin1};
      my $nbin2=$p->{nbin2};
      my $bin_interval1=$p->{bin_interval1};
      my $bin_interval2=$p->{bin_interval2};
      foreach my $exp (keys (%$data))
	{
	  foreach my $rec (keys (%{$data->{$exp}}))
	    {
	      my $v1=$data->{$exp}{$rec}{$bin_field1};
	      my $v2=$data->{$exp}{$rec}{$bin_field2};
	      
	      my $bin1=int($v1/$bin_interval1);
	      $bin1=($bin1<$nbin1)?$bin1:$nbin1-1;

	      my $bin2=int($v2/$bin_interval2);
	      $bin2=($bin2<$nbin2)?$bin2:$nbin2-1;
	      
	      my $bin=($bin2*$nbin1)+$bin1;
	      $data->{$exp}{$rec}{bin}=$bin;
	    
	      print "$v1 -> $bin1 $v2 -> $bin2 ===> $bin\n";
	    }
	}
      return $data;
    }
sub data2bin1d
    {
      my $data=shift;
      my $p=shift;

      my $bin_field=$p->{bin_field};
      
      my $bin_mode=$p->{bin_mode};
      my $nbin=$p->{nbin};
      my $bin_interval=$p->{bin_interval};
      
      if ( !$bin_field){$bin_field=$p->{bin_field}="value::good_pair";}
      if ( !$nbin){$nbin=$p->{nbin}=10;}
      if ( !$bin_mode){$bin_mode="fixed_interval";}
      if ( !$bin_interval){$p->{bin_interval}=$bin_interval=100;}
      
      
     
      if ($bin_mode eq "fixed_interval")
	{
	  foreach my $exp (keys (%$data))
	    {
	      foreach my $rec (keys (%{$data->{$exp}}))
		{
		  my $value=$data->{$exp}{$rec}{$bin_field};
		  my $bin=int($value/$bin_interval);
		  if ($nbin!=-1){$bin=($bin<$nbin)?$bin:$nbin-1;}
		  $data->{$exp}{$rec}{bin}=$bin;
		}
	    }
	}
      
      return $data;
    }
    
sub data2normalize
    {
      my $data=shift;
      my $p=shift;
      
      my $t_good;
      my $n_good;
      my $m_good;
      my $norm;
      
      my $count={};
      
      if ($p->{normalize} eq "no"){return $data;}
      elsif (!$p->{normalize} )      {$p->{normalize}="mediane";}
      elsif ( $p->{normalize} eq "1"){$p->{normalize}="mediane";}
      
      if (!$p->{bin_field}){$p->{bin_field}="value::good_pair";}
      
      my $mode= $p->{normalize};
      my $field=$p->{bin_field};
      
      foreach my $exp (keys (%$data))
	{
	  foreach my $rec (keys (%{$data->{$exp}}))
	    {
	      if ($data->{$exp}{$rec}{$field})
		{
		  $count->{$data->{$exp}{$rec}{$field}}++;
		  $t_good+=$data->{$exp}{$rec}{$field};
		  $n_good++;
		}
	    }
	}
      
      if    ($mode eq "mediane"){$norm=(sort {$count->{$a}<=>$count->{$b}}(keys(%$count)))[-1];}
      elsif ($mode eq "averge" ){$norm=$t_good/=$n_good;} 
      
      $p->{normalize_factor}=$norm;
      
      foreach my $exp (keys (%$data))
	{
	  foreach my $rec (keys (%{$data->{$exp}}))
	    {
	       foreach my $f (keys (%{$data->{$exp}{$rec}}))
		 {
		   if ($f=~/value::/)
		     {
		       $data->{$exp}{$rec}{$f}*=1000;
		       $data->{$exp}{$rec}{$f}/=$norm;
		       $data->{$exp}{$rec}{$f}=int($data->{$exp}{$rec}{$f});
		     }
		 }
	     }
	}
      return $data;
    }

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                    Parameters                                    #
#                                                                  #
#                                                                  #
####################################################################

sub check_parameters 
  {
    my $p=shift;
    my $rp={};

    $rp->{file}=1;
    $rp->{bin_mode}=1;
    $rp->{dim}=1;
    
    $rp->{bin_field}=1;
    
    $rp->{bin_interval}=1;
    $rp->{nbin}=1;
    
    $rp->{bin_field1}=1;
    $rp->{bin_interval1}=1;
    $rp->{nbin1}=1;

    $rp->{bin_field2}=1;
    $rp->{bin_interval2}=1;
    $rp->{nbin2}=1;
    
    
    
    $rp->{normalize}=1;
    
    if (!$p)
      {
	foreach my $k (keys (%$rp))
	  {
	    print STDERR "-$k <value> ";
	  }
	print "\n";
	die;
      }
    
    foreach my $k (keys(%$p))
      {
	if (!$rp->{$k})
	  {
	    print STDERR "\n****ERROR: $k is an unknown pararmeter[FATAL]***\n";
	    die;
	  }

      }
    return $p;
  }
sub display_param
  {
    my $p=shift;
    my $F=shift;
    
    print $F "************** PARAM::START *****************\n";
    
    foreach my $v (keys(%$p))
      {
	print $F "PARAM: $v ---> $p->{$v}\n";
      }
    printf $F "Command line:\ndata2bin.pl ";
    foreach my $v (keys(%$p))
      {
	print $F "-$v $p->{$v} ";
      }
    print $F "\n************** PARAM::START *****************\n";
  }
sub process_param
  {
    my @arg=@_;
    my $cl=join(" ", @arg);
    
    my @commands=split (/\s\-+/,$cl);
    my $param={};
    if ($cl)
      {
	foreach my $c (@commands)
	  {
	    if (!($c=~/\S/)){next;}
	    $c=~/(\w+)\s*(.*)\s*/;
	    my $k=$1;
	    if (!$2){$param->{$k}=1;}
	    else {$param->{$k}=$2;}
	    $param->{$k}=~s/\s*$//;
	  }
	return check_parameters ($param);
      }
    else
      {
	check_parameters();
      }
  }
