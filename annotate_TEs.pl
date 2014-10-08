#!perl -w

use strict;

#this path may need to be altered
###########################################################
my $ref = 'TE_positions';
###########################################################

my $dir = shift(@ARGV) or die;
my $left_contig_file = shift(@ARGV) or die;
my $pos_est = shift(@ARGV) or die; #the estimated position of the TE insertion from teclust
my $chrom_num = shift(@ARGV) or die; #the chrom_num file
my $output = shift(@ARGV) or die; #the standard output file
#The standard output is in the following format

#Not currently printing the BED file
#my $bed_out = shift(@ARGV) or die;  #the bed formatted output file

$left_contig_file =~ s/$dir//go;

my $event = $left_contig_file;

$event =~ s/\.left\.fasta\.contigs//;

my $id = -1;

if($event =~ m/short/) {

    my @event = split(/\./, $event);

    $id = $event[2];

}else{

    my @event = split(/\./, $event);

    $id = $event[1];
}

#make the rest of the file names
$left_contig_file = $dir . $left_contig_file;
my $left_chrom_file = $left_contig_file . ".chromosome.output";

my $right_contig_file = $left_contig_file;
$right_contig_file =~ s/left/right/go;

my $right_chrom_file = "";
if(-e $right_contig_file) {

    $right_chrom_file = $right_contig_file . ".chromosome.output";

}else{

    if ($left_contig_file =~ m/short/) {

	$right_contig_file =~ s/short\.//go;

	if(-e $right_contig_file) {

	    $right_chrom_file = $right_contig_file . ".chromosome.output";

	}else {

	    die;

	}
    }elsif ($left_contig_file !~ m/short/) {

	$right_contig_file =~ s/$dir//go;

	$right_contig_file = $dir . "short." . $right_contig_file;

	if(-e $right_contig_file) {

	    $right_chrom_file = $right_contig_file . ".chromosome.output";

	}else {

	    die;

	}
    }

}

my %contig = (); #keeps track of contigs -stores length of contig
my %novel_breakpoint = (); #keeps track of contigs - determines if they contain a breakpoint or not
my %shared_breakpoint = ();
my %contig_chrom = (); #does the contig match a unique position
my %contig_te = (); #does the contig match a TE
my %contig_unique = ();
my %unique_start = ();
my %unique_stop = ();
my %multi_start = ();
my %multi_stop = ();
my %contig_rev = ();
my %shared_name = ();
my %novel_name = ();

#this keeps track of the string/number for each chromosome
my %chrom = ();

#these keep track of TEs in the reference
my %te = ();
my %techrom = ();
my %testart = ();
my %testop = ();
my %tename = ();
my %tefreq = ();

#get list of TE positions in reference

open(A, "<$ref");

my $count = 0;

while(my $line = <A>) {
    
    chomp $line;
    
    my @x = split(/\t/, $line);
    
    $te{$count} = $count;
    $techrom{$count} = $x[0];
    $testart{$count} = $x[1];
    $testop{$count} = $x[2];
    $tename{$count} = $x[3];
    $tefreq{$x[3]} = $x[4];
    ++$count;
}

close A;

my @contig_files = ($left_contig_file, $right_contig_file); #find the lengths of the contigs

foreach my $contig_file (@contig_files) {
    
    my $contig_side = -1;
    
    open(CONTIG, "<$contig_file");
    
    if ($contig_file =~ m/left/) {
	
	$contig_side = 0;
	
    }elsif ($contig_file =~ m/right/) {
	
	$contig_side = 1;	
    }
    
    my $current = 0;
    my $length = 0;
    
  LOOP:while(my $contig_line = <CONTIG>) {
      
      chomp $contig_line;
      
      if($contig_line =~ m/Contig/) { #if it is a new contig - set length to 0 set contig value to contig add to hash
	  
	  my @a = split(/\./, $contig_line);
	  
	  $current = $contig_side . ":" . $a[4];
	  $contig{$current} = 0;
	  $shared_breakpoint{$current} = 0;
	  $novel_breakpoint{$current} = 0;
	  $contig_chrom{$current} = -1;
	  $contig_unique{$current} = -1;
	  $contig_te{$current} = -1;
	  $unique_start{$current} = -1;
	  $unique_stop{$current} = -1;
	  $multi_start{$current} = -1;
	  $multi_stop{$current} = -1;
	  $contig_rev{$current} = -1;
	  $shared_name{$current} = "";
	  $length = 0;
	  next LOOP;
      }
      
      if($contig_line !~ m/Contig/) { #else if not a new contig continue to add length to current contig
	  
	  $length = $length + length($contig_line);
	  $contig{$current} = $length;
      }
  }
    
    close CONTIG;
}

open(CHROM, "<$chrom_num"); #get the number / chrom name conversion

while( my $ch_line = <CHROM>) {

    chomp $ch_line;

    $ch_line =~ s/\s+/\t/;
    
    my @ch = split(/\t/, $ch_line);

    $chrom{$ch[0]} = $ch[1];

}
close CHROM;

#now open the blast output files and annotate each position in the contigs

my @chrom_files = ($left_chrom_file, $right_chrom_file);

my $shared = 0;
my $TE = 0;

foreach my $file (@chrom_files) {
    
    my $distance = 500;
    my $side = -1;
    
    my @z = split(/\./, $file);
    
    my $est_chrom = 0;
    my $event_id = $z[1];
    my $est1 = 0;
    my $est2 = 0;
#find the estimated position of the insertion according to teclust
    open(E, "<$pos_est");
    
    while(my $posline = <E>) {
	
	chomp $posline;
	
	my @g = split(/\t/, $posline);
	
	if ($g[1] == $event_id){
	    
	    $est1 = $g[2];
	    $est2 = $g[3];
	    $est_chrom = $g[0];	    
	}
    }
    
    close E;
    
    if ($file =~ m/left/) {
	
	$side = 0;
	
    }elsif ($file =~ m/right/) {
	
	$side = 1;
	
    }
#read the output file and examine each HSP    
    open(B, "<$file"); 
    
    my $sstart = 0;
    my $sstop = 0;
    
  HSPLOOP:while(my $line2 = <B>) { #one HSP
      
      chomp $line2;
      
      my @e = split(/\t/, $line2);
      
      my @f = split(/\./, $e[0]);
      
      if($e[3] < 100) {
	  
	  next HSPLOOP; #not considering anything that is too short
	  
      }
      if($e[2] < 90) {
	  
	  next HSPLOOP;
      }
      
      if($e[8] < $e[9]) {
	  
	  $sstart = $e[8];
	  $sstop = $e[9];
	  $contig_rev{$side . ":" . $f[4]} = 0;
	  
      }else {
	  
	  $sstart = $e[9];
	  $sstop = $e[8];
	  $contig_rev{$side . ":" . $f[4]} = 0;
      } 
#first - discover if the HSP matches the region identified by teclust
      if( (($e[1] == $est_chrom) and (abs($sstart - $est1) < $distance)) or (($e[1] == $est_chrom) and (abs($sstop - $est1) < $distance)) or (($e[1] == $est_chrom) and (abs($sstart - $est2) < $distance)) or (($e[1] == $est_chrom) and (abs($sstop - $est2) < $distance)) or (($sstart > $est1) and ($sstart < $est2) and ($e[1] == $est_chrom)) or (($sstop > $est1) and ($sstop < $est2) and ($e[1] == $est_chrom))){#found a match
	  
	  $contig_chrom{$side . ":" . $f[4]} = $e[1];#record the chromosome
	  $contig_unique{$side . ":" . $f[4]} = 1;#mark contig as having unique sequence
	  
	  if($multi_start{$side . ":" . $f[4]} == -1) {
	      
#no multi yet
	      
	  }elsif($multi_start{$side . ":" . $f[4]} != -1) {
	      
#already a multi - check to see if the new unique overrides the multi
	      
	      if($e[2] >= 98) {
		  
		  $multi_start{$side . ":" . $f[4]} = -1;
		  $multi_stop{$side . ":" . $f[4]} = -1;
		  $contig_te{$side . ":" . $f[4]} = -1;
		  
	      }
	      
	  }
	  
	  if($unique_start{$side . ":" . $f[4]} == -1) {
	      
	      $unique_start{$side . ":" . $f[4]} = $e[6];
	      $unique_stop{$side . ":" . $f[4]} = $e[7];
	      
	      if($side == 0) {#left side 
		  
		  $novel_breakpoint{$side . ":" . $f[4]} = $sstop;
		  
	      }elsif($side == 1) {#right side
		  
		  $novel_breakpoint{$side . ":" . $f[4]} = $sstart;
		  
	      }
	      
	  }elsif($unique_start{$side . ":" . $f[4]} != -1){#this is a gap - still want to get the whole unique region
	      
	      if(($side == 0) and ($contig_rev{$side . ":" . $f[4]} == 0)) { #left e[8] < e[9]
		  
		  if($e[6] < $unique_start{$side . ":" . $f[4]}) {
		      
		      $unique_start{$side . ":" . $f[4]} = $e[6];
		      $novel_breakpoint{$side . ":" . $f[4]} = $sstop;
		      
		  }
		  if($e[7] > $unique_stop{$side . ":" . $f[4]}) {
		      
		      $unique_stop{$side . ":" . $f[4]} = $e[7];
		      
		  }
		  
	      }elsif(($side == 0) and ($contig_rev{$side . ":" . $f[4]} == 1)) {#left e[8] > e[9]
		  
		  if($e[6] < $unique_start{$side . ":" . $f[4]}) {
		      
		      $unique_start{$side . ":" . $f[4]} = $e[6];
		      
		  }
		  if($e[7] > $unique_stop{$side . ":" . $f[4]}) {
		      
		      $unique_stop{$side . ":" . $f[4]} = $e[7];
		      $novel_breakpoint{$side . ":" . $f[4]} = $sstop;
		      
		  }
		  
	      }elsif(($side == 1) and ($contig_rev{$side . ":" . $f[4]} == 0)) {#right e[8] < e[9]
		  
		  if($e[6] < $unique_start{$side . ":" . $f[4]}) {
		      
		      $unique_start{$side . ":" . $f[4]} = $e[6];
		      $novel_breakpoint{$side . ":" . $f[4]} = $sstart;
		      
		  }
		  if($e[7] > $unique_stop{$side . ":" . $f[4]}) {
		      
		      $unique_stop{$side . ":" . $f[4]} = $e[7];
		      
		  }
		  
	      }elsif(($side == 1) and ($contig_rev{$side . ":" . $f[4]} == 1)) {#right e[8] > e[9]
		  
		  if($e[6] < $unique_start{$side . ":" . $f[4]}) {
		      
		      $unique_start{$side . ":" . $f[4]} = $e[6];
		      
		      
		  }
		  if($e[7] > $unique_stop{$side . ":" . $f[4]}) {
		      
		      $unique_stop{$side . ":" . $f[4]} = $e[7];
		      $novel_breakpoint{$side . ":" . $f[4]} = $sstart;
		      
		  }
		  
	      }	  
	  }
	  
#second - if the contig matches the region - see if it contains any shared insertions
	  
	SHAREDLOOP:while((my $key, my $value) = each (%te)) { 
	    
	    if(($techrom{$key} == $e[1]) and ($est_chrom == $e[1])) {#if the te is on the right chrom
		
		if(((($sstart) >= ($testart{$key} + 25)) and (($sstop) <= ($testart{$key} - 25))) 
		   or ((($sstop) >= ($testart{$key} + 25)) and (($sstart) <= ($testart{$key} - 25))) 
		   or ((($testop{$key} -  25) >= $sstart) and (($testop{$key}  + 25) <= ($sstop))) 
		   or ((($testop{$key} - 25) >= ($sstop)) and (($testop{$key} + 25) <= $sstart))) { #if HSP and TE overlap at all - found a shared TE - contig must span the start or stop position of the shared TE by 25 bp on either side
		    
		    $shared_breakpoint{$e[1] . "\t" . ($testart{$key} -1)} = $side . ":" . $f[4];
		    $shared_name{$e[1] . "\t" . ($testart{$key} -1)} = $tename{$key};
		    
		    if(($sstart >= ($testart{$key} -1)) and ($sstop <= ($testop{$key} -1))) {
			    $contig_unique{$side . ":" . $f[4]} = -1; #have to reset this if unique contig is identified as a shared
		    }
		    
		    if(($side == 0) and ($sstop > ($testop{$key} -1))) {
			
		    }else{
			
			$contig_unique{$side . ":" . $f[4]} = -1; #have to reset this if unique contig is identified as a shared and the contig is entirely within a TE range
			
		    }
		    
		    if(($side == 1) and ($sstart < ($testart{$key} -1))) {
			
		    }else{
			
			$contig_unique{$side . ":" . $f[4]} = -1; #have to reset this if unique contig is identified as a shared and the contig is entirely within a TE range
			
		    }
		    if(!(exists( $shared_breakpoint{$e[1] . "\t" . ($testart{$key} -1)}))){
			
			$shared_breakpoint{$e[1] . "\t" . ($testart{$key} -1)} = $side . ":" . $f[4];
			++$shared;
			$shared_name{$e[1] . "\t" . ($testart{$key} -1)} = $tename{$key};
			if(($sstart >= ($testart{$key} -1)) and ($sstop <= ($testop{$key} -1))) {
			    $contig_unique{$side . ":" . $f[4]} = -1; #have to reset this if unique contig is identified as a shared
			}
			
			if(($side == 0) and ($sstop > ($testop{$key} -1))) {
			    
			}else{
			    
			    $contig_unique{$side . ":" . $f[4]} = -1; #have to reset this if unique contig is identified as a shared and the contig is entirely within a TE range
			    
			}
			
			if(($side == 1) and ($sstart < ($testart{$key} -1))) {
			    
			}else{
			    
			    $contig_unique{$side . ":" . $f[4]} = -1; #have to reset this if unique contig is identified as a shared and the contig is entirely within a TE range
			    
			}
		    }
		}elsif(((($testart{$key} - 25) < $sstart) and (($testop{$key} + 25) > $sstop)) or ((($testart{$key} - 25) < $sstop) and (($testop{$key} + 25) > $sstart)) ){
		    $contig_te{$side . ":" . $f[4]} = 1; #HSP is a TE
		    $contig_unique{$side . ":" . $f[4]} = -1; #HSP is a TE - matches to the right region, but does not span breakpoint therefore it is not going to be called a shared TE - this is because it could be from somewhere else - must have unique upstream or downstream data to confirm as shared
		}elsif((($sstart < ($testop{$key} + 1000)) and ($sstart > $testop{$key}) and ($sstop > $testop{$key})) or (($sstart < ($testart{$key}) and ($sstop > ($testart{$key} - 1000)) and ($sstop < $testart{$key})))) {
		    $contig_te{$side . ":" . $f[4]} = -1;
		    $contig_unique{$side . ":" . $f[4]} = -1;	 #this contig is near the predicted area which is also near, but not crossing a shared TE boundary - it may either be an absence - which can be confirmed later or a shared TE with not enough sequence crossing the position to call it   
		}
	    }else {
		
		next SHAREDLOOP;
		
	    }
	}
	  
	  next HSPLOOP; #found a unique contig - go to next HSP
     }else{#the HSP does not match the area in question - now is it a TE at all? - also make sure only to call this if this region of the contig has not already been identified as unique

#first check to see if the HSP is a subset of another HSP that was previously annotated as unique - this arises when there are shared insertions as part of the HSP aligns to many regions - however a small amount of overlap may be due to a few base pairs being the same in the TE and unique regions for a novel event - so a little overlap is ok
	 
	 $multi_start{$side . ":" . $f[4]} = $e[6];
	 $multi_stop{$side . ":" . $f[4]} = $e[7];
	 
	 if(($multi_start{$side . ":" . $f[4]} >= $unique_start{$side . ":" . $f[4]}) and ($multi_stop{$side . ":" . $f[4]} <= $unique_stop{$side . ":" . $f[4]})){ #the HSP is within a unique HSP - skip
	     $multi_start{$side . ":" . $f[4]} = -1;
	     $multi_stop{$side . ":" . $f[4]} = -1;
	     next HSPLOOP;
	     
	 }elsif(($unique_start{$side . ":" . $f[4]} < $multi_start{$side . ":" . $f[4]}) and ($side == 0) and ($contig_rev{$side . ":" . $f[4]} == 0)) { #if there is too much overlap between fragments - skip
	     
	     if(($unique_stop{$side . ":" . $f[4]} - $multi_start{$side . ":" . $f[4]}) > 25) {
		 $multi_start{$side . ":" . $f[4]} = -1;
		 $multi_stop{$side . ":" . $f[4]} = -1;
		 next HSPLOOP;
		 
	     }
	 }elsif(($multi_start{$side . ":" . $f[4]} < $unique_start{$side . ":" . $f[4]}) and ($side == 0) and ($contig_rev{$side . ":" . $f[4]} == 1)) {
	     
	     if(($multi_stop{$side . ":" . $f[4]} - $unique_start{$side . ":" . $f[4]}) > 25) {
		 $multi_start{$side . ":" . $f[4]} = -1;
		 $multi_stop{$side . ":" . $f[4]} = -1;
		 next HSPLOOP;
		 
	     }
	 }elsif(($unique_start{$side . ":" . $f[4]} > $multi_start{$side . ":" . $f[4]}) and ($side == 1) and ($contig_rev{$side . ":" . $f[4]} == 1)) { 
	     
	     if(($multi_stop{$side . ":" . $f[4]} - $unique_start{$side . ":" . $f[4]}) > 25) {
		 $multi_start{$side . ":" . $f[4]} = -1;
		 $multi_stop{$side . ":" . $f[4]} = -1;
		 next HSPLOOP;
		 
	     }
	 }elsif(($multi_start{$side . ":" . $f[4]} > $unique_start{$side . ":" . $f[4]}) and ($side == 1) and ($contig_rev{$side . ":" . $f[4]} == 0)) { 
	     
	     if(($unique_stop{$side . ":" . $f[4]} - $multi_start{$side . ":" . $f[4]}) > 25) {
		 $multi_start{$side . ":" . $f[4]} = -1;
		 $multi_stop{$side . ":" . $f[4]} = -1;
		 next HSPLOOP;
		 
	     }
	 }	  
	 
       TELOOP:while((my $key2, my $value2) = each(%te)) {
	   
	   if((($sstart >= ($testart{$key2} -1)) and ($sstart <= ($testop{$key2} -1))) or (($sstop >= ($testart{$key2} -1)) and ($sstop <= ($testop{$key2} -1))) or ((($testart{$key2} -1) >= $sstart) and (($testart{$key2} -1) <= $sstop)) or ((($testop{$key2} -1) >= $sstart) and (($testop{$key2} -1) <= $sstop))) { #if HSP and TE overlap at all
	       
	       if (($e[3] >= 100) and ($e[2] > 90)) {
		   
		   $contig_te{$side . ":" . $f[4]} = 1; #HSP is a TE
		   
		   push (@{$novel_name{$side . ":" . $f[4]}}, $tename{$key2} . ":" . $e[10]);

	       }else{
		   
		   next TELOOP;
		   
	       }	    
	   }  
       }
     }
  }
}
open(BED, ">>$bed_out");
open(OUT, ">>$output");
    
#loop over all contigs searching for insertions and print out    
my $left_novel = 0;
my $right_novel = 99999999999999999999;	
my $novelchrom = -1; 
my $novelleft = 0;
my $novelright = 0;
my $contig_id = 0;
my $novel = 0;
my @left_novel = ();
my @right_novel = ();
my $left_type = 0;
my $right_type = 0;

while(( my $contig_key, my $contig_value) = each (%contig)) {#for each contig
    
    if(($contig_te{$contig_key} == 1) and ($contig_unique{$contig_key} == 1)) {#this finds contigs that contains breakpoints

	my @nov_side = split(/:/, $contig_key);

	my $nov_side = $nov_side[0];

	if($nov_side == 0) {
	
	    $novelleft = 1;#we found a novel insertion
	    $contig_id = $contig_key;#mark
	    $novelchrom = $contig_chrom{$contig_key};#can just assign b/c this was confirmed earlier

	    $left_type = 1;

	    if($novel_breakpoint{$contig_key} > $left_novel) {#largest number should be the breakpoint
		
		$left_novel = $novel_breakpoint{$contig_key}; #position of the breakpoint

		foreach my $left_list (@{$novel_name{$contig_key}}) {

		    push @left_novel, $left_list;
		
		}
	    }

	}

	if($nov_side == 1) {
	
	    $novelright = 1;#we found a novel insertion
	    $contig_id = $contig_key;#mark
	    $novelchrom = $contig_chrom{$contig_key};#can just assign b/c this was confirmed earlier

	    $right_type = 1;
	
	    if($novel_breakpoint{$contig_key} < $right_novel) {#largest number should be the breakpoint
		
		$right_novel = $novel_breakpoint{$contig_key}; #position of the breakpoint
	
		foreach my $right_list (@{$novel_name{$contig_key}}) {

		    push @right_novel, $right_list;
		
		}
	    }

	}

    }
}

my $anyteleft = 0;
my $anyuniqueleft = 0;

if($novelleft == 0) {#did not find a contig with the breakpont = see if there are contigs w/te and w/unique

    while((my $key1, my $value1) = each (%contig)) {

	my @sideleft = split(/:/, $key1);

	my $current_left = $sideleft[0];

	if($current_left == 0) {

	    if($contig_te{$key1} == 1) {

		$anyteleft = 1;

		foreach my $left_list2 (@{$novel_name{$key1}}) {

		    push @left_novel, $left_list2;
		
		}

	    }
	    if($contig_unique{$key1} == 1) {
		
		$anyuniqueleft = 1;
		$novelchrom = $contig_chrom{$key1};

		if($left_novel == 0) {

		    $left_novel = $novel_breakpoint{$key1};

		}elsif($left_novel != 0) {

		    if($novel_breakpoint{$key1} > $left_novel) {

			$left_novel = $novel_breakpoint{$key1};	
			
		    }

		}
	    }
	}
    }
	
    if(($anyteleft == 1) and ($anyuniqueleft == 1)) {
	
	$novelleft = 1;
	$left_type = 2;
    }
	    
}

my $anyteright = 0;
my $anyuniqueright = 0;

if($novelright == 0) {#did not find a contig with the breakpont = see if there are contigs w/te and w/unique   

    while((my $key2, my $value2) = each (%contig)) {
	
	my @sideright = split(/:/, $key2);

	my $current_right = $sideright[0];
	
	if($current_right == 1) {

	    if($contig_te{$key2} == 1) {

		$anyteright = 1;

		foreach my $right_list2 (@{$novel_name{$key2}}) {

		    push @right_novel, $right_list2;
		
		}		

	    }
	    if($contig_unique{$key2} == 1) {

		$anyuniqueright = 1;
		$novelchrom = $contig_chrom{$key2};

		if($right_novel == 99999999999999999999) {
		    
		    $right_novel = $novel_breakpoint{$key2};
		    
		}elsif($right_novel != 99999999999999999999) {
		    
		    if($novel_breakpoint{$key2} < $right_novel) {
			
			$right_novel = $novel_breakpoint{$key2};	
			
		    }
		    
		}
	    }
	}
    }

    if(($anyteright == 1) and ($anyuniqueright == 1)) {
	
	$novelright = 1;
	$right_type = 2;
    }
        
}

my $left_string = "empty_left";
my $right_string = "empty_right";

if(@left_novel) {

    my %left_hits = ();
    my %left_norm = ();
    my $leftsumx = 0;
    my %left_prob = ();
    my %left_prob_adj = ();
    my $leftprobsum = 0;

    foreach my $lhits (@left_novel) {#count hits and get the sum of frequencies for TEs hit

	my @lhits = split(/:/, $lhits);

	if (!(exists($left_hits{$lhits[0]}))) {

	    if($lhits[1] == 0) {

		$lhits[1] = 1e-300;
	    }

	    $left_hits{$lhits[0]} = 1/$lhits[1];
	    $leftsumx = $leftsumx + $tefreq{$lhits[0]}; #sum frequencies the first time each name is seen

	}elsif(exists($left_hits{$lhits[0]})) {

	    if($lhits[1] == 0) {

		$lhits[1] = 1e-300;
	    }

	    $left_hits{$lhits[0]} = $left_hits{$lhits[0]} + (1/$lhits[1]);

	}
	
    }

    while((my $lh_key, my $lh_value) = each(%left_hits)) { #should contain each name once

	$left_norm{$lh_key} = $tefreq{$lh_key} / $leftsumx;

    }
 
    while((my $lp_key, my $lp_value) = each(%left_hits)) {
	
	$left_prob{$lp_key} = $lp_value * $left_norm{$lp_key};

	$leftprobsum = $leftprobsum + $left_prob{$lp_key};

    }

    while(( my $lp2_key, my $lp2_value) = each(%left_prob)) {

	$left_prob_adj{$lp2_key} = $left_prob{$lp2_key} / $leftprobsum;

	if($left_prob_adj{$lp2_key} > 0.05) {
	    
	    $left_string = $left_string . $left_prob_adj{$lp2_key} . ";" . $lp2_key . ",";
	}

    }

}else{

    $left_string = "empty left";

}


if(@right_novel) {

    my %right_hits = ();
    my %right_norm = ();
    my $rightsumx = 0;
    my %right_prob = ();
    my %right_prob_adj = ();
    my $rightprobsum = 0;

    foreach my $rhits (@right_novel) {#count hits and get the sum of frequencies for TEs hit

	my @rhits = split(/:/, $rhits);

	if (!(exists($right_hits{$rhits[0]}))) {

	    if($rhits[1] == 0) {

		$rhits[1] = 1e-300;
	    }

	    $right_hits{$rhits[0]} = 1/$rhits[1];
	    $rightsumx = $rightsumx + $tefreq{$rhits[0]}; #sum frequencies the first time each name is seen

	}elsif(exists($right_hits{$rhits[0]})) {

	    if($rhits[1] == 0) {

		$rhits[1] = 1e-300;
	    }

	    $right_hits{$rhits[0]} = $right_hits{$rhits[0]} + (1/$rhits[1]);

	}
	
    }

    while((my $rh_key, my $rh_value) = each(%right_hits)) { #should contain each name once

	$right_norm{$rh_key} = $tefreq{$rh_key} / $rightsumx;

    }
    
    while((my $rp_key, my $rp_value) = each(%right_hits)) {
	
	$right_prob{$rp_key} = $rp_value * $right_norm{$rp_key};

	$rightprobsum = $rightprobsum + $right_prob{$rp_key};

    }

    while(( my $rp2_key, my $rp2_value) = each(%right_prob)) {

	$right_prob_adj{$rp2_key} = $right_prob{$rp2_key} / $rightprobsum;

	if($right_prob_adj{$rp2_key} > 0.05) {
	    
	    $right_string = $right_string . $right_prob_adj{$rp2_key} . ";" . $rp2_key . ",";
	
	}

    }



}else{

    $right_string = "empty_right";

}

if(($novelleft == 1) and ($novelright == 1)){ #this will go if both of the sides has a contig that spans the breakpoint - add one to all positions to correct that our files have chromosomes start at 0 and flybase etc start at position 1

    if(($left_string ne "empty_left") and ($right_string ne "empty_right")) {
    
	print OUT $novelchrom, "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t", $id, ":",  $left_string, "!", $right_string,  "\t", "novel", "\t", $left_type, "\t", $right_type, "\n";
	
	print BED $chrom{$novelchrom}, "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t",$id, ":",  $left_string, "!", $right_string,  "\t", "1000", "\t", "+", "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t", "\n";
    }
}
if(($novelleft == 1) and ($anyuniqueright == 1) and ($novelright == 0)){ #this will go if both of the sides has a contig that spans the breakpoint
    
    if(($left_string ne "empty_left") and ($right_string ne "empty_right")) {   

	print OUT $novelchrom, "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t",$id, ":",  $left_string, "!", $right_string, "\t", "novel", "\t", $left_type, "\t", $right_type, "\n";
    
	print BED $chrom{$novelchrom}, "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t",$id, ":",  $left_string, "!", $right_string, "\t", "1000", "\t", "+", "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t", "\n";
    }
}

if(($anyuniqueleft == 1) and ($novelleft == 0) and ($novelright == 1)){ #this will go if both of the sides has a contig that spans the breakpoint
    
    if(($left_string ne "empty_left") and ($right_string ne "empty_right")) {

	print OUT $novelchrom, "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t", $id, ":",  $left_string, "!", $right_string,  "\t", "novel",  "\t", $left_type, "\t", $right_type, "\n";
	
	print BED $chrom{$novelchrom}, "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t", $id, ":",  $left_string, "!", $right_string, "\t", "1000", "\t", "+", "\t", ($left_novel + 1), "\t", ($right_novel + 1), "\t", "\n";
	
    }
}


while(( my $shared_key, my $shared_value) = each(%shared_breakpoint)) { #this prints out the starting locations of shared insertions
    
    if ($shared_value ne 0) {
	
	my @sh = split(/\t/, $shared_key);
	
	print OUT $sh[0], "\t", ($sh[1] + 1), "\t", ($sh[1] + 1), "\t", $id, ":", $shared_name{$shared_key}, "\t", "shared", "\n";
	
	print BED $chrom{$sh[0]}, "\t", ($sh[1] + 1), "\t", ($sh[1] + 1), "\t", $id, ":", $shared_name{$shared_key}, "\t", "1000", "\t", "+", "\t", $sh[1], "\t", $sh[1], "\n";
    }
}


close OUT;
close BED;
