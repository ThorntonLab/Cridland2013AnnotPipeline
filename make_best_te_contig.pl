#!perl -w

use strict;
my $input = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
my $tes = shift(@ARGV) or die; #must use the output_trimmed_restricted_complete file - keep straight first pass and second pass event ids
#find the annotated breakpoint, identify which contig that is and then find the TE contig that matches and select the best TE

my %posid = ();
my %posstart = ();
my %posstop = ();
my %posch = ();

open(X, "<$tes");
while(my $line = <X>) {
    chomp $line;
    my @x = split(/\t/, $line);
    my @y = split(/:/, $x[3]);
    $posid{$y[0]} = $y[0];
    $posch{$y[0]} = $x[0];
    $posstart{$y[0]} = $x[1];
    $posstop{$y[0]} = $x[2];

}
close X;
my %te = ();
my %bestlen = ();
my %beste = ();
my %contig = ();
open(B, "<$input");
while(my $line2 = <B>) {

    chomp $line2;
    my @b = split(/\t/, $line2);
    my @c = split(/\./, $b[0]); 
  #  print $c[1], "\t", $b[1], "\n";
  #  print $posch{$c[1]}, "\n";
    if($posch{$c[1]} == $b[1]) { #on the right chrom
	if(($posstart{$c[1]} == ($b[8] + 1)) or ($posstart{$c[1]} == ($b[9] + 1)) or ($posstop{$c[1]} == ($b[8] + 1)) or ($posstop{$c[1]} == ($b[9] + 1))) {#annotated contig start or stop matches the chromosome.output file position - it is the correct contig
	    print "Found contig \n";
	    $contig{$b[0]} = $b[0]; #store the whole first string
	    
	}
    
	if(!(exists($te{$b[0]}))) {
	    
	    $te{$b[0]} = $b[1];
	    $bestlen{$b[0]} = $b[3];
	    $beste{$b[0]} = $b[11];
	}elsif(exists($te{$b[0]})){
	    
	    if($b[3] > $bestlen{$b[0]}) {
		if($b[11] < $beste{$b[0]}) { 
		    
		    $te{$b[0]} = $b[1];
		    $bestlen{$b[0]} = $b[3];
		    $beste{$b[0]} = $b[11];
		}
	    }
	}
    }
}
close B;
open(A, ">>$output");
while((my $key, my $value) = each(%te)) {

    if(exists($contig{$key})) {
	print A $key, "\t", $value, "\n";
    }
}
close A;
