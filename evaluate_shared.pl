#!perl -w

use strict;

my $input = shift(@ARGV) or die;
my $out = shift(@ARGV) or die;
my $pos = shift(@ARGV) or die;
my $freq = shift(@ARGV) or die;
my $dist = shift(@ARGV) or die;
my %techr = ();
my %tel = ();
my %ter = ();
my %tename = ();
my %name = ();

open(X, "<$freq");
while(my $line3 = <X>) {

    chomp $line3;
    my @x = split(/\t/, $line3);
    $name{$x[0] . ":" . $x[1] . ":" . $x[2]} = $x[3];
}
close X;
open(A, "<$pos");#te_pos_estimate
while(my $line = <A>) {

    chomp $line;
    my @a = split(/\t/, $line);

    if(exists($name{$a[0] . ":" . $a[2] . ":" . $a[3]})){
	$techr{$a[1]} = $a[0];
	$tel{$a[1]} = $a[2];
	$ter{$a[1]} = $a[3];
	$tename{$a[1]} = $name{$a[0] . ":" . $a[2] . ":" . $a[3]};
    }
}
close A;

my $parse = $input; #parse the input file line for TE
$parse =~ s/short\.//go;
my @p = split(/\./, $parse);
my $chr = $p[0];
my $event = $p[1];
my %contig = ();
my %before = ();
my %after = ();
my $shared = 0;

open(B, "<$input");
open(C, ">>$out");
LOOP:while(my $line2 = <B>) {

    chomp $line2;
    #print $line2, "\n";
    my @b = split(/\t/, $line2);
    #print $b[0], "\n";
    my @c = split(/\./, $b[0]);
    #print $c[0], "\n";

    if(!(exists($techr{$event}))){

	last LOOP;

    }
    if(!(exists($contig{$c[4]}))) {

	$contig{$c[4]} = 1;
	$before{$c[4]} = 0;#initialize
	$after{$c[4]} = 0;#these to 0
    }    
    if($b[2] < 90) {#impose quality filter
	
	next LOOP;
    }
    if($b[1] != $techr{$event}) {#check if correct chrom

	next LOOP;
    }
    if($b[3] < 50) {#long enough segment
	
	next LOOP;
    }
#sort positions low to high
    my $start = 0;
    my $stop = 0;
    
    if($b[8] < $b[9])  {
	
	$start = $b[8];
	$stop = $b[9];
    }else{
	
	$start = $b[9];
	$stop = $b[8];
    }
    
#first check if it is shared- must cross breakpoint by 25 bp
    
    if(($start < ($tel{$event} - 25)) and ($stop > $tel{$event} +25)) {#crosses left TE breakpoint
	
	if(exists($techr{$event})) {

	    $shared = 1;
	    print C $techr{$event}, "\t", $tel{$event}, "\t", $tel{$event}, "\t", $event, ":", $tename{$event}, "\t", "shared", "\n";
	}
    }
    if(($start < ($ter{$event} - 25)) and ($stop > $ter{$event} +25)) {#crosses right TE breakpoint
	
	if(exists($techr{$event})) {

	    $shared = 1;
	    print C $techr{$event}, "\t", $tel{$event}, "\t", $tel{$event}, "\t", $event, ":", $tename{$event}, "\t", "shared", "\n";	
	}
    }   

#check if there is an absence breakpoint

    if(($stop < ($tel{$event} + 15)) and ($stop > ($tel{$event} - $dist))){
	
	$before{$c[4]} = 1;
	
    }
    if(($start > ($ter{$event} - 15)) and ($start < ($ter{$event} + $dist))){

	$after{$c[4]} = 1;

    }
}
close B;

if($shared == 0){

    while((my $key, my $value) = each(%contig)) {
	
	if(($before{$key} == 1) and ($after{$key} == 1)) {
	    
	    print C $techr{$event}, "\t", $tel{$event}, "\t", $tel{$event}, "\t", $event, ":", $tename{$event}, "\t", "absent", "\n";
	    
	}
    }
}
close C;
