#!perl -w

use strict;

my $input = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;

unlink(qq{$output});

my %hash = ();

open(A, "<$input");
open(B, ">>$output");

while(my $line = <A>){

    chomp $line;

    my @a = split(/\t/, $line);

    if(!(exists($hash{$a[0] . ":" . $a[1] . ":" . $a[2] . ":" . $a[4]}))) {

	$hash{$a[0] . ":" . $a[1] . ":" . $a[2] . ":" . $a[4]} = $line;

    }
}
close A;

while((my $key, my $value) = each(%hash)) {

    my @c = split(/:/, $key);

    if($c[3] eq "absent") {

	if(!(exists($hash{$c[0] . ":" . $c[1] . ":" . $c[2] . ":shared"}))) {

	    print B $value, "\n";

	}
    }else{

	print B $value, "\n";

    }
}
close B;
