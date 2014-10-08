#!perl -w

use strict;
my $input = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
#find the annotated breakpoint, identify which contig that is and then find the TE contig that matches and select the best TE

my %te = ();
my %bestlen = ();
my %beste = ();
my %contig = ();
open(B, "<$input");
while(my $line2 = <B>) {

    chomp $line2;
    my @b = split(/\t/, $line2);

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
close B;
open(A, ">>$output");
while((my $key, my $value) = each(%te)) {

    print A $key, "\t", $value, "\t", $bestlen{$key}, "\n";
    
}
close A;
