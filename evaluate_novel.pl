#!perl -w
use strict;
my $inline = shift(@ARGV) or die;
my $pos = shift(@ARGV) or die;
my $out = shift(@ARGV) or die;
unlink(qq{$out});
my %ch = ();
my %start = ();
my %stop = ();
open(A, "<$pos");
#storing the positions of all novel breakpoints
while(my $line = <A>) {#have to sort this way because I sorted this way when assigning left/rightness to fasta files
    chomp $line;
    my @a = split(/\t/, $line);
    $ch{$a[2]} = $a[0];
    
    if($a[3] <= $a[4]) {
	$start{$a[2]} = $a[3];
	$stop{$a[2]} = $a[4];
    }else{
	$start{$a[2]} = $a[4];
	$stop{$a[2]} = $a[3];
    }
}
close A;

opendir DIR, ".";
my %lhash = ();#store id number and if it is reconstructed or not
my %rhash = ();
my @all = grep {/\d+\.\d+\.left\.fasta\.contigs\.chromosome\.output/} readdir DIR;
closedir DIR;
open(C, ">>$out");
my $start = 0;
my $stop = 0;
foreach my $all (@all) {
    open(B, "<$all");
    my $found = 0;
    my $current = 0;
  LOOP:while(my $line2 = <B>) {
      chomp $line2;
      my @b = split(/\t/, $line2);
      my @c = split(/\./, $b[0]);
      $current = $c[1]; #the id of the te
      if($b[1] != $ch{$c[1]}){
	  next LOOP;
      }
      if($b[8] < $b[9]) {
	  $start = $b[8];
	  $stop = $b[9];
      }else{
	  $start = $b[9];
	  $stop = $b[8];
      }
      if((($start{$c[1]} - $start) > 25) and (($stop - $start{$c[1]}) > 25)) { 
	  
	  $lhash{$current} = "reconstructed";	
	  $found = 1;
	  last LOOP; 
      }
  }
    if($found == 0) {
#did not reconstruct - is it a novel TE? - only concerned with complete bp - now call annotate_TEs.pl
	my $dir = "/genomics/julie/founder_bwa_te/pull_reads/line" .  $inline . "_verify/";
	my $contig = $all;
	$contig =~ s/\.chromosome\.output//go;
	$contig = $dir . $contig;
	print $contig, "\n";
	my $posest = "/genomics/julie/founder_bwa_te/pull_reads/combined_founders_novel_pos";
	my $chnum = "/genomics/julie/founder_bwa_te/chrom_num";
	my $o1 = "/genomics/julie/founder_bwa_te/pull_reads/line" . $inline . "_new_novel";
	my $o2 = "/genomics/julie/founder_bwa_te/pull_reads/line" . $inline . "_new_novel_bed";

	system(qq{perl /home/julie/pipeline/annotate_TEs_novel.pl $dir $contig $posest $chnum $o1 $o2});
	$lhash{$current} = "different";	
	$found = 0; 
    }
    close B;
}

opendir DIR2, ".";
my @all2 = grep {/\d+\.\d+\.right\.fasta\.contigs\.chromosome\.output/} readdir DIR2;
#my @all = grep {/13.3381.left.fasta.contigs.chromosome.output/} readdir DIR;
closedir DIR2;
foreach my $all2 (@all2) {
    open(D, "<$all2");
    my $found2 = 0;
    my $current2 = 0;
  LOOP2:while(my $line3 = <D>) {
      chomp $line3;
      my @e = split(/\t/, $line3);
      my @f = split(/\./, $e[0]);
      $current2 = $f[1];
      if($e[1] != $ch{$f[1]}){
	  next LOOP2;
      }
      if($e[8] < $e[9]) {
	  $start = $e[8];
	  $stop = $e[9];
      }else{
	  $start = $e[9];
	  $stop = $e[8];
      }
      if((($stop{$f[1]} - $start) > 25) and (($stop - $stop{$f[1]}) > 25)) { 
	  
	  $rhash{$current2} = "reconstructed";	
	  $found2 = 1;
	  last LOOP2; 
      }
  }
    if($found2 == 0) {

	$rhash{$current2} = "different";
#	print $current2, "\n";	
	$found2 = 0; 
    }
    close D;
}
while((my $key, my $value) = each(%ch)){
    
    #print $key, "\n";
    #print $rhash{$key}, "\n";
    #print $lhash{$key}, "\n";

    if(!(exists($rhash{$key}))) {

	$rhash{$key} = "void";
    }
    if(!(exists($lhash{$key}))) {
	$lhash{$key} = "void";
    }
    
    if(($rhash{$key} eq "reconstructed") and ($lhash{$key} eq "reconstructed")) {
	print C $ch{$key}, "\t", $key, "\t", $rhash{$key}, "\t", $lhash{$key}, "\n";
	#print C $key, "\t both \n";
    }elsif(($rhash{$key} eq "reconstructed") and ($lhash{$key} ne "reconstructed")) {
	print C $ch{$key}, "\t", $key, "\t", $rhash{$key}, "\t", $lhash{$key}, "\n";
	#print C $key, "\t right \n";
    }elsif(($rhash{$key} ne "reconstructed") and ($lhash{$key} eq "reconstructed")) {
	print C $ch{$key}, "\t", $key, "\t", $rhash{$key}, "\t", $lhash{$key}, "\n";
	#print C $key, "\t left \n";
    }elsif(($rhash{$key} eq "different") and ($lhash{$key} eq "different")) {
	print C $ch{$key}, "\t", $key, "\t", $rhash{$key}, "\t", $lhash{$key}, "\n";
	#print C $key, "\t neither \n";
    }else{
# }elsif(($rhash{$key} eq "void") and ($lhash{$key} eq "void")) {
	print C $ch{$key}, "\t", $key, "\t", $rhash{$key}, "\t", $lhash{$key}, "\n";
	#print C $key, "\t void \n";
    }
    
}
close C;
