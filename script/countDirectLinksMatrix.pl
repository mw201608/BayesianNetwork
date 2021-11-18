#!/usr/bin/perl -w
use strict;
# count links in a series of graphs

#the node file
my $in = $ARGV[0];

#the pattern of the series of graphs
my $pat = $ARGV[1]; #result file prefix
my $start = $ARGV[2];
my $end =  $ARGV[3];
my $output1 = $ARGV[4];
my $output2 = $ARGV[5];
my $cutoff = $ARGV[6]; #minimal proportion of networks that the links are present
my %node;
my @Node;

srand(1000);
#get list of nodes
my$seedF = 1;
my $n=0;
open(IN, $in) || die;
while(<IN>) {
	my @tmp = split(/[ \t]+/,substr($_,0,256));
	$node{$tmp[0]}=$n;
	$Node[$n] = $tmp[0];
	$n++;
}
close(IN);

#initial count
my @count;
for(my $i=0; $i<$n; $i++) {
	for(my $j=0; $j<$n; $j++) {
		$count[$i][$j] = 0;
	}
}
print "initiation is done....\n";

#check the series of graphs
my $total = 0;
for(my $i= $start; $i<=$end; $i++) {
	my $file = $pat.'.'.$i;
	if( -e "$file" ) {
		$total++;
		open(FILE, $file) || die "$file is not found.";
		while(<FILE>) {
			if(/->/) {
				chop($_);
				my @tmp=split('->');
				if(length($node{$tmp[0]}) ==0 || length($node{$tmp[1]}) ==0 ) {
					print "$_\n";
					exit(0);
				}
				$count[$node{$tmp[0]}][$node{$tmp[1]}] ++;
			}
		}
		close(FILE);
	}
}
#clear
%node =();
#output
if($total==0){
	die "Zero BN file $pat.* found\n";
}
print "outputing links...\n";
my $hcutoff = $cutoff/2;
open(OUT1, ">$output1") || die;
print OUT1 "digraph G {\n";
for(my $i=0; $i<$n; $i++) {
	for(my $j=0; $j<$n; $j++) {
		my $r = $count[$i][$j] /$total;
		my $R = $count[$j][$i] /$total;
		if($cutoff>0) {
			if($r>=$hcutoff && ($R+$r)>=$cutoff && $r>=$R) {
				if($r==$R) {
					$r += (rand(100)/10000);
				}
				printf(OUT1  "%s->%s [label=%f];\n", $Node[$i],$Node[$j],$r);
			}
		}
	}
}
print OUT1 "}\n";
close(OUT1);

print "outputing link matrix...\n";
open(OUT2, ">$output2") || die;
for(my $i=0; $i<$n; $i++) {
	for(my $j=0; $j<$n; $j++) {
		my $r = $count[$i][$j] /$total;
		my $R = $count[$j][$i] /$total;
		if($cutoff>0) {
			if($r>=$hcutoff && ($R+$r)>=$cutoff && $r>=$R) {
				if($r==$R) {
					$r += (rand(100)/10000);
				}
				printf( OUT2 "%f\t", $r);
			}
			else {
				print OUT2 "0\t";
			}
		}
		else {
			printf( OUT2 "%f\t", $r);
		}
	}
	print OUT2 "\n";
}
close(OUT2);

print STDERR "Total $total\n";
