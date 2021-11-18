#!/usr/bin/perl

#get a histogram from a set of data

#min value
$min = $ARGV[0];

#max value
$max  = $ARGV[1];

#increment
$incr = $ARGV[2];

#input file
$file = $ARGV[3];

#the column to scan
$col = $ARGV[4];

$absF = $ARGV[5];

#calculate dimension according to min, max and increment
$N = int(($max-$min)/$incr)+1;

#initialize count
for($i=0; $i<$N; $i++) {
    $count[$i] = 0;
}
$sum =0;

#scan file
open(IN, $file) || die;
while(<IN>) {
    if(!/^\#/) {
		chop($_);
		@tmp = split(/[ \t]+/);
		if($absF) {
			$value = abs($tmp[$col]);
		}
		else{
			$value = $tmp[$col];
		}

		if($value<$min) {
			#print "min $value\n";
			$n = 0;
		}
		elsif($value>$max) {
			#print "max $value\n";
			$n = $N-1;
		}
		else {
			$n = int(($value-$min)/$incr);
		}
		$count[$n]++;
		$sum++;
    }
}
close(IN);

#output
$s=0;
for($i=0; $i<$N; $i++) {
    $l = $min+$i*$incr;
    $v = $count[$i]/$sum;
    $s+= $v;
    print "$l\t$v\t$s\t$count[$i]\n";
}
