#!/usr/bin/perl

#update banned matrix based on prior information and selected cutoff
#if a node has too many candidate parents, then only the top ones are selected.


$data = $ARGV[0];
$priorFile = $ARGV[1];
$ncol = $ARGV[2];
$cutoff1 = $ARGV[3];
$cutoff2 = $ARGV[4];
$perCutoff = $ARGV[5];  #percentage cutoff
$banFile = $ARGV[6];

#read data
open(IN, $data) || die;
$n=0;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    $name{$tmp[0]}= $n;
    $n++;
}
close(IN);

#maximum of candidate parents
$nMax = int($n*$perCutoff);

#initialization
for($i=0; $i<$n; $i++) {
    for ($j=0; $j<$n; $j++) {
		$v[$i][$j] = 0;
    }
    $v[$i][$i] = 1;
    $c[$i] =0;  #node-specific cutoff 
}

#read prior file
open(IN, $priorFile) || die;
$count=0;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $i=$name{$tmp[0]};
    $j=$name{$tmp[2]};
    if($tmp[$ncol] >=$cutoff1) {
		$v[$i][$j] = 1.0 ;
    }
    $count++;
    if($count%100000) {
		print "\015$count";
    }
}
print "\n";
close(IN);

#get node specific cutoff
for ($j=0; $j<$n; $j++) {
    #print "processing node $j\n";
    $s=0;
    for ($i=0; $i<$n;$i++) {
		$t = $v[$i][$j];
		$s += $t;
    }
    if($s>$nMax) {
		$c[$j] = $cutoff2;
    }
    else {
		$c[$j] = $cutoff1;
    }
}
    
#reset the matrix
for($i=0; $i<$n; $i++) {
    for ($j=0; $j<$n; $j++) {
		$v[$i][$j] = 0;
    }
    $v[$i][$i] = 1;
}

#rescan the prior
open(IN, $priorFile) || die;
$count=0;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $i=$name{$tmp[0]};
    $j=$name{$tmp[2]};
    if($tmp[$ncol] <$c[$j]) {
		$v[$i][$j] = 1 ;
    }
    $count++;
    if($count%100000) {
		print "\015$count";
    }
}
print "\n";
close(IN);

#read in old ban matrix
if (-s $banFile >0) {
    print "reading old ban matrix....\n";
    open(IN, $banFile) || die;
    $i=0;
    while(<IN>) {
		chop($_);
		@tmp = split(/[ \t]+/);
		for ($j=0; $j<$n; $j++) {
			if( $tmp[$j]==1) {
				$v[$i][$j] = 1;
			}
		}
		$i++;
    }
    close(IN);
}

#output ban matrix
open(OUT, ">$banFile") || die "$banFile can not be created";
for ($i=0; $i<$n;$i++) {
    for ($j=0; $j<$n; $j++) {
		printf(OUT "%d\t", $v[$i][$j]);
    }
    print OUT "\n";
}
close(OUT);
