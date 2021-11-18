#! /usr/bin/perl

#update banned matrix based on prior information and selected cutoff

$data = $ARGV[0];
$banFile = $ARGV[1];
$cisFile = $ARGV[2];

#read data
open(IN, $data) || die;
$n=0;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    $name{$tmp[0]}= $n;
    $n++;
}
close(IN);

#initialization
for($i=0; $i<$n; $i++) {
    for ($j=0; $j<$n; $j++) {
        $ban[$i][$j] = 0;
    }
    $ban[$i][$i] = 1;
}

#read banned file
open(IN, $banFile) || die;
$i=0;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    for ($j=0; $j<@tmp; $j++) {
	$ban[$i][$j] = int($tmp[$j]);
    }
    $i++;
}
close(IN);

#initialize cis vector
for ($i=0; $i<$n; $i++) {
    $cis[$i] = 0;
}

#read cis
open(IN, $cisFile) || die;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $cis[$name{$tmp[0]}] =1;
    $chrom[$name{$tmp[0]}] = $tmp[1];
    $pos[$name{$tmp[0]}] = $tmp[2];
}
close(IN);

#update banned matrix based on cis
for ($i=0; $i<$n; $i++) {
    for ($j=0; $j<$n; $j++) {
	if ($cis[$j] ==1 && $cis[$i] ==0) {
	    $ban[$i][$j] = 1;
	}
	elsif($cis[$j] ==1 && $cis[$i] ==1 && 
	      ($chrom[$i] != $chrom[$j] || abs($pos[$i] -$pos[$j])>30000)) {
	    $ban[$i][$j] = 1;
	}
    }
}

#with high quality CIS
if (@ARGV>3) {
    print "updating based on cis_HQ...\n";
    #read cis_HQ
    open(IN, $ARGV[3]) || die;
    while(<IN>) {
	chop($_);
	@tmp = split(/[ \t]+/);
	$j= $name{$tmp[0]};
	for ($i=0; $i<$n; $i++) {
	    $ban[$i][$j] = 1;
	}
    }
    close(IN);   
}

#output ban matrix
print "Updating matrix ....\n";
open(OUT, ">$banFile") || die;
for ($i=0; $i<$n;$i++) {
    for ($j=0; $j<$n; $j++) {
	printf(OUT "%d\t", $ban[$i][$j]);
    }
    print OUT "\n";
}
close(OUT);
