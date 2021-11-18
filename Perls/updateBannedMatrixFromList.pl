#! /usr/bin/perl

#update banned matrix based on prior information and selected cutoff

$data = $ARGV[0];
$banFile = $ARGV[1];
$listFile = $ARGV[2];

#read data
open(IN, $data) || die;
$n=0;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    $name{$tmp[0]}= $n;
    $n++;
}
close(IN);

#read old banned matrix
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

#read list file
open(IN, $listFile) || die;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $i=$name{$tmp[0]};
    $j=$name{$tmp[1]};
    $ban[$i][$j] = 1;
}
close(IN);

#output ban matrix
open(OUT, ">$banFile") || die;
for ($i=0; $i<$n;$i++) {
    for ($j=0; $j<$n; $j++) {
	printf(OUT "%d\t", $ban[$i][$j]);
    }
    print OUT "\n";
}
close(OUT);
