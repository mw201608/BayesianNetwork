#!/usr/bin/perl

#update prior
$data = $ARGV[0];
$p1 = $ARGV[1];
$p2 = $ARGV[2];
$offset = $ARGV[3];
$scale = $ARGV[4];

#read data
open(IN, $data) || die;
$n=0;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    $name{$tmp[0]}= $n;
    $n++;
}
close(IN);

#partial correlation result
open(P2, $p2) || die;
while(<P2>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $key = "$tmp[0] $tmp[2]";
    $V{$key} = $tmp[4];
}
close(P2);

#read in priors
open(P1, $p1) || die;
$out = 'A'.rand();
open(OUT, ">$out") || die;
while(<P1>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $key = "$tmp[0] $tmp[2]";
    if ($V{$key} >0) {
		$t = $V{$key};
		$v = log($t/$offset)*log($scale);
    }
    else {
		$v = $tmp[3];
    }

    print OUT "$tmp[0] $tmp[1] $tmp[2] $v $tmp[4]\n";
}
close(P1);
close(OUT);

rename($out, $p1);
