#!/usr/bin/perl

#update prior by adding up prior evidence
$data = $ARGV[0];
$list = $ARGV[1];
$prior = $ARGV[2];

#read data
open(IN, $data) || die;
$n=0;
$ndim = 0;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    if($ndim==0) {
		$ndim = @tmp;
    }
    $name{$tmp[0]}= $n;
    $n++;
}
close(IN);
$scale = 2*log($ndim);

#read in list
open(IN, $list) || die;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $key = "$tmp[0] $tmp[1]";
    if(@tmp<3) {
		$V{$key} = $scale;
    }
    else {
		if(length($V{$key})==0) {
			$V{$key} = $tmp[2];
		}
		else {
			if ($V{$key}<$tmp[2])  {
				$V{$key} = $tmp[2];
			}
		}
    }
}
close(IN);

#read in priors
open(IN, $prior) || die;
$out = 'A'.rand();
open(OUT, ">$out") || die;
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $key = "$tmp[0] $tmp[2]";
    
    if (length($V{$key}) >0)  {
		$v = $tmp[3]+ $V{$key};
		#print "$key $tmp[3]  $v\n";
    }
    else {
		$v = $tmp[3];
    }

    print OUT "$tmp[0] $tmp[1] $tmp[2] $v $tmp[4]\n";
}
close(IN);
close(OUT);

rename($out, $prior);
