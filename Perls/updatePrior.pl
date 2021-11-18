#!/usr/bin/perl

#update prior in a target prior file given prior from a source file
$data = $ARGV[0];
$in0 = $ARGV[1]; #source prior
$in1 = $ARGV[2]; #prior file to be updated
$offset = $ARGV[3];
$scale = $ARGV[4];

open(IN, $data) || die;
$n=0;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    $names{$tmp[0]} = $n;
    $n++;
}
close(IN);
print "nNode $n\n";

#output
$out = "A".rand();
unlink("$out");
$delta =10000;
for ($I=0; $I<int(($n+1)/$delta+1); $I++) {
    #initializing
    print "initializing...\n";
    my @v = ();
    for ($i=$I*$delta; $i<($I+1)*$delta & $i<$n; $i++) {
		#print "\015$i";
		for ($j=0; $j<$n; $j++) {
			$v[$i][$j] = 0;
		}
    }
    print "initialization is done!\n";

    $k=1;
    #read in priors
    open(IN0, $in0) || die;
    while(<IN0>) {
		chop($_);
		@tmp = split(/[ \t]+/);
		$i= $names{$tmp[0]};
		$j = $names{$tmp[2]};
		if ($tmp[4] ==0) {
			$tmp[4] = 0.01;
		}
		if ($i>=($I*$delta) && $i<($I+1)*$delta) {
			#print "$i $j\n";
			$v[$i][$j] = log(abs($tmp[4])/$offset)*log($scale)*$k;
		}
    }
    close(IN1);

    print "outputing...\n";
    open(OUT, ">>$out") || die;   
    #read in priors
    open(IN1, $in1) || die;
    while(<IN1>) {
		chop($_);
		@tmp = split(/[ \t]+/);
		$i= $names{$tmp[0]};
		$j = $names{$tmp[2]};
		if ($i>=($I*$delta) && $i<($I+1)*$delta) {
			printf(OUT "%s %s %s %f %f\n", $tmp[0], $tmp[1], $tmp[2], $v[$i][$j], $tmp[4]);
		}
	}
    close(IN1);
    close(OUT);
}
rename($out, $in1);
