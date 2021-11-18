#!/usr/bin/perl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A package for Bayesian network reconstruction
#
# Author:  Dr. Jun Zhu, Amgen Inc., Thousand Oaks, CA, 2002
#          Dr. Jun Zhu, Rosetta Informatics, a wholly owned subsidiary of Merck & CO., Seattle, WA, 2004, 2008
#
# Acknowledge:  Thanks to Amgen Inc. and Merck & CO. for their generous supports
#
# If you use this package, please cite the following references
#
# (1) Zhu J, Lum PY, Lamb J, GuhaThakurta D, Edwards SW, et al. An integrative genomics approach to the
#     reconstruction of gene networks in segregating populations. Cytogenet Genome Res 105: 363-374 (2004)
# (2) Zhu J, Wiener MC, Zhang C, Fridman A, Minch E, Lum PY, Sachs JR, & Schadt EE Increasing the power to
#     detect causal associations by combining genotypic and expression data in segregating populations
#     PLoS Comput Biol 3, e69. (2007)
# (3) Zhu, J., Zhang, B., Smith, E.N., Drees, B., Brem, R.B., Kruglyak, L., Bumgarner, R.E. and Schadt, E.E.
#     Integrating large-scale functional genomic data to dissect the complexity of yeast regulatory networks,
#     Nat Genet, 40, 854-861 (2008)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

$BNBIN='/hpc/users/wangm05/bin';
$BN_PERL = "/hpc/users/wangm05/work/bayesian/RIMBANet/Perls";

$workDir = $ENV{PWD};

$dim = $ARGV[0];
$nNodes = $ARGV[1];
$init = $ARGV[2];
$data = $ARGV[3];
$ban = $ARGV[4];
$prior = $ARGV[5];
if(@ARGV>6) {
    $sleeptime = $ARGV[6];
}
else {
   $sleeptime = 10;
}

if(@ARGV>7) {
    $label = $ARGV[7];
}
else {
    $label = "Label";
}


$NS = 1;

$interactive =1;
$prog = "$BNBIN/testBN -f 0 -M 5000000 ";

#the values for new ban matrix
$qratio = 1/($nNodes+1000);
$te=0;
$N=1000;
$alpha = 0.65-($dim/100)*0.015;
$r=1;

#calculate trylist
system("touch Trylist");
&submit2Que(0, $N);
#check whether all results have finished.
$n=0;
while($n<$N) {
    sleep(30);

    #count try lists to be concatenated
    $n=0;
    for($i=0; $i<$N; $i++) {
	if(-s "trylist.$i" >0) {
	    $n++;
	}
    }

    if($n>2) {
	$cmd="perl $BN_PERL/updateTrylist.pl";
	print "$cmd\n";
	system($cmd);
    }

    #count prematured terminated results
    $n=0;
    $time = time();
    for($i=0; $i<$N; $i++) {
	if(!-e "result.out.$i" || (-s "result.out.$i"==0 && -e "junkK.$i")) {
	    $file = "junkK.$i";
	    $mtime = (stat($file))[9];
	    $diff = ($time-$mtime)/60; 
	    if ($diff>180) {
		system("ls -l $file");
		unlink("result.out.$i");
		unlink("junkK.$i");
		$n++;
	    }
	}
    }
    if($n>0) {
	&submit2Que(0, $N);
    }


    #count finished results
    $n=0;
    for($i=0; $i<$N; $i++) {
        if(-s "result.out.$i" >0) {
            $n++;
        }
    }
    system("pwd");
    print "waiting for $n/$N\n";
}

#create structure
$cmd = "perl $BN_PERL/createStructure.pl $dim $nNodes $init $data \"$label\" ";
print "$cmd\n";
system($cmd);

sub submit2Que  {
    $start = $_[0];
    $end = $_[1];
    
    $time =time();

    my($njobs)=0;
    my($M);
    for($j=$start; $j<$end; $j+=$NS) {
	$cmd = "";
	for ($i=$j; $i<$j+$NS; $i++) {
	    if( !-e "result.out.$i") {
		system("touch result.out.$i");
		$time++;
		if(-s "Trylist" ==0) {
		    $up ="-U trylist.$i";
		}
		else {
		    $up = "";
		}
		$up = "";
		$cmd .= "cd $workDir; $prog -s $time -b $init -d $data -t -1 -T $te -D $dim -r $r -P $prior -a $alpha -q $qratio -l Trylist -g $ban $up -o result.out.$i>junkK.$i ;";

	    }
	}
	if(length($cmd) >0) {
	    print "$cmd\n";
	    $njobs++;
	
	    if($interactive) {
		system($cmd);
	    }
	    else {
		if($njobs<600) {  
		    $q = "short64";
		}
		else {
		    $q = "long64";
		}
		$q="short64";
		system("bsub -o out.$i -e err.$i   -R  \"rusage[mem=1950]\" \"$cmd \" ");
	    }
	    if ($njobs%100 ==0) {
		sleep($sleeptime);
	    }
	}
    }
}
