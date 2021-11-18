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
$BNBIN = "/hpc/users/wangm05/bin";
$BN_PERL = "/hpc/users/wangm05/work/bayesian/RIMBANet/Perls";

#r of pvalue=0.01
$corrFile = "$BN_PERL/corrP01.txt";
open(IN, $corrFile) || die ("$corrFile is not found.");
while(<IN>) {
    chop($_);
    @tmp = split(/[ \t]+/);
    $Offset[$tmp[0]] = $tmp[1];
}
close(IN);

$prog = "$BNBIN/testBN -f 0 ";
if(@ARGV<3) {
    print "$0 step BD BC\n";exit(0);
}
$step = $ARGV[0];
$BD = $ARGV[1];
$BC = $ARGV[2];

$data = "data.txt";
$init = "init.xml";
$prior = "qtl.prior";
$ban = "banned.matrix";
$cis = "cisqtl.txt";
$cis_HQ = "cisqtl_HQ.txt";
$causal ="bannedList_causal.txt";
$partialCorr =  "partialCorr_list.txt";
$mi_cutoff = "mi.cutoff1";
$qtloverlap = "QTLOverlap.txt";
$cisTransOverlap = "cisTransOverlap.txt";
$fixedPriors ="fixed.priors";
$fixedPriors_TF_PPI ="fixed.TF_PPI";

#get data dimension
open(IN, "$BC/$data") || die;
$nNode=0;
while(<IN>){
    chop($_);
    @tmp = split(/[ \t]+/);
    $nSample = @tmp-1;
    $nNode++;
}
close(IN);

#step 0: create ban matrix
#step 1: initialize priors
#step 2: update priors
#step 3: add fixed links
#step 4: calculate structure
if ($step==0) {
    #prepare for continuous data
    if(-s "$BC/$prior"==0) {
		$cmd = "cd $BC; perl $BN_PERL/generateBIF.pl $data \"continuous \" -1 >$init";
		print "$cmd\n";
		system($cmd);
		$cmd = "cd $BC; $prog -b $init -d $data -t 0 -T 0 -D $nSample -o junk -r 1 -L 1  >junk ; egrep  \\> junk > $prior";
		print "$cmd\n";
		system($cmd);
		unlink("$BC/junk");
    }

    #prepare for discrete data data
    $cmd = "cd $BD; perl $BN_PERL/generateBIF.pl $data \"discrete \"  >$init";
    print "$cmd\n";
    system($cmd);

    if(-s "$BD/$prior"==0) {
		$cmd = "cd $BD; $prog -b $init -d $data -t 0 -T 0 -D $nSample -o junk -r 1 -L 1  >junk ; egrep  \\> junk > $prior";
		print "$cmd\n";
		system($cmd);
		unlink("$BD/junk");
    }

    #goto L1;
    #determined cutoff value at X percentile
    $mi_cutoffR =0.8; #mi cutoff rate
    if(!-e "$BD/$mi_cutoff" || -s "$BD/$mi_cutoff" ==0) {
		print "calculating mi cutoff...\n";
		@res = `cd $BD; perl $BN_PERL/getHistogram.pl 0 0.3 0.001 $prior 4 | egrep -v "max|min" `;
		$ti = 0;
		$p=0;
		for ($i=0; $i<@res; $i++) {
			@tmp = split(/[ \t]+/, $res[$i]);
			if ($p<$mi_cutoffR && $tmp[2]>$mi_cutoffR) {
				break;
			}
			else {
				$ti = $tmp[0];
				$p = $tmp[2];
			}
		}
		#upper bound
		$ti2 = 0;
		$p=0;
		$mi_cutoffR2 =0.8; #mi cutoff rate
		for ($i=0; $i<@res; $i++) {
			@tmp = split(/[ \t]+/, $res[$i]);
			if ($p<$mi_cutoffR2 && $tmp[2]>$mi_cutoffR2) {
				break;
			}
			else {
				$ti2 = $tmp[0];
				$p = $tmp[2];
			}
		}
    }
    else {#cutoff based on permuation
        open(IN, "$BD/$mi_cutoff") || die;
		$_=<IN>;
		chop($_);
		$ti = $_;
		close(IN);
	}

    print "ti=$ti\n";
    print "ti2=$ti2\n";
    
    #create the banned matrix
    $r= (1-$mi_cutoffR)*1.2; #percentage of nodes that could be parent for a node (comment added by Minghui Wang)
    $cmd = "cd $BD; perl $BN_PERL/createBannedMatrix.pl $data $prior 4 $ti $ti2 $r $ban";
    print "$cmd\n";
    system($cmd);

    #update banned matrix based on cis-QTL
    $cmd = "cd $BD; perl $BN_PERL/updateBannedMatrix4Cis.pl $data $ban $cis $cis_HQ";
    print "$cmd\n";
    system($cmd);
 
    #update banned matrix based on causal-list
    if (-e "$BD/$causal") {
		$cmd = "cd $BD; perl $BN_PERL/updateBannedMatrixFromList.pl $data $ban $causal ";
		print "$cmd\n";
		system($cmd);
    }
}   
if ($step ==1) {
    #update prior based on continuous data
    $offset = $Offset[$nSample];
    if ($offset ==0) {
		print "Error: offset is 0 for $nSample!\n";
		exit(0);
    }

    $cmd="cd $BD; perl $BN_PERL/updatePrior.pl $data ../$BC/$prior $prior $offset $nNode";
    print "$cmd\n";
    system($cmd);

    #update based on partial correlation
    if (-e "$partialCorr") {
		$cmd="cd $BD; perl $BN_PERL/partialCorr2Prior.pl $data $prior $partialCorr $offset $nNode";
		print "$cmd\n";
		system($cmd);
    }
}

if($step ==2) {
    #scan prior base on causality
    if (-e "$BD/prior.causality") {
        $cmd="cd $BD; perl $BN_PERL/addPrior.pl $data prior.causality $prior ";
        print "$cmd\n";
        system($cmd);
    }
    
    #scan prior base on KEGG
    if (-e "$BD/prior.KEGG") {
        $cmd="cd $BD; perl $BN_PERL/addPrior.pl $data prior.KEGG $prior ";
        print "$cmd\n";
        system($cmd);
    }


    #scan prior base on PPI
    if (-e "$BD/prior.PPI") {
		$cmd="cd $BD; perl $BN_PERL/addPrior.pl $data prior.PPI $prior ";
		print "$cmd\n";
		system($cmd);
    }
 
    #scan prior base on TF_PPI
    if (-e "$BD/prior.TF_PPI") {
		$cmd="cd $BD; perl $BN_PERL/addPrior.pl $data prior.TF_PPI $prior ";
		print "$cmd\n";
		system($cmd);
    }
}


if($step==3) {
    #add fixed priors
    if(-s "$BD/$fixedPriors" >0) {
		$cmd = "cd $BD; perl $BN_PERL/addPriorFromFile.pl $fixedPriors $init >junk ";
		print "$cmd\n";
		system($cmd);
    }
}

if($step==4) {
    # calculate structure
    $cmd="cd $BD; perl $BN_PERL/run_w_genetics.pl $nSample $nNode $init $data $ban $prior &";
    print "$cmd\n";
    system($cmd);
}
