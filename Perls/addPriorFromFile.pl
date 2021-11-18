#!/usr/bin/perl

# add prior information to a network file

#prior files
$in = $ARGV[0];
#initial network 
$network = $ARGV[1];
$output = "A".rand();

#open initial network file and output old information
open(IN, $network) || die;
open(OUT, ">$output") || die;
while(<IN>) {
    if(!/\/NETWORK/ && !/\/BIF/) {
	print OUT $_;
    }
}
close(IN);

#append new prior information to the end

#open prior information file
#format is : parent child prior
open(IN, $in) || die;
while(<IN>) {
    if(!/^\#/) {#skip comment line 
        chop($_);
	@tmp = split(/[ \t]+/);
	$pname = $tmp[0];
	$cname = $tmp[1];
	$prior = $tmp[2];
	
	print OUT
	    "<PRIOR> \n",
	    "    <FOR>", $cname, "</FOR>\n",
	    "    <GIVEN>", $pname,"</GIVEN>\n",
	    "    <PROPERTY>",$prior,"</PROPERTY>\n",
	    "</PRIOR>\n";    
    }
}

#print footer
print OUT
"</NETWORK>\n",
"</BIF>\n";

close(IN);
close(OUT);
rename($output, $network);
