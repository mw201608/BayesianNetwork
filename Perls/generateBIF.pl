#!/usr/bin/perl

# generate a network with a collection of unconnected nodes


#node files
$in = $ARGV[0];
#network name
$networkName = $ARGV[1];
#n descrete states
$nstate = $ARGV[2];
if($nstate==0) {
    $nstate=3;
}


#open node file
open(IN, $in) || die;
$n=0;
while(<IN>) {
    if(!/^\#/) {#skip comment line 
        chop($_);
	@tmp = split(/[ \t]+/);
	$name[$n] = $tmp[0];
	$n++;
    }
}
close(IN);

$header = 
"<?xml version=\"1.0\"?>\n\n".

"<!-- DTD for the BIF format -->\n".
"<!DOCTYPE BIF [ \n".
"      <!ELEMENT BIF ( NETWORK )*> \n".
"      <!ELEMENT PROPERTY (#PCDATA)> \n".
"      <!ELEMENT TYPE (#PCDATA)> \n".
"      <!ELEMENT VALUE (#PCDATA)> \n".
"      <!ELEMENT NAME (#PCDATA)> \n".
"      <!ELEMENT NETWORK \n".
"          ( NAME, ( PROPERTY | VARIABLE | PROBABILITY |LIKELIHOOD)* )> \n".
"      <!ELEMENT VARIABLE ( NAME, TYPE, ( VALUE |  PROPERTY )* ) > \n".
"      <!ELEMENT PROBABILITY \n".
"          ( FOR | GIVEN | TABLE | ENTRY | DEFAULT | PROPERTY )* > \n".
"      <!ELEMENT PRIOR ( FOR | GIVEN | PROPERTY )* > \n".
"      <!ELEMENT FOR (#PCDATA)> \n".
"      <!ELEMENT GIVEN (#PCDATA)> \n".
"      <!ELEMENT TABLE (#PCDATA)> \n".
"      <!ELEMENT DEFAULT (TABLE)> \n".
"      <!ELEMENT ENTRY ( VALUE* , TABLE )> \n".
"]> \n".
"<BIF> \n".
"<NETWORK size=\"$n\"> \n".
"<NAME>$networkName</NAME>\n".
"<!-- Variables --> \n";

#print header
print $header;

for($i=0; $i<$n; $i++) {
    if($nstate<0) {
	print 
	    "<VARIABLE> \n",
	    "    <NAME>", $name[$i], "</NAME>\n",
	    "    <TYPE>continuous</TYPE>\n",
	    "</VARIABLE>\n";   
    }
    elsif($nstate==3) {
   	print 
	    "<VARIABLE> \n",
	    "    <NAME>", $name[$i], "</NAME>\n",
	    "    <TYPE>discrete</TYPE>\n",
	    "    <VALUE>down</VALUE>\n",
	    "    <VALUE>no</VALUE>\n",
	    "    <VALUE>up</VALUE>\n",
	    "</VARIABLE>\n";    
    }
    else {
	print 
	    "<VARIABLE> \n",
	    "    <NAME>", $name[$i], "</NAME>\n",
	    "    <TYPE>discrete</TYPE>\n";
        for($j=0; $j<$nstate; $j++) {  
	    print
	    "    <VALUE>S$j</VALUE>\n";
	}
	print
	    "</VARIABLE>\n";    
    }
}

#print footer
print 
"</NETWORK>\n",
"</BIF>\n";
