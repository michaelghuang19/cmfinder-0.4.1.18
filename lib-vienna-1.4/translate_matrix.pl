#!/usr/bin/perl -w

$unpaired = 0;
$paired = 0;
$header1=0;
@col_header1=();
$header2=0;
@col_header2=();
%pairs=();

while (<>){
    if (/RIBOSUM/){
	$unpaired = 1;
	next;
    }
    if ($unpaired){
	next if (/^\s*$/ && !$header1);
	if (/^\s*$/ && $header1){
	    $unpaired = 0;
	    $paired = 1;
	    next;
	}
	@fields = split ;
	if (!$header1){
	    print "Read unpaired\n";
	    @col_header1=@fields;	
	    $header1=1;
	}
	else{
	    $syma = shift @fields;
	    for($i=0; $i< scalar @fields; $i++){
		$symb= $col_header1[$i];
		$pairs{$syma.$symb}=$fields[$i];
	    }
	    
	}
    }
    if ($paired){
	next if (/^\s*$/ && !$header2);
	@fields = split ;
	if (!$header2){
	    print "Read paired\n";
	    @col_header2=@fields;	    
	    $header2=1;
	}
	else{
	    $syma = shift @fields;
	    for($i=0; $i< scalar @fields; $i++){
		$symb= $col_header2[$i];
		$pairs{$syma.$symb}=$fields[$i];
	    }	    
	}
    }
}

@syms=("", "A",   "C",   "G",   "U",  "GC", "CG",   "AU",  "UA",  "GU",  "UG", "N");
$single_gap=;
$double_gap=; 

foreach $a (@syms){
    foreach $b (@syms){
	if ($a eq ""){
	    if (length($b)==1){
		
	    }	    
	    
	}
    }
}
