@bases = ("A","C","G", "U");
$header="";
for($i=0; $i<4; $i++){
    for($j=0; $j < 4;$j++){
	$header= $header." &".$bases[$i].$bases[$j];	
    }
}


while(<>){
    if (/Frequency/){
	print "Freq $header\\\\\n";
	$l = <>;
	$l =~ s/ +/&/g;
	print $l, "\n";
    }
    if (/Rate Matrix/){
	for($i=0; $i<4; $i++){
	    for($j=0; $j < 4;$j++){	
		$l = <>;
		$l =~ s/ +/&/g;
		$h = $bases[$i].$bases[$j];
		print "$h & $l\\\\\n";
	    }
	}
	
    }
}
