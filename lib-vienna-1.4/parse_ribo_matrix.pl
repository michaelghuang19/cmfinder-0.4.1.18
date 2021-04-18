#!/usr/bin/perl -w

$single = 0;
$double = 0;
%alphabet =();
$alphabet{"Null"}=0;
$alpha_count = 1;
%pair_score = ();
while(<>){
    if (/^\s+$/){
	if (!$single) {
	    $single=1;
	    last;
	}
    }
}
exit if (!$single);

while(<>){
    if (/\s+(\S+.*\S+)\s*$/){
	@single_header = split /\s+/,$1;		
	for($i=0; $i < scalar @single_header; $i++){
	    $alphabet{$single_header[$i]}= $alpha_count;
	    $alpha_count++;
	}
	$line = 0;
	while(<>){
	    last if ($line >= scalar @single_header);
	    if (/(\S+)\s+(\S+.*\S+)\s*$/){
		$a = $1;
		@scores = split /\s+/,$2;
		for($i=0; $i < scalar @scores; $i++){
		    $b = $single_header[$i];
		    $pair_score{"$a-$b"}=$scores[$i];
		}
	    }
	    $line++;
	}
    }
    last;
}


while(<>){
    if (/^\s+$/){
	if (!$double) {
	    $double=1;
	    last;
	}
    }
}
exit if (!$double);
while(<>){
    if (/\s+(\S+.*\S+)\s*$/){
	@double_header = split /\s+/,$1;		
	for($i=0; $i < scalar @double_header; $i++){
	    $alphabet{$double_header[$i]}= $alpha_count;
	    $alpha_count++;
	}
	$line = 0;
	while(<>){
	    last if ($line >= scalar @double_header);
	    if (/(\S+)\s+(\S+.*\S+)\s*$/){
		$a = $1;
		@scores = split /\s+/,$2;
		for($i=0; $i < scalar @scores; $i++){
		    $b = $double_header[$i];
		    $pair_score{"$a-$b"}=$scores[$i];
		}
	    }
	    $line++;
	}
	last;
    }	
}    

#compute indelscore
$pair_score{"Null-Null"}=0;
foreach $a1 (sort {$alphabet{$a} <=> $alphabet{$b}} keys %alphabet){
    if ($a1 eq "Null"){
	next;
    }
    $temp_sum = 0;
    $temp_count = 0;
    $min =0;
    foreach $a2 (sort {$alphabet{$a} <=> $alphabet{$b}} keys %alphabet){
	$pair = "$a1-$a2";
	if ($alphabet{$a1} < $alphabet{$a2}){
	    $pair = "$a2-$a1";
	}
	if (exists $pair_score{$pair} && $pair_score{$pair} < 0){
	    if ($pair_score{$pair} < $min){
		$min = $pair_score{$pair};
	    }
	    $temp_sum += $pair_score{$pair};
	    $temp_count++;
	}
    }
    $score = $temp_sum/$temp_count;
    $pair_score{"$a1-Null"}=$score;    
}

print "/*  ";
foreach $a1 (sort {$alphabet{$a} <=> $alphabet{$b}} keys %alphabet){
    printf("%4s, ",$a1);
}
print "*/\n";

foreach $a1 (sort {$alphabet{$a} <=> $alphabet{$b}} keys %alphabet){
    print "{   ";
    $first=1;
    foreach $a2 (sort {$alphabet{$a} <=> $alphabet{$b}} keys %alphabet){
	if (!$first){ 
	    print ",  ";
	}
	else{
	    $first = 0;
	}
	$pair = "$a1-$a2";
	if ($alphabet{$a1} < $alphabet{$a2}){
	    $pair = "$a2-$a1";
	}
	if (exists $pair_score{$pair}){
	    printf("%.2f",$pair_score{$pair});	    
	}	
	else{
	    print "INF";
	}

    }
    print "},\n";
}


