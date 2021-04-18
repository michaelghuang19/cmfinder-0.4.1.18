use Class::Struct;

use Getopt::Long;
$path= $ENV{CMfinder};
$bin_path = "$path/bin";
do "$bin_path/io.pl";


$dump_file = shift @ARGV;
@files = @ARGV;


sub match_file{  
    my @files = ();
    $f = shift @_;
    while(<$f*>) {
	push @files, $_;
    }
    return @files;
}

sub better{ # to be related to regression model
    my($stat1, $stat2) = @_;
    return 1 if ($stat1->{"Energy"} < 0 && $stat2->{"Energy"} > 0);
    return 0 if ($stat1->{"Energy"} > 0 && $stat2->{"Energy"} < 0);
    
    return 1 if ($stat1->{"Weight"} >= $stat2->{"Weight"} && $stat1->{"Score"} - $stat2->{"Score"} > 20);
    return 0 if ($stat1->{"Weight"} <= $stat2->{"Weight"} && $stat2->{"Score"} - $stat1->{"Score"} > 20);

    return 1 if ($stat1->{"Weight"} >= $stat2->{"Weight"} && $stat1->{"BP"} > $stat2->{"BP"});
    return 0 if ($stat1->{"Weight"} <= $stat2->{"Weight"} && $stat1->{"BP"} < $stat2->{"BP"});
    
    return ($stat1->{Weight} * $stat1->{"Score"} > $stat2->{"Weight"} * $stat2->{"Score"});
}

@align_files = ();
foreach $f (@files) {
    push @align_files, match_file($f);
}



my $i;

%alignments=();
%all_stat = ();
%all_score=();
%to_remove=();

foreach $f (@align_files) {
    # print STDERR $f, "\n";
    $alignments{$f} = read_stockholm($f);
    $result = `$path/summarize $f`;
    # print "$result\n";
    my %stat = $result=~ /\s*(\S+)=(\S+)\s*/g;
    $all_stat{$f}=\%stat;
    $all_score{$f}=0;
    $all_score{$f} = $stat{"Weight"} * $stat{"Score"};
}



while(scalar @align_files > 0) 
{    
    $f1 = shift @align_files;    
    next if (exists $to_remove{$f1}) && $to_remove{$f1};
    for($j = 0; $j < scalar @align_files; $j++) {
	$f2 = $align_files[$j];
	next if (exists $to_remove{$f2}) && $to_remove{$f2};
	$align1 = $alignments{$f1}->seqs;
	$align2 = $alignments{$f2}->seqs;
		
	$olap_score1=0;
	$olap_score2=0;
	#detect overlap
	foreach $id (keys %$align1) {
	    next if (!exists $align2->{$id});
	    $motif1 = $align1->{$id};
	    $motif2 = $align2->{$id};
	    $len1 = abs($motif1->end - $motif1->start);
	    $len2 = abs($motif2->end - $motif2->start);
	    $overlap1 = 0;
	    $overlap2 = 0;
	    if ($motif1->start > $motif2->start) {
		if ($motif1->end < $motif2->end) { #overlap 
		    $overlap1 =1
		}
		else{ 
		    $overlap = $motif2->end - $motif1->start;
		    if ( $overlap > 0.7 * $len1 && $len1-$overlap < 25){
			$overlap1=1;
		    }
		    if ($overlap > 0.7 * $len2 && $len2-$overlap < 25){
			$overlap2=1;
		    }
		}		    
	    }
	    else{
		if ($motif1->end < $motif2->end  ) { 
		    $overlap = $motif1->end - $motif2->start;
		    if ( $overlap > 0.7 * $len1 && $len1-$overlap < 25){
			$overlap1=1;
		    }
		    if ($overlap > 0.7 * $len2 && $len2-$overlap < 25){
			$overlap2=1;
		    }
		}
		else{ #overlap 
		    $overlap2 =1;
		}
	    }
	    $olap_score1 += $overlap1 * $motif1->score * $motif1->weight;
	    $olap_score2 += $overlap2 * $motif2->score * $motif2->weight;		    	    	
	}
	#print "$f1 ", $all_score{$f1}, " $olap_score1 \t",  "$f2 ", $all_score{$f2}, " $olap_score2 \n";	
	
	if ($olap_score1 > 0.7 * $all_score{$f1} &&  better($all_stat{$f1}, $all_stat{$f2})){
	    $to_remove{$f1} = 1;
	    print "remove $f1\n";
	}		    
	if ($olap_score2 > 0.7 * $all_score{$f2} &&  better($all_stat{$f1}, $all_stat{$f2})){
	    $to_remove{$f2} = 1;	
	    print "remove $f2\n";	
	}		    
    }
    
}
#print join " ", (sort keys %to_remove), "\n";
foreach $f (sort keys %to_remove){
    `echo $f >> $dump_file`;
    `cat $f >> $dump_file`;
    `rm -f $f`;
}
