#!/usr/bin/perl -w

use Class::Struct;
use Getopt::Long;
$path= $ENV{CMfinder};
$bin_path = "$path/bin";
do "$bin_path/io.pl";

print STDERR "CMfinder $path\n";

$MAX_GAP = 100;
$len_energy_threshold = 0.1;
$min_overlap=2;
$min_num=2.5;
$max_len = 200;

struct 'MergeMotif' => {motif1=> '$', motif2=>'$', num_seq=>'$', score=>'$', gap=>'$', 
			overlap=>'$', weight=>'$'};



$cluster_file = "";
$tree_file = "";
$hmm = 0;
$banded = 0;
$weight_option="";
$default=0;
$verbose=0;
GetOptions(
	   "r=s" => \$cluster_file,
	   "t=s" => \$tree_file,
	   "o=i" => \$min_overlap,
	   "n=i" => \$min_num,
           "w=s" => \$weight_option,
	   "hmm" => \$hmm,
           "hbanded"=>\$banded,
           "def" =>\$default, 	 
	   "v" => \$verbose,
	   "isAReRun" => \$isAReRun
	   );
$tree_option = "";
if ($tree_file ne ""){
    $tree_option = "--t $tree_file";
}
$hmm_option="";
if ($hmm){
    $hmm_option = "--hmm";
}
$cluster_option="";
if ($cluster_file ne ""){
    $cluster_option="-r $cluster_file";
}

$banded_option="";
if ($banded){
   $banded_option = "--hbanded";
}

if ($weight_option ne ""){
   $weight_option = "-w $weight_option";
}

# $verbose_flag=$verbose ? "-v" : "";

sub my_strcmp
{
    ($s1, $s2)= @_;
    @t1 = split /\./, $s1;
    @t2 = split /\./, $s2;
    $prefix="";
    for( $i = 0; $i < scalar @t1 && $i < scalar @t2; $i++) {
	last if ($t1[$i] ne $t2[$i]);	
	if ( $prefix eq "") {
	    $prefix =$t1[$i];	
	}
	else{
	    $prefix ="$prefix.".$t1[$i];	
	}
    }
    $suffix1 = join ".", @t1[$i..($#t1)];
    $suffix2 = join ".", @t2[$i..($#t2)];    
    return ($prefix, $suffix1, $suffix2);
}	  

sub match_file{  
    my @files = ();
    $f = shift @_;
    while(<$f*>) {
	if ($isAReRun) {
	    if ((/[.]html$/) || (/[.]sto$/) || (/[.]temp$/) || (/[.]gap$/)) {
		next;
	    }
	    if (!(/motif[.]h[12]_[0-9]+$/)) { # only simple motifs
		next;
	    }
	}
	push @files, $_;
	print "match_file: $_\n";
    }
    return @files;
}


$seq_file = shift @ARGV;
@files = @ARGV;
@align_files = ();
@alignments = ();

foreach $f (@files) {
    push @align_files, match_file($f);
}

if ($verbose) { print "align_files=\n".join("",map {"\t$_\n"} @align_files); }

if ($verbose) { print "about to read alignments\n"; }

@all_files = @align_files;
%all_stats=();
foreach $f (@align_files) {
    if ($verbose) { print "\treading $f\n"; }
    $align = read_stockholm($f);  
    if ($align->weight >= $min_num){
	if ($verbose) { print "Running $bin_path/summarize $f\n"; }
	$summary = `$bin_path/summarize $f`;
	my %stat = $summary=~ /\s*(\S+)=(\S+)\s*/g;
	$all_stats{$f} = \%stat;
	$alignments{$f}= $align;
    }
}

if ($verbose) { print "\tdone reading alignments\n"; }

%merge_motif=();


sub try_merge{
    my $f1 = shift @_;
    my $f2 = shift @_;    
    my $index = join ".", my_strcmp($f1, $f2);	    
    print "try_merge $f1 $f2\n";
    if (exists $merge_motif{$index}) {
	if ($verbose) { print "\texists \$merge_motif{$index}\n"; }
	return;
    }
    my $num_overlap = 0;
    my $start1 = 0;
    my $start2 = 0;
    my $start1_score=0;
    my $start2_score=0;
    my $gap1=0;
    my $gap2=0;
    my $overlap1 = 0;
    my $overlap2 = 0;
	    
    #detect overlap
    if (!exists $alignments{$f1} || !exists $alignments{$f2}) {
	if ($verbose) { print "!exists \$alignments{$f1} || !exists \$alignments{$f2}\n";}
	return;
    }
    $align1 = $alignments{$f1}->seqs;
    $align2 = $alignments{$f2}->seqs;
    foreach $id (keys %$align1) {		
	next if (!exists $align2->{$id});
	$motif1 = $align1->{$id};
	$motif2 = $align2->{$id};
	$len1 = abs($motif1->end - $motif1->start);
	$len2 = abs($motif2->end - $motif2->start);
	#print "$id ", $motif1->start, "-", $motif1->end, "\t", $motif2->start, "-", $motif2->end, "\n";
	if ($motif1->start > $motif2->start) {
	    if ($motif1->end < $motif2->end) { #overlap 
		$num_overlap += $motif1->weight * $motif2->weight;
	    }
	    else{ 
		$overlap = $motif2->end - $motif1->start;
		#print "$id $overlap\n";
		if ($overlap > - $MAX_GAP){
		    if(($len1 - $overlap < 25 || $len2 - $overlap < 25) &&
		       ($overlap > 0.9 * $len1 || $overlap > 0.9 * $len2)){
			$num_overlap+=  $motif1->weight * $motif2->weight;
		    }
		    else{
			$start2 +=  ($motif1->weight + $motif2->weight)/2;
			$start2_score += $motif1->score * $motif1->weight + $motif2->score* $motif2->weight;
			if ($overlap > 0){
			    $overlap2 += $overlap *  $motif1->weight * $motif2->weight;
			}
			elsif ( - $overlap < $MAX_GAP){
			    $gap2 += - $overlap * $motif1->weight * $motif2->weight;
			}
		    }
		}		
	    }    
	}
	else{
	    if ($motif1->end < $motif2->end) { 
		$overlap = $motif1->end - $motif2->start;
		#print "$id $overlap\n";
		if ($overlap > - $MAX_GAP){
		    if(($len1 - $overlap < 25 || $len2 - $overlap < 25) &&
		       ($overlap > 0.9 * $len1 || $overlap > 0.9 * $len2)){
			$num_overlap+=  $motif1->weight * $motif2->weight;
		    }
		    else{
			$start1+=  ($motif1->weight + $motif2->weight)/2;
			$start1_score += $motif1->score * $motif1->weight + $motif2->score * $motif2->weight;
			if ($overlap > 0){
			    $overlap1 += $overlap *  $motif1->weight * $motif2->weight;
			}
			elsif ( - $overlap1 < $MAX_GAP){
			    $gap1 += - $overlap *  $motif1->weight * $motif2->weight;
			}
		    }
		}
	    }
	    else{ #overlap 
		$num_overlap +=  $motif1->weight * $motif2->weight;
	    }
	}		
    }
    if ($verbose) { print "\t( start1,start2,num_overlap = $start1,$start2,$num_overlap )\n"; }
    if ($num_overlap > $start1 && $num_overlap > $start2) {
	if ($verbose) { print "\$num_overlap > \$start1 && \$num_overlap > \$start2: $num_overlap > $start1 && $num_overlap > $start2\n"; }
	return;
    }
    if ($start1_score > $start2_score) {
	if ($start1 < $min_overlap) {
	    if ($verbose) { print "\t\$start1 < \$min_overlap : $start1 < $min_overlap\n"; }
	    return;
	}
	#motif1 is before motif2
	$merge_motif{$index} = 
	    MergeMotif->new(motif1 => $f1, 
			    motif2 => $f2, 
			    num_seq=> $start1, 
			    score => $start1_score, 
			    gap => $gap1/$start1, 
			    overlap=>$overlap1/$start1, 
			    weight=> $start1_score - $gap1/2 - $overlap1, 
			    );
	print "\taccept id=$index: $f1 - $f2: $start1_score  $gap1/2  $overlap1  $start1\n";
    }
    else{
	if ($start2 < $min_overlap) {
	    if ($verbose) { print "\t\$start2 < \$min_overlap : $start2 < $min_overlap\n"; }
	    return;
	}
	#motif2 is before motif1
	$merge_motif{$index} = 
	    MergeMotif->new(motif1 => $f2, 
			    motif2 => $f1, 
			    num_seq => $start2,
			    score => $start2_score, 
			    gap => $gap2/$start2, 
			    overlap=>$overlap2/$start2,
			    weight=> $start2_score - $gap2/2 - $overlap2, 
			    );
	print "\taccept id=$index: $f2 - $f1: $start2_score   $gap2/2  $overlap2  $start2\n";
    }		     
}

foreach $f1 (@align_files){ 
    foreach $f2 (@align_files){ 
	if ($f1 lt $f2){
	    try_merge($f1,$f2);
	}
    }
}

%processed=();
while(scalar keys %merge_motif > 0){ 
    if ($verbose) { print "entering for my \$id, with list:\n".join("",map {"\t$_\n"} reverse sort {$merge_motif{$a}->weight <=> $merge_motif{$b}->weight} keys %merge_motif); }
    foreach $id (reverse sort {$merge_motif{$a}->weight <=> $merge_motif{$b}->weight } keys %merge_motif) {
	if ($verbose) { print "while keys merge_motif : id=$id\n"; }
	$m = $merge_motif{$id};
	delete($merge_motif{$id});
	if (exists $processed{$id}) {
	    if ($verbose) { print "\twhile keys merge_motif NEXT : exists \$processed{$id}\n"; }
	    next;
	}
	$processed{$id}=$m;
	$f1 = $m->motif1;
	$f2 = $m->motif2;	
			
	if ($m->weight <= 0) {
	    if ($verbose) { print "\twhile keys merge_motif NEXT : weight <= 0\n"; }
	}
	if ($m->weight > 0) {
	    if (exists $merged_files{$f1} || exists $merged_files{$f2}) {
		if ($verbose) { print "\twhile keys merge_motif NEXT : exists \$merged_files{$f1}=".(exists $merged_files{$f1}?"true":"false")." || exists \$merged_files{$f2}=".(exists $merged_files{$f2}?"true":"false")."\n"; }
		next;
	    }
	    if ($m->gap > $MAX_GAP) {
		my $m_gap=$m->gap;
		if ($verbose) { print "\twhile keys merge_motif NEXT : \$m->gap > \$comb_max_gap : $m_gap > $MAX_GAP\n"; }
		next;
	    }
	    $f = $id;
	    $found = 0;
	    foreach $tmp (@all_files){
		if ($tmp eq $f){
		    $found = 1;
		    last;
		}
	    }
	    if ($found) {
		if ($verbose) { print "\twhile keys merge_motif NEXT : not found.  all_files = ".join(" ",@all_files)."\n"; }
		next;
	    }

	    print "( near merge_motif: $id, ",$m->num_seq,"\t$f1\t$f2\t",$m->weight,"\t",$m->gap, "\t",$m->overlap," )\n";
	    print ("perl $bin_path/merge_motif.pl $seq_file $f1 $f2 $f.temp\n");
	    $status=system("perl $bin_path/merge_motif.pl $seq_file $f1 $f2 $f.temp");
	    if ($status!=0) {
		# die "problem with merge_motif.pl";
	    }
	    $cm_file = $f;
	    $cm_file =~ s/\.motif/\.cm/;

	    if ($default){
		$tmp = `$bin_path/summarize $f.temp`;
		print "$f.temp $tmp\n";
		if ($tmp =~ /Len=(\S+)/){
		    $len =$1;
		    if ($len < 50){
			$hmm_option = "";
			$banded_option = "";
			$cluster_option = "";
		    }
		    else{
			$hmm_option = "--hmm";
			$banded_option = "--hbanded";
		    }
		}			
	    }	    
	    $cmd= "$bin_path/cmfinder $cluster_option $tree_option $hmm_option $banded_option $weight_option -o $f -a $f.temp $seq_file $cm_file";
	    print "$cmd\n";
	    $status = system($cmd);	    
	    if ($status != 0){
		if (-e $f){
		    system("rm $f");
		}
		if (-e $cm_file){
		    system("rm $cm_file");
		}
	    }
	    system("rm -f $f.temp*");
	    if (!-s $f) {
		if ($verbose) { print "\twhile keys merge_motif NEXT : !-s $f\n"; }
		next;
	    }
	    $summary = `$bin_path/summarize $f`;	    
	    my %stat = $summary=~ /\s*(\S+)=(\S+)\s*/g;
	    $all_stats{$f} = \%stat;	    
	    
	    if ($stat{"Weight"} < $min_num || 
		$stat{"Score"} < $all_stats{$f1}->{"Score"} && $stat{"Weight"} <= $all_stats{$f1}->{"Weight"}  ||
		$stat{"Score"} < $all_stats{$f2}->{"Score"} && $stat{"Weight"} <= $all_stats{$f2}->{"Weight"}  ||
		$stat{"BP.org"} < $all_stats{$f1}->{"BP.org"} + 3 && $stat{"Weight"} <= $all_stats{$f1}->{"Weight"}  ||
		$stat{"BP.org"} < $all_stats{$f2}->{"BP.org"} + 3 && $stat{"Weight"} <= $all_stats{$f2}->{"Weight"}){
		if ($verbose) { print "\twhile keys merge_motif NEXT : score/weight/BP.org not good enough\n"; }
		print "\t\tRemove $f\n";
		print "\t\tremove-info $f: ", $stat{"Weight"}, "\t", $stat{"Score"}, "\t", $stat{"BP"}, "\n";
		print "\t\tremove-info $f1:", $all_stats{$f1}->{"Weight"}, "\t", $all_stats{$f1}->{"Score"},"\t", $all_stats{$f1}->{"BP"}, "\n";		
		print "\t\tremove-info $f2:", $all_stats{$f2}->{"Weight"}, "\t", $all_stats{$f2}->{"Score"},"\t", $all_stats{$f2}->{"BP"}, "\n";		
	    	system("rm -f $f");
		system("rm -f $cm_file");
	    }
	    else{
		if ($verbose) { print "\t\tset merged_files{$f1}=1\n";}
		if ($verbose) { print "\t\tset merged_files{$f2}=1\n";}
		$merged_files{$f1} = 1;
		$merged_files{$f2} = 1;	
		if ($stat{"Energy"} > $stat{"Len"} * $len_energy_threshold ||$stat{"Len"} > $max_len) {
		    if ($verbose) { print "\twhile keys merge_motif NEXT : \$stat{\"Energy\"} > \$stat{\"Len\"} * \$comb_len_energy_threshold || \$stat{\"Len\"} > \$comb_max_len\n"; }
		    next;
		}
		$merge_align = read_stockholm($f);
		$alignments{$f}= $merge_align;
		foreach $f3 (@align_files){
		    if (!exists $merged_files{$f3}){
			if ($verbose) { print "\twhile keys merge_motif : try_merge in context of id=$id\n"; }
			try_merge($f, $f3);
		    }
		    else {
			if ($verbose) { print "\twhile keys merge_motif : !exists \$merged_files{$f3}, so skipping\n"; }
		    }
		}
		push @align_files, $f;
		last;
	    }		
	}
    }
}
    


