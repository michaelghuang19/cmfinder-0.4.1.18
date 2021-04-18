#!/usr/bin/perl -w

use Class::Struct;
$path= $ENV{CMfinder};
$bin_path = "$path/bin";
do "$bin_path/io.pl";

$MAX_GAP = 100;

sub count_helix{
    my ($pos1, $pt, $len)= @_;
    return 0 if !exists $pt->{$pos1};
    my $pos2 = $pt->{$pos1};
    if ($pos1 > $pos2) {
	$temp = $pos1;
	$pos1 = $pos2;
	$pos2 = $temp;
    }
    $helix_outer = 0;
    for($i=$pos2 + 1; $i < $len; $i++) {
	last if !exists $pt->{$i};
	last if $pt->{$i} != $pt->{$i-1} - 1;	
	$helix_outer++;
    }
    $helix_inner =0;
    for($i=$pos2 - 1; $i < $len; $i++) {
	last if !exists $pt->{$i};
	last if $pt->{$i} != $pt->{$i+1} + 1;	
	$helix_inner++;
    }
    return ($helix_outer, $helix_inner);
}


sub resolve_overlap{
    my ($alignment1, $alignmen2)=@_;
    my $cost1=0;
    my $cost2=0;
    my $align1 = $alignment1->seqs;
    my $align2 = $alignment2->seqs;

    foreach $id (keys %$align1) {
	if (exists $align2->{$id}) {
	    $motif1 = $align1->{$id};
	    $motif2 = $align2->{$id};
	    $seq1 = $motif1->align_seq;
	    $seq2 = $motif2->align_seq;
	    $ss1 = $motif1->align_ss;
	    $ss2 = $motif2->align_ss;    
    my $map1 = $motif1->align_map;
    my $map2 = $motif2->align_map;

	    if ($motif2->start <= $motif1->end) {
		$olap_start = $motif2->start;
		$olap_end = $motif1->end;
		%pt2 = pair_table($ss2);
		$conflict2 =0;		
		for($i=0; $i < length($seq2); $i++) {
#		    if (defined($map2->{$i})) {			die;		    }		    print "map2->i = \"".(defined($map2->{$i})?$map2->{$i}:"undef")."\"\n";
		    next if !exists $map2->{$i};
		    last if ($map2->{$i} > $olap_end);
		    next if !exists $pt2{$i};
		    if ($pt2{$i} >= 0) {
			$conflict2++; 
		    }
		}
		$conflict1 =0;
		%pt1 = pair_table($ss1);
		for($i= length($seq1) -1; $i>= 0; $i--) {
#		    if (defined($map1->{i})) { die; }  print "map1->i = \"".(defined($map1->{$i})?$map1->{$i}:"undef")."\"\n";
		    next if !exists $map1->{$i};
		    last if ($map1->{$i} < $olap_start);
		    next if !exists $pt1{$i};
		    if ($pt1{$i} >= 0) {
			$conflict1++;
		    }		    
		}
		$cost1 += $conflict1 * $motif1->weight;
		$cost2 += $conflict2 * $motif2->weight;		
	    }
	}
    }
    return 1 if ($cost1 < $cost2);	
    return 2;
}

sub merge_motif{
    #assume that $motif1 should be before $motif2;
    my ($motif1, $motif2, $whole_seq, $olap_own) = @_;
    my  $merge_motif;
    my ($seq1, $gap_seq, $seq2);
    my ($ss1, $gap_ss, $ss2);
    my ($start, $end);

    if ($motif1 eq "") {
       return ($motif2->start, $motif2->end, "", "", "", "", $motif2->align_seq, $motif2->align_ss);
    }
    elsif ($motif2 eq "") {
    	return ($motif1->start, $motif1->end,  $motif1->align_seq, $motif1->align_ss, "", "","", "");
    }
    
    # Two motifs can't be merged.    
    if ($motif1 eq "" || $motif2 eq "" || $motif1->start > $motif2->start){
	return (-1, -1, "", "", "", "", "", "") ;
    }
    
    $seq1 = $motif1->align_seq;
    $seq2 = $motif2->align_seq;
    $ss1 = $motif1->align_ss;
    $ss2 = $motif2->align_ss;    
    $start = $motif1->start;
    $end = $motif2->end;

    #motif2 is contained in motif1
    if ($motif2->end < $motif1->end) {  
	#return ($motif1->start, $motif1->end, $seq1, $ss1, "", "", "","");
	return (-1, -1, "", "", "", "", "", "") ;
    }

    my $map1 = $motif1->align_map;
    my $map2 = $motif2->align_map;

    if ($motif2->start <= $motif1->end) {
	$olap_start = $motif2->start;
	$olap_end = $motif1->end;	
	#assign the overlap region to motif1, remove it from motif2
	if($olap_own == 1){
	    %pt = pair_table($ss2);
	    for($i=0; $i < length($seq2); $i++) {
		next if !exists $map2->{$i};
		last if ($map2->{$i} > $olap_end);
		substr $seq2, $i, 1, '.';
		substr $ss2, $i, 1, '.';
		next if !exists $pt{$i};
		if ($pt{$i} >= 0) {
		    substr $ss2, $pt{$i}, 1, '.';
		}
	    }
	}
	#assign the overlap region to motif2, remove it from motif1
	else{
	    %pt = pair_table($ss1);
	    for($i= length($seq1) -1; $i>= 0; $i--) {
		next if !exists $map1->{$i};
		last if ($map1->{$i} < $olap_start);
		substr $seq1, $i, 1, '.';
		substr $ss1, $i, 1, '.';
		next if !exists $pt{$i};
		if ($pt{$i} >= 0) {
		    substr $ss1, $pt{$i}, 1, '.';
		}
	    } 
	}
	$gap_seq="";
	$gap_ss ="";
    }
    #No overlap
    else{
	if ($motif2->start == $motif1->end + 1){
	    $gap_seq="";
	    $gap_ss="";
	}
	else{

	    $gap_seq = substr($whole_seq, $motif1->end, $motif2->start - $motif1->end - 1);

	    if (length($gap_seq) > $MAX_GAP) {
		if($motif1->score > $motif2->score){
		    $start = $motif1->start;
		    $end   = $motif1->end;
		    $gap_seq ="";
		    $gap_ss = "";
		    $seq2 =  make_string(length($seq2), '.');
		    $ss2 =   make_string(length($seq2), '.');
		}
		else{
		    $start = $motif2->start;
		    $end   = $motif2->end;
		    $gap_seq ="";
		    $gap_ss = "";
		    $seq1 =  make_string(length($seq1), '.');
		    $ss1 =  make_string(length($seq1), '.');
		    
		}
	    }
	    $gap_ss= make_string(length($gap_seq), '.');
	}
    }   
    $gap_seq =~ s/[tT]/U/g;  
    return ($start, $end, $seq1, $ss1, $gap_seq, $gap_ss, $seq2, $ss2);
}


sub merge_alignment{
    my ($alignment1, $alignment2, $seqs) = @_;
    my ($align1, $align2, $ss_cons1, $ss_cons2, $rf1, $rf2) = 
	($alignment1->seqs, $alignment2->seqs, 
	 $alignment1->ss_cons, $alignment2->ss_cons,
	 $alignment1->rf, $alignment2->rf);
	 
    my %motif_overlap = ();
    my $max_gap_len = 0;
    my $max_m1_len = 0;
    my $max_m2_len = 0;
    my $avg_gap_len=0;
    my ($start, $end, $seq1, $ss1, $gap_seq, $gap_ss, $seq2, $ss2);
    my %merged_motif = ();
    my $merged_ss_cons;
    my $merged_rf;
    my $olap_own = resolve_overlap($alignment1, $alignment2);

    my %ids = ();
    $max_id = 0;
    foreach $id (keys %$align1) {	
	align_seq_map($align1->{$id});
	$ids{$id} = $align1->{$id}->id;
	if ($ids{$id} > $max_id) {
	    $max_id = $ids{$id};
	}
    }

    foreach $id (keys %$align2) {
	align_seq_map($align2->{$id});
	next if ($ids{$id});
	$max_id++;
	$ids{$id} = $max_id;
    }
    foreach $id ( keys %ids) {
	if (exists $align2->{$id}) {
	    $acc = $align2->{$id}->acc;

	}
	else{
	    $acc = $align1->{$id}->acc;
	}		
	if (!exists $seqs->{$id}){
	    die "$id does not exist in input sequences ";
	}
	if (!exists $align2->{$id}) {
	    $weight = $align1->{$id}->weight;
	    ($start, $end, $seq1, $ss1, $gap_seq, $gap_ss, $seq2, $ss2) = 
		merge_motif($align1->{$id}, "", $seqs->{$id}->seq, $olap_own);	   
	}
	elsif (!exists $align1->{$id}) {
	    $weight = $align2->{$id}->weight;
	    ($start, $end, $seq1, $ss1, $gap_seq, $gap_ss, $seq2, $ss2) = 
		merge_motif("", $align2->{$id}, $seqs->{$id}->seq, $olap_own);	   
	}		
	else{
	    $weight = $align2->{$id}->weight + $align1->{$id}->weight;
	    ($start, $end, $seq1, $ss1, $gap_seq, $gap_ss, $seq2, $ss2) = 
		merge_motif($align1->{$id}, $align2->{$id}, $seqs->{$id}->seq, $olap_own);	   
	}
	if ($start >= 0 && $end >=0) {	    
	    $motif_overlap{$id} = {weight => $weight/2, start=> $start,     end => $end, 
				   seq1 => $seq1,      ss1=> $ss1, 
				   gap_seq=> $gap_seq, gap_ss => $gap_ss, 
				   seq2=> $seq2,       ss2=> $ss2};	    
	    if ($max_gap_len < length($gap_seq)){
		$max_gap_len = length($gap_seq);
	    }
	    if ($max_m1_len < length($seq1)){
		$max_m1_len = length($seq1);
	    }
	    if ($max_m2_len < length($seq2)){
		$max_m2_len = length($seq2);
	    }
	}
    }

    return if (scalar keys %motif_overlap < 2) ;
    $avg_gap_len = $avg_gap_len / (scalar keys %motif_overlap);

    open OUT, ">$out.gap";	
    $gap_count=0;
    foreach $id (keys %motif_overlap) {
	$gap = $motif_overlap{$id}{gap_seq};
	if (length($gap) > 0){
	    print OUT ">$id\n";
	    print OUT "$gap\n";
	    $gap_count ++;
	}
    }
    close OUT;
    
    #system("fold_seq $out.gap $out.gap.cand");
    #system("cands -n 1 -f 1 $out.gap $out.gap.cand");
    #system("canda $out.gap.cand_1 $out.gap $out.gap.align");
    if ($gap_count > 1){
	`$bin_path/clustalw -infile=$out.gap -outfile=$out.gap.aln`;    
	system("$bin_path/sreformat stockholm $out.gap.aln > $out.gap.align");    
	$gap_sto = read_stockholm("$out.gap.align");    
	$gap_align  = $gap_sto ->seqs;
	foreach $id (keys %motif_overlap) {
	    if (exists $gap_align->{$id}){
		$gap_seq = $gap_align->{$id}->align_seq;
		$motif_overlap{$id}{gap_seq}=$gap_seq;
		$motif_overlap{$id}{gap_ss} = "";
		if ($max_gap_len < length($gap_seq)){
		    $max_gap_len = length($gap_seq);
		}
	    }	    
	}
	system("rm $out.gap*");
    }

    $merged_ss_cons = $ss_cons1.pad_string("", $max_gap_len, '.',1).$ss_cons2;
    $merged_rf = $rf1.pad_string("", $max_gap_len, '.',1).$rf2;

    foreach $id (keys %motif_overlap) {

	$gap_seq = pad_string($motif_overlap{$id}{gap_seq}, $max_gap_len, '.',1);
	$gap_ss =  pad_string($motif_overlap{$id}{gap_ss},  $max_gap_len, '.',1);
	$align_seq = pad_string($motif_overlap{$id}{seq1},$max_m1_len, '.',1).$gap_seq.pad_string($motif_overlap{$id}{seq2}, $max_m2_len, '.',1);
	$align_ss  = pad_string($motif_overlap{$id}{ss1}, $max_m1_len, '.',1).$gap_ss.pad_string($motif_overlap{$id}{ss2},  $max_m2_len, '.',1);
	$score1 = $score2 = 0;
	$score1 = $align1->{$id}->score if (exists $align1->{$id});
	$score2 = $align2->{$id}->score if (exists $align2->{$id});
	$start = $motif_overlap{$id}{start};
	$end = $motif_overlap{$id}{end};
	$score = $score1 + $score2 - length(remove_gap($gap_seq));
	$weight = $motif_overlap{$id}{weight};
	$desc = sprintf "%3d..%3d\t %.3f", $start, $end, $score;
	$merged_motif{$id} = AlignSeq->new(acc=>$id, id=> $ids{$id}, 
					   start => $start, end=>$end, desc=>$desc,
					   score=>$score, weight => $weight, 
					   align_seq=>$align_seq, align_ss=>$align_ss);		

    }    
    return Alignment->new(seqs=>\%merged_motif, 
			  ss_cons=>$merged_ss_cons, 
			  rf=>$merged_rf, 
			  flags=>$alignment1->flags);    
}

($seq_file, $ali_f1, $ali_f2, $out) = @ARGV;


$seqs = read_fasta($seq_file);
$alignment1= read_stockholm($ali_f1);
$alignment2= read_stockholm($ali_f2);


$new_alignment = merge_alignment($alignment1, $alignment2, $seqs);
write_stockholm($new_alignment, $out);





