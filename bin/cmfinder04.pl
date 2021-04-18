#!/usr/bin/perl -w

use strict;
use Class::Struct;
use Getopt::Long qw(:config no_ignore_case);
use FindBin qw($Bin);
my $bin_path=$Bin;
do "$bin_path/io.pl";

my $CMFINDER_PACKAGE_VERSION="0.4.1.18";

#default parameters

my $CAND=40;
my $MAXSPAN1=100;
my $MINSPAN1=30;
my $MAXSPAN2=100;
my $MINSPAN2=40;
my $CLUSTER = 3;
my $FRACTION=0.8;
my $SINGLE = 5;
my $DOUBLE = 5;
my $verbose = 0;
my $help = 0;
my $COMBINE = 0;
my $cand_weight_option="";
my $cmfinderBaseExe="cmfinder04";
my ($likeold,$skipClustalw,$useOldCmfinder,$simpleMotifsAlreadyDone,$justGetCmfinderCommand,$copyCmfinderRunsFromLog,$amaa,$version,$filterNonFrag,$fragmentary,$commaSepEmFlags,$commaSepSummarizeFlags,$commaSepCandfFlags,$saveTimer,$allCpus,$cpu,$candsParallel,$outFileSuffix,$columnOnlyBasePairProbs,$fragmentaryCommaSepEmFlags);
my $emulate_apparent_bug_in_resolve_overlap=1;

# parameters from comb_motif.pl
my $comb_max_gap=100; # motif instances within this distance of each other can be considered close enough for merging
my $comb_len_energy_threshold = 0.1;
my $comb_min_overlap=2;
my $comb_min_num=2.5;
my $comb_max_len = 200;
my $output_more_like_old_cmfinder_pl=0;
my $motifList;
my $minCandScoreInFinal=0; # be more like old cmfinder for now
my $emSeq;

struct 'MergeMotif' => {motif1=> '$', motif2=>'$', num_seq=>'$', score=>'$', gap=>'$', overlap=>'$', weight=>'$'};

if (!GetOptions(
	 "version" => \$version,
	 "h" => \$help,
	 "v" => \$verbose,
	 "w=s" => \$cand_weight_option,
	 "c=i" => \$CAND,	
	 "minspan1=i" => \$MINSPAN1,
	 "maxspan1=i" => \$MAXSPAN1,
	 "minspan2=i" => \$MINSPAN2,
	 "maxspan2=i" => \$MAXSPAN2,
	 "f=f" => \$FRACTION,
	 "s1=i" => \$SINGLE,
	 "s2=i" => \$DOUBLE,
	 "combine" =>\$COMBINE,
	 "o=i" => \$comb_min_overlap,
	 "n=i" => \$comb_min_num,
	 "skipClustalw" => \$skipClustalw,
	 "skipClustalW" => \$skipClustalw,
	 "likeold" => \$likeold,
	 "useOldCmfinder" => \$useOldCmfinder,
	 "simpleMotifsAlreadyDone" => \$simpleMotifsAlreadyDone,
	 "justGetCmfinderCommand" => \$justGetCmfinderCommand,
	 "copyCmfinderRunsFromLog=s" => \$copyCmfinderRunsFromLog,
	 "amaa" => \$amaa,
	 "motifList=s" => \$motifList,
	 "minCandScoreInFinal=s" => \$minCandScoreInFinal,
	 "emSeq=s" => \$emSeq,
	 "filterNonFrag" => \$filterNonFrag,
		 "fragmentary" => \$fragmentary,
		 "fragmentary-commaSepEmFlags" => \$fragmentaryCommaSepEmFlags,
	 "commaSepEmFlags=s" => \$commaSepEmFlags,
	 "commaSepSummarizeFlags=s" => \$commaSepSummarizeFlags,
	 "commaSepCandfFlags=s" => \$commaSepCandfFlags,
	 "saveTimer=s" => \$saveTimer,
	 "allCpus" => \$allCpus,
	 "cpu=s" => \$cpu,
	 "candsParallel" => \$candsParallel,
	 "outFileSuffix=s" => \$outFileSuffix,
	 "columnOnlyBasePairProbs" => \$columnOnlyBasePairProbs
    )){
    print STDERR "Invalid options\n";
    print_help();
    exit(1);
}

if ($fragmentaryCommaSepEmFlags) {
	die "please separate the flags -fragmentary and -commaSepEmFlags by a space; they're separate flags";
}

if ($help) {
    print_help();
    exit(0);
}

if ($version) {
    print "CMFINDER_PACKAGE_VERSION=$CMFINDER_PACKAGE_VERSION.\n";
    exit(0);
}

my $SEQ=shift @ARGV;
my $seqForExpectationMaximization=$SEQ;
if ($emSeq) {
    $seqForExpectationMaximization=$emSeq;
}

if (!$outFileSuffix) {
    $outFileSuffix="";
}

my $dummyCmfileParamForCmfinder="";
if ($useOldCmfinder) {
    $cmfinderBaseExe="cmfinder";
    $dummyCmfileParamForCmfinder="$SEQ.temp.cm";
}

my @cmfinder_inf11FlagsList=();
my @summarizeFlagsList=();
if ($allCpus) {
    $cpu=-1; # cmfinder interprets this as use all CPUs
}
if ($cpu) {
    push @cmfinder_inf11FlagsList,"--cpu $cpu";
}
if ($likeold && !$useOldCmfinder) {
    push @cmfinder_inf11FlagsList,"--enone";
    push @cmfinder_inf11FlagsList,"--p56";
    push @cmfinder_inf11FlagsList,"--degen-rand";
    push @cmfinder_inf11FlagsList,"--ints-like-03";
    push @cmfinder_inf11FlagsList,"--min-seq-weight 0";
    push @cmfinder_inf11FlagsList,"--no-elim-iden-seq";
    push @cmfinder_inf11FlagsList,"--no-elim-iden-subseq";
    push @cmfinder_inf11FlagsList,"--min-cand-score-in-final 0";
}
push @cmfinder_inf11FlagsList,"--min-cand-score-in-final $minCandScoreInFinal";
if ($amaa) {
    push @cmfinder_inf11FlagsList,"--amaa";
}
if ($filterNonFrag) {
    push @cmfinder_inf11FlagsList,"--filter-non-frag";
}
if ($fragmentary) {
    push @cmfinder_inf11FlagsList,"--fragmentary";
    push @summarizeFlagsList,"--fragmentary";
}
my $cmfinder_inf11Flags=join(" ",@cmfinder_inf11FlagsList);
if ($commaSepEmFlags) {
    $commaSepEmFlags =~ s/^x//g;
    my $spacey=$commaSepEmFlags;
    $spacey =~ s/,/ /g;
    $cmfinder_inf11Flags .= " $spacey";
}
if ($columnOnlyBasePairProbs) {
    $cmfinder_inf11Flags .= " --column-only-base-pair-probs";
}
my $summarizeFlagsStr=join(" ",@summarizeFlagsList);
if ($commaSepSummarizeFlags) {
    $commaSepSummarizeFlags =~ s/^x//g;
    my $spacey=$commaSepSummarizeFlags;
    $spacey =~ s/,/ /g;
    $summarizeFlagsStr .= " $spacey";
}
my $candfExtraFlags="";
if ($commaSepCandfFlags) {
    $commaSepCandfFlags =~ s/^x//g;
    my $spacey=$commaSepCandfFlags;
    $spacey =~ s/,/ /g;
    $candfExtraFlags .= " $spacey";
}
if ($cand_weight_option ne ""){
    $cand_weight_option = "-w $cand_weight_option";
}
my $saveTimerFlag="";
my $saveTimer03Flag="";
if ($saveTimer) {
    unlink($saveTimer); # make sure nothing's in the file at the beginning
    $saveTimerFlag="--timer-append $saveTimer";
    $saveTimer03Flag="-t $saveTimer";
}

if ($justGetCmfinderCommand) {
    my $cmfinder_cmd="$bin_path/$cmfinderBaseExe $cmfinder_inf11Flags $cand_weight_option";
    print "-justGetCmfinderCommand:$cmfinder_cmd.\n";
    exit(0);
}

if ($motifList) {
    if (!open(MOTIFLIST,">$motifList")) {
	die "cannot open $motifList";
    }
}

if ($copyCmfinderRunsFromLog) {
    my $dir=$SEQ;
    $dir =~ s/\/[^\/]+$//g;
    print "dir=$dir\n";
    my $gotCmfinderCmd=0;
    if (!open(LOG,$copyCmfinderRunsFromLog)) {
	die "cannot open $copyCmfinderRunsFromLog";
    }
    while (<LOG>) {
	s/[\r\n]//g;
	my $line=$_;
	if ($line =~ /^\S*cmfinder(|_inf11)\s/) {
	    my $version=$1;
	    my ($inputMsa,$inputFasta,$outputFileBase);
	    if ($line =~ /-a\s(\S+)\s/) {
		$inputMsa=$1;
	    }
	    else {
		die;
	    }
	    if ($line =~ /-o\s(\S+)\s/) {
		$outputFileBase=$1;
		$outputFileBase =~ s/^.*\///g;
	    }
	    else {
		die;
	    }
	    my @spaceSepList=split /\s+/,$line;
	    my $cm=pop @spaceSepList;
	    if ($cm =~ /[.]cm[.]/) {
		if ($version eq "_inf11") {
		    die "cmfinder_inf11 should not have cmfile";
		}
		$inputFasta=pop @spaceSepList;
	    }
	    else {
		if (length($version)==0) {
		    die "cmfinder (0.3) should have cmfile";
		}
		$inputFasta=$cm; 
	    }
	    print "ORIGINAL CMD: $line\n";
	    if (!-e $inputMsa) {
		print "\tSKIPPING: input msa \"$inputMsa\" doesn't exist, so I assume it isn't good, and anyway it'd be a hassle to run.\n";
		next;
	    }
	    my $output_motif_file="$bin_path/$cmfinderBaseExe";
	    RunCmfinder("$output_motif_file $cmfinder_inf11Flags $cand_weight_option -o $dir/$outputFileBase -a $inputMsa $inputFasta",$output_motif_file);
	    $gotCmfinderCmd=1;
	}
    }
    close(LOG);
    if (!$gotCmfinderCmd) {
	die "didn't find any cmfinder commands";
    }
    exit(0);
}

if (!defined($SEQ)){
    print STDERR "No sequence file is specified\n";
    print_help();
    exit(1);
}

my $unaligned_seqs=read_fasta($seqForExpectationMaximization);

my $tempFileListFileName="$SEQ.file-list";

my %alignments = ();

my @cands=();
if (!$simpleMotifsAlreadyDone) {
    my $runCandsParallel=0;
    if ((defined($cpu) && $cpu>=2) || $candsParallel || $allCpus) { # for -allCpus, we'll assume that there are >=2 CPUs.  The modules to determine how many actual CPUs there are in Perl don't seem to be standard, and the direct methods seem hacky, so I'd rather just assume >=2 CPUs.
	$runCandsParallel=1;
    }
    my @candsJobs=();
    if ($SINGLE) {
	push @candsJobs,{candf=>"$bin_path/candf $candfExtraFlags $saveTimer03Flag -c $CAND -o $SEQ.cand.h1$outFileSuffix -M $MAXSPAN1 -m $MINSPAN1 -s 1 -S 1 $SEQ",
			 cands=>"$bin_path/cands $saveTimer03Flag -l $tempFileListFileName.single$outFileSuffix -n $SINGLE -f $FRACTION $SEQ $SEQ.cand.h1$outFileSuffix"};
    }
    if ($DOUBLE){
	push @candsJobs,{candf=>"$bin_path/candf $candfExtraFlags $saveTimer03Flag -c $CAND -o $SEQ.cand.h2$outFileSuffix -M $MAXSPAN2 -m $MINSPAN2 -s 2 -S 2 $SEQ",
			 cands=>"$bin_path/cands $saveTimer03Flag -l $tempFileListFileName.double$outFileSuffix -n $DOUBLE -f $FRACTION $SEQ $SEQ.cand.h2$outFileSuffix"};
    }

    # run candf jobs in serial
    my $job;
    if ($runCandsParallel) {
	for my $stepName (qw(candf cands)) {
	    # run candf/s in parallel
	    print "parallel running $stepName\n";
	    for $job (@candsJobs) {
		$job->{childPid}=undef;
		my $childPid=fork();
		if (!defined($childPid)) {
		    die "fork failed";
		}
		if ($childPid==0) {
		    RunCmd($job->{$stepName});
		    exit(0);  # exit child process
		}
		else {
		    $job->{childPid}=$childPid;
		    print "ran childPid=$childPid\n";
		}
	    }
	    # collect children
	    for $job (@candsJobs) {
		print "about to collect childPid=$job->{childPid}\n";
		my $waitpidResult=waitpid($job->{childPid},0);
		my $status=$?;
		if ($waitpidResult==-1) {
		    die "waitpid returned -1 -- where did the child process go?";
		}
		if ($waitpidResult==0) {
		    die "waitpid returned 0 -- but I thought we weren't blocking";
		}
		if ($waitpidResult!=$job->{childPid}) {
		    die "waitpid returned something other than the expected childPid";
		}
		if ($status!=0) {
		    die "child exited with status=$status";
		}
		# okay, at this point, the child exited successfully
	    }
	}
    }
    else {
	# run cands jobs in serial, as requested
	print "running serial\n";
	for $job (@candsJobs) {
	    RunCmd($job->{candf});
	    RunCmd($job->{cands});
	}
    }

    push @cands,GetCandFiles("$tempFileListFileName.single$outFileSuffix");
    push @cands,GetCandFiles("$tempFileListFileName.double$outFileSuffix");
}

my @motifFiles=();
for my $cand (@cands){   
    if ($cand =~ /$SEQ(\.\d)?\.cand.*\_\d+$/){
	my $pref = $1;
	my $align =$cand;
	$align =~ s/cand/align/;
	my $motif = $cand;
	$motif =~ s/cand/motif/;
	if ($simpleMotifsAlreadyDone) {
	    if (-e $motif) {
		push @motifFiles,$motif;
	    }
	}
	else {
	    RunCmd("$bin_path/canda $saveTimer03Flag $cand $SEQ $align");
	    if (RunCmfinder("$bin_path/$cmfinderBaseExe $saveTimerFlag $cmfinder_inf11Flags $cand_weight_option -o $motif -a $align $seqForExpectationMaximization $dummyCmfileParamForCmfinder",$motif)) {
		# produced acceptable output
		push @motifFiles,$motif;
		if ($motifList) {
		    print MOTIFLIST "$motif\n";
		}
	    }
	    else {
		# oh well, nothing acceptable
	    }
	}
    }
    else {
	die "weird cand file name \"$cand\"";
    }
}

if ($COMBINE) {
    CombMotif($cand_weight_option,$seqForExpectationMaximization,\@motifFiles);
}

if (!$verbose){
    system("rm -f $SEQ.*cand*");
    system("rm -f $SEQ.*align*");
    if ($useOldCmfinder) {
	system("rm -f $SEQ.*-*.cm");
    }
}

if ($motifList) {
    close(MOTIFLIST);
}

sub count_helix{
    my ($pos1, $pt, $len)= @_;
    return 0 if !exists $pt->{$pos1};
    my $pos2 = $pt->{$pos1};
    if ($pos1 > $pos2) {
	($pos1,$pos2)=($pos2,$pos1);
    }
    my $helix_outer = 0;
    my $i;
    for($i=$pos2 + 1; $i < $len; $i++) {
	last if !exists $pt->{$i};
	last if $pt->{$i} != $pt->{$i-1} - 1;	
	$helix_outer++;
    }
    my $helix_inner =0;
    for($i=$pos2 - 1; $i < $len; $i++) {
	last if !exists $pt->{$i};
	last if $pt->{$i} != $pt->{$i+1} + 1;	
	$helix_inner++;
    }
    return ($helix_outer, $helix_inner);
}


sub resolve_overlap {
    my ($alignment1, $alignment2)=@_;
    my $i;
    my $cost1=0;
    my $cost2=0;
    my $align1 = $alignment1->seqs;
    my $align2 = $alignment2->seqs;

    for my $id (keys %$align1) {
	if (exists $align2->{$id}) {
	    my $motif1 = $align1->{$id};
	    my $motif2 = $align2->{$id};
	    my $seq1 = $motif1->align_seq;
	    my $seq2 = $motif2->align_seq;
	    my $ss1 = $motif1->align_ss;
	    my $ss2 = $motif2->align_ss;
	    my ($map1,$map2);
	    if ($emulate_apparent_bug_in_resolve_overlap) {
		# skip these inits
	    }
	    else {
		$map1 = $motif1->align_map;
		$map2 = $motif2->align_map;
	    }

	    if ($motif2->start <= $motif1->end) {
		my $olap_start = $motif2->start;
		my $olap_end = $motif1->end;
		my %pt2 = pair_table($ss2);
		my $conflict2 =0;
		for($i=0; $i < length($seq2); $i++) {
		    next if !exists $map2->{$i};
		    last if ($map2->{$i} > $olap_end);
		    next if !exists $pt2{$i};
		    if ($pt2{$i} >= 0) {
			$conflict2++;
		    }
		}
		my $conflict1 =0;
		my %pt1 = pair_table($ss1);
		for($i= length($seq1) -1; $i>= 0; $i--) {
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
    
    my $i;
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
	my $olap_start = $motif2->start;
	my $olap_end = $motif1->end;	
	#assign the overlap region to motif1, remove it from motif2
	if($olap_own == 1){
	    my %pt = pair_table($ss2);
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
	    my %pt = pair_table($ss1);
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

	    if (length($gap_seq) > $comb_max_gap) {
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
    my ($alignment1, $alignment2, $seqs, $out) = @_;
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
    my $max_id = 0;
    for my $id (keys %$align1) {	
	align_seq_map($align1->{$id});
	$ids{$id} = $align1->{$id}->id;
	if ($ids{$id} > $max_id) {
	    $max_id = $ids{$id};
	}
    }

    for my $id (keys %$align2) {
	align_seq_map($align2->{$id});
	next if ($ids{$id});
	$max_id++;
	$ids{$id} = $max_id;
    }
    for my $id ( keys %ids) {
	my ($acc,$weight);
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

    if (!open(OUT,">$out.gap")) {
	die "cannot open $out.gap";
    }
    my $gap_count=0;
    for my $id (keys %motif_overlap) {
	my $gap = $motif_overlap{$id}{gap_seq};
	if (length($gap) > 0){
	    print OUT ">$id\n";
	    print OUT "$gap\n";
	    $gap_count ++;
	}
    }
    close OUT;
    
    if ($gap_count > 1) {
	if ($skipClustalw) {
	    if ($output_more_like_old_cmfinder_pl) {
		print "Can't exec \"clustalw\": No such file or directory at $bin_path/merge_motif.pl line 290.\n";
		print "\n";
		print "FATAL: Alignment file $out.gap.aln could not be opened for reading\n";
		print "Illegal division by zero at $bin_path/io.pl line 353.\n";
	    }
	    else {
		print "\t\taborting merge_motif since -skipClustalw\n";
	    }
	    return undef;
	}
	else {
	    RunCmd("$bin_path/clustalw -infile=$out.gap -outfile=$out.gap.aln");
	    RunCmd("$bin_path/sreformat stockholm $out.gap.aln > $out.gap.align");    
	    my $gap_sto = read_stockholm("$out.gap.align");    
	    my $gap_align  = $gap_sto ->seqs;
	    for my $id (keys %motif_overlap) {
		if (exists $gap_align->{$id}){
		    $gap_seq = $gap_align->{$id}->align_seq;
		    $motif_overlap{$id}{gap_seq}=$gap_seq;
		    $motif_overlap{$id}{gap_ss} = "";
		    if ($max_gap_len < length($gap_seq)){
			$max_gap_len = length($gap_seq);
		    }
		}	    
	    }
	    system("rm -f $out.gap*");
	}
    }

    $merged_ss_cons = $ss_cons1.pad_string("", $max_gap_len, '.',1).$ss_cons2;
    $merged_rf = $rf1.pad_string("", $max_gap_len, '.',1).$rf2;

    my ($align_seq,$align_ss,$score1,$score2,$score,$weight,$desc);
    for my $id (keys %motif_overlap) {

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

sub CombMotif {
    my ($cand_weight_option,$seq_file,$motifFilesRef)=@_;

    my @align_files = @$motifFilesRef;

    my @all_files = @align_files;
    my %all_stats=();
    for my $f (@align_files) {
	my $align = read_stockholm($f);  
	if ($align->weight >= $comb_min_num){
	    my %stat=RunSummarize($f);
	    $all_stats{$f} = \%stat;
	    $alignments{$f}= $align;
	}
    }

    # try merging all pairs of motifs, and see how they fit together
    my %merge_motif=();
    for my $f1 (@align_files){ 
	for my $f2 (@align_files){ 
	    if ($f1 lt $f2){
		try_merge($f1,$f2,\%merge_motif);
	    }
	}
    }

    my %processed=();
    my %merged_files=();
    while (scalar keys %merge_motif > 0) {
	if ($verbose) { print "entering for my \$id, with list:\n".join("",map {"\t$_\n"} reverse sort {$merge_motif{$a}->weight <=> $merge_motif{$b}->weight} keys %merge_motif); }
	for my $id (reverse sort {$merge_motif{$a}->weight <=> $merge_motif{$b}->weight} keys %merge_motif) { # take motifs whose combination has the biggest weight first
	    if ($verbose) { print "while keys merge_motif : id=$id\n"; }

	    my $m = $merge_motif{$id};
	    delete($merge_motif{$id});

	    if (exists $processed{$id}) {
		if ($verbose) { print "\twhile keys merge_motif NEXT : exists \$processed{$id}\n"; }
		next;
	    }
	    $processed{$id}=$m;
	    my $f1 = $m->motif1;
	    my $f2 = $m->motif2;	

	    if ($m->weight <= 0) {
		if ($verbose) { print "\twhile keys merge_motif NEXT : weight <= 0\n"; }
	    }
	    if ($m->weight > 0) { # has to be favorable
		if (exists $merged_files{$f1} || exists $merged_files{$f2}) {
		    # ??  apparently each motif can only appear in one merging?
		    if ($verbose) { print "\twhile keys merge_motif NEXT : exists \$merged_files{$f1}=".(exists $merged_files{$f1}?"true":"false")." || exists \$merged_files{$f2}=".(exists $merged_files{$f2}?"true":"false")."\n"; }
		    next;
		}

		if ($m->gap > $comb_max_gap) {
		    my $m_gap=$m->gap;
		    if ($verbose) { print "\twhile keys merge_motif NEXT : \$m->gap > \$comb_max_gap : $m_gap > $comb_max_gap\n"; }
		    next;
		}

		my $f = $id;
		my $found = 0;
		for my $tmp (@all_files){
		    if ($tmp eq $f){
			$found = 1;
			last;
		    }
		}
		if ($found) {
		    if ($verbose) { print "\twhile keys merge_motif NEXT : not found.  all_files = ".join(" ",@all_files)."\n"; }
		    next; # file has already been made
		}

		print "( near merge_motif: $id, ",$m->num_seq,"\t$f1\t$f2\t",$m->weight,"\t",$m->gap, "\t",$m->overlap," )\n";

		# this used to be a call to the merge_motif.pl script
		if ($output_more_like_old_cmfinder_pl) {
		    if ($verbose) { print "$bin_path/merge_motif.pl $SEQ $f1 $f2 $f.temp\n"; }
		}
		else {
		    if ($verbose) { print "\tcall to merge_motif $f1 $f2 $f.temp\n"; }
		}
		if (!defined($alignments{$f1}) || !defined($alignments{$f2})) {
		    die "internal error";
		}
		my $f_temp="$f.temp";
		my $new_alignment=merge_alignment($alignments{$f1},$alignments{$f2},$unaligned_seqs,$f_temp);
		if (!defined($new_alignment)) {
		    if (!$skipClustalw) { die "unexpected"; }
		    if ($output_more_like_old_cmfinder_pl) {
			if ($verbose) {
			    my $cmfile=$seq_file;
			    $cmfile =~ s/[.]motif[.]/.cm./g;
			    print "$bin_path/cmfinder      -o $f -a $f_temp $seq_file $cmfile\n";
			    print "\n";
			    print "FATAL: Alignment file $f_temp could not be opened for reading\n";
			    print "while keys merge_motif NEXT : !-s $f\n";
			}
		    }
		    else {
			if ($verbose) { print "\twhile keys merge_motif NEXT : -skipClustalw\n"; }
		    }
		    next;
		}
		write_stockholm($new_alignment,$f_temp);

		if (!RunCmfinder("$bin_path/$cmfinderBaseExe $saveTimerFlag $cmfinder_inf11Flags $cand_weight_option -o $f -a $f_temp $seq_file $dummyCmfileParamForCmfinder",$f)) {
		    # couldn't produce acceptable output
		    next;
		}
		if (!$verbose) {
		    system("rm -f $f.temp*");
		}
		if (!-s $f) {
		    die "internal error - unless exit code is 2, cmfinder_inf11 should produce a file";
		}

		# add summary stats on the merged motif
		my %stat=RunSummarize($f);
		$all_stats{$f} = \%stat;	    
		
		# check if we should reject the merged motif
		if ($stat{"Weight"} < $comb_min_num ||  # same kind of check that cmfinder does, but in general comb_min_num could be different from cmfinder's parameter
		    ($stat{"Score"} < $all_stats{$f1}->{"Score"} && $stat{"Weight"} <= $all_stats{$f1}->{"Weight"})  ||  # the score and weight of the merged motif is worse than one of the original motifs
		    ($stat{"Score"} < $all_stats{$f2}->{"Score"} && $stat{"Weight"} <= $all_stats{$f2}->{"Weight"})  ||
		    ($stat{"BP.org"} < $all_stats{$f1}->{"BP.org"} + 3 && $stat{"Weight"} <= $all_stats{$f1}->{"Weight"})  ||  # the BP.org and weight of the merged motif is worse that that of one of the original motifs, or not much better on the BP.org side (by 3), where BP.org is avg_bp in summarize.c and just counts the number of base pairs, so we're saying we need to predict at least 3 additional base pairs for it to be worth it.
		    ($stat{"BP.org"} < $all_stats{$f2}->{"BP.org"} + 3 && $stat{"Weight"} <= $all_stats{$f2}->{"Weight"})) {

		    # reject merged motif
		    if ($verbose) { print "\twhile keys merge_motif NEXT : score/weight/BP.org not good enough\n"; }
		    print "\t\tRemove $f\n";
		    print "\t\tremove-info $f: ", $stat{"Weight"}, "\t", $stat{"Score"}, "\t", $stat{"BP"}, "\n";
		    print "\t\tremove-info $f1:", $all_stats{$f1}->{"Weight"}, "\t", $all_stats{$f1}->{"Score"},"\t", $all_stats{$f1}->{"BP"}, "\n";		
		    print "\t\tremove-info $f2:", $all_stats{$f2}->{"Weight"}, "\t", $all_stats{$f2}->{"Score"},"\t", $all_stats{$f2}->{"BP"}, "\n";		
		    system("rm -f $f");
		}
		else {
		    # merged motif is good, let's keep it
		    if ($verbose) { print "\t\tset merged_files{$f1}=1\n";}
		    if ($verbose) { print "\t\tset merged_files{$f2}=1\n";}
		    $merged_files{$f1} = 1;
		    $merged_files{$f2} = 1;	
		    if ($stat{"Energy"} > $stat{"Len"} * $comb_len_energy_threshold || $stat{"Len"} > $comb_max_len) {
			# if energy per nuc is not good enough, or the length is greater than the maximum we'll allow
			if ($verbose) { print "\twhile keys merge_motif NEXT : \$stat{\"Energy\"} > \$stat{\"Len\"} * \$comb_len_energy_threshold || \$stat{\"Len\"} > \$comb_max_len\n"; }
			next;
		    }

		    # keep the new merged MSA in memory
		    my $merge_align = read_stockholm($f);
		    $alignments{$f}= $merge_align;

		    if ($motifList) {
			print MOTIFLIST "$f\n";
		    }


		    # see if we can merge this MSA with further MSAs
		    for my $f3 (@align_files){
			if (!exists $merged_files{$f3}){
			    if ($verbose) { print "\twhile keys merge_motif : try_merge in context of id=$id\n"; }
			    try_merge($f, $f3, \%merge_motif);
			}
			else {
			    if ($verbose) { print "\twhile keys merge_motif : !exists \$merged_files{$f3}, so skipping\n"; }
			}
		    }

		    # add the list of files we know about
		    push @align_files, $f;

		    # loop to find the next-best merge option, taking into account the merged motifs that we added
		    last;
		}		
	    }
	}
    }
}

sub try_merge {
    my ($f1,$f2,$merge_motif_ref)=@_;
    if (!defined($merge_motif_ref)) { die;}
    my $index = join(".", my_strcmp($f1, $f2));
    print "try_merge $f1 $f2\n";
    if (exists $$merge_motif_ref{$index}) {
	if ($verbose) { print "\texists \$\$merge_motif_ref{$index}\n"; }
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
	    
    if (!exists $alignments{$f1} || !exists $alignments{$f2}) {
	# one of the alignments didn't satify an earlier criterion ($align->weight >= $min_num)
	# so don't proceed
	if ($verbose) { print "!exists \$alignments{$f1} || !exists \$alignments{$f2}\n";}
	return;
    }
    #detect overlap
    my $align1 = $alignments{$f1}->seqs;
    my $align2 = $alignments{$f2}->seqs;
    for my $id (keys %$align1) {

	next if (!exists $align2->{$id});
	# below, we know that hit id $id is common to both alignments $f1 and $f2

	# figure out if the hits for id=$id in the two motifs overlap
	my $motif1 = $align1->{$id}; # motif1 is info on hit id=$id in align1
	my $motif2 = $align2->{$id};
	my $len1 = abs($motif1->end - $motif1->start);
	my $len2 = abs($motif2->end - $motif2->start);
	if ($motif1->start > $motif1->end) { die "unexpected";} # we only look for motifs on the forward strand
	if ($motif2->start > $motif2->end) { die "unexpected";}
	#print "$id ", $motif1->start, "-", $motif1->end, "\t", $motif2->start, "-", $motif2->end, "\n";

	if ($motif1->start > $motif2->start) {
	    if ($motif1->end < $motif2->end) {
		# overlap
		# motif1 contained within motif2
		$num_overlap += $motif1->weight * $motif2->weight;
	    }
	    else{ 
		my $overlap = $motif2->end - $motif1->start;
		#print "$id $overlap\n";
		if ($overlap > - $comb_max_gap) { # they either overlap, or have no more than $comb_max_gap nucs between them
		    if(($len1 - $overlap < 25 || $len2 - $overlap < 25) &&
		       ($overlap > 0.9 * $len1 || $overlap > 0.9 * $len2)) {
			# the motif instances almost entirely overlap (have up leave at most 25 nucs and 10% of each instance un-overlapped
			$num_overlap+=  $motif1->weight * $motif2->weight;
		    }
		    else {
			# the motif instances are near to one another or somewhat overlapping, let's record the weight for this orientation (motif2 instance coords > motif1 instance coords)
			$start2 +=  ($motif1->weight + $motif2->weight)/2;
			$start2_score += $motif1->score * $motif1->weight + $motif2->score* $motif2->weight;
			if ($overlap > 0) {
			    # the motif instances actually overlap
			    $overlap2 += $overlap *  $motif1->weight * $motif2->weight;
			}
			elsif ( - $overlap < $comb_max_gap) { # I think this is always true because of the test some lines up
			    $gap2 += - $overlap * $motif1->weight * $motif2->weight;
			}
		    }
		}
	    }
	}
	else {
	    # this 'else' case is the mirror image of the 'if' case
	    if ($motif1->end < $motif2->end) {
		my $overlap = $motif1->end - $motif2->start;
		#print "$id $overlap\n";
		if ($overlap > - $comb_max_gap){
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
			elsif ( - $overlap1 < $comb_max_gap){
			    $gap1 += - $overlap *  $motif1->weight * $motif2->weight;
			}
		    }
		}
	    }
	    else {
		#overlap
		# motif2 contained within motif1
		$num_overlap +=  $motif1->weight * $motif2->weight;
	    }
	}		
    }
    if ($verbose) { print "\t( start1,start2,num_overlap = $start1,$start2,$num_overlap )\n"; }
    if ($num_overlap > $start1 && $num_overlap > $start2) {
	# the motifs are almost identical, in that they overlap a lot
	if ($verbose) { print "\$num_overlap > \$start1 && \$num_overlap > \$start2: $num_overlap > $start1 && $num_overlap > $start2\n"; }
	return;
    }
    if ($start1_score > $start2_score) {  # is motif1 more often to the left of motif2 or vice versa?  Where's the bias?

	if ($start1 < $comb_min_overlap) {
	    # not enough sequences favor this orientation in absolute terms (compared to the constant $comb_min_overlap) for it to be worthwhile
	    if ($verbose) { print "\t\$start1 < \$comb_min_overlap : $start1 < $comb_min_overlap\n";}
	    return;
	}

	#motif1 is before motif2
	$$merge_motif_ref{$index} = 
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
    else {
	# this 'else' case mirrors the 'if' case
	if ($start2 < $comb_min_overlap) {
	    if ($verbose) { print "\t\$start2 < \$comb_min_overlap : $start2 < $comb_min_overlap\n"; }
	    return;
	}
	#motif2 is before motif1
	$$merge_motif_ref{$index} = 
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

sub my_strcmp { # for file names $s1,$s2, return ($prefix,$suffix1,$suffix2) where $prefix is the common part of there names (split by dots), $suffix1 is the unique suffix of $s1 and $suffix2 is the unique suffix of $s2
    my ($s1, $s2)= @_;
    my @t1 = split /\./, $s1;
    my @t2 = split /\./, $s2;
    my $prefix="";
    my $i;
    for( $i = 0; $i < scalar @t1 && $i < scalar @t2; $i++) {
	last if ($t1[$i] ne $t2[$i]);	
	if ( $prefix eq "") {
	    $prefix =$t1[$i];	
	}
	else{
	    $prefix ="$prefix.".$t1[$i];	
	}
    }
    my $suffix1 = join ".", @t1[$i..($#t1)];
    my $suffix2 = join ".", @t2[$i..($#t2)];    
    return ($prefix, $suffix1, $suffix2);
}	  


sub RunCmd {
    my ($cmd)=@_;
    if ($output_more_like_old_cmfinder_pl) {
	print "$cmd\n";
    }
    else {
	print "Running: $cmd\n";
    }
    my $status=system($cmd);
    if ($status!=0) {
	die "\tcommand failed: $? $!";
    }
}

sub RunCmfinder {
    my ($cmd,$output_motif_file)=@_;
    if (!defined($output_motif_file)) { die "must pass output_motif_file";}
    if ($output_more_like_old_cmfinder_pl) {
	print "$cmd\n";
    }
    else {
	print "Running: $cmd\n";
    }
    my $status=system($cmd);
    if ($status==0) {
	# made an alignment
	return 1;
    }
    else {
	if ($status==(2<<8)) {
	    # couldn't produce acceptable output
	    return 0;
	}
	else {
	    if ($useOldCmfinder) {
		print "cmfinder returned error, but I'm ignoring it because I think it's benign\n";
		unlink($output_motif_file);
		return 0;
	    }
	    else {
		die "problem running cmd: $? $!";
	    }
	}
    }
}    

sub RunSummarize {
    my ($f)=@_;
    if (0) {
	my $cmd="$bin_path/summarize $f";
	print "Running $cmd\n";
	my $summary = `$cmd`;
	if ($? != 0) {
	    die "command \"$cmd\" failed: $? $!";
	}
    }
    my $cmd2="$bin_path/cmfinder04 $saveTimerFlag --summarize $summarizeFlagsStr $f";
    print "Running $cmd2\n";
    my $summary2 = `$cmd2`;
    if ($? != 0) {
	die "command \"$cmd2\" failed: $? $!";
    }
    #if ($summary ne $summary2) {	die "summary=\"$summary\"\nsummary2=\"$summary2\"";    }

    my %stat = $summary2=~ /\s*(\S+)=(\S+)\s*/g;
    return %stat;
}

sub GetCandFiles {
    my ($tempFileListFileName)=@_;
    my @l=();
    if (!open(F,$tempFileListFileName)) {
	die "cannot open $tempFileListFileName";
    }
    while (<F>) {
	s/[\r\n]//g;
	push @l,$_;
    }
    if (!close(F)) {
	die "problem closing $tempFileListFileName";
    }
    return @l;
}


sub print_help{
    print STDERR "bin_path=$bin_path\n";
    print STDERR <<HELP;
CMFINDER [options] SEQ
Options:
    -c <number>      
     The maximum number of candidates in each sequence. Default 40. No bigger than 100.
    -m <number>      
     The minimum length of candidates. Default 30
    -M <number>      
     The maximum length of candidates. Default 100
    -f <number>      
     The fraction of the sequences expected to contain the motif. Default 0.80
    -s1 <number>     
     The max number of output single stem-loop motifs
    -s2 <number>    
     The max number of output double stem-loop motifs    
    -minspan1 <number>
     minimum span of a candidate sub-sequence in the heuristics to come up with an initial alignment for single-hairpin (h1) motifs
    -maxspan1 <number>
     like -minspan1, but maximum
    -minspan2 <number>
     like -minspan1, but for double-hairpin (h2) motifs
    -maxspan2 <number>
     like -minspan2, but maximal
    -combine         
     Combine the output motifs. Default False
    -motifList <file> 
     Produce a list of motifs generated, one motif per line.
    -o <number>
     Minimum overlap for combining motifs
    -n <number>      
     Minimum number of sequences (weighted) for combining motifs
    -emSeq <file>
     Use the sequences in this fasta file for the expectation maximization step (i.e., the C executable cmfinder), but not for the earlier steps related to finding candidate motifs.  The reason for this distinction is that it is somewhat easier to add weighting to the cmfinder program, than the various canda, candf, cands and align programs.
    -likeold         
     Behave as much as possible like the old CMfinder, e.g., passing --enone, --p56 and --degen-rand to cmfinder_inf11.  It's not possible to produce identical results to CMfinder 0.3, but these flags make it more similar.
    -fragmentary
     Pass --fragmentary for cmfinder
    -amaa            
     Pass --amaa to cmfinder (align max align accuracy)
    -useOldCmfinder  
     Run the old cmfinder executable, mainly to test whether we get different results because of this perl script, or the cmfinder_inf11 executable.
    -skipClustalw    
     Do not run clustalw, like older installations lacking this program.
    -justGetCmfinderCommand    
     Print the command to run for the cmfinder executable, with appropriate partial flags.  This can be used to realign an existing .sto file, for example.
    -copyCmfinderRunsFromLog <log-file> 
     For debugging.  Reads a log file that contains cmfinder commands, and re-runs them with new CMfinder.
    -commaSepEmFlags x<flags>
     List of flags and arguments to pass to the EM-step cmfinder exe.  There's an 'x' at the beginning of the flags, so that perl doesn't interpret the flags as flags for it.  It's comma-separated where on the command line it would be space separated.  I think commas are safe, and mean that I don't have to worry about quoting stuff.  e.g., -commaSepEmFlags x--fragmentary,--filter-non-frag,--filter-non-frag-pad,10 would pass this to the cmfinder program: "--fragmentary --filter-non-frag --filter-non-frag-pad 10", i.e., just replace commas with spaces.
    -commaSepSummarizeFlags x<flags>
     Flags to pass to the --summarize command.  Same syntax as for --commaSepEmFlags
    -commaSepCandfFlags x<flags>
     Flags to pass to the candf command.  Same syntax as for --commaSepEmFlags
    -minCandScoreInFinal <number>    
     Pass --min-cand-score-in-final <number> to cmfinder.  WARNING: there's a difference between using this flag (where even intermediate motifs will avoid these hits) and taking the low-scoring instances out of the final alignments (which might be combinations of motifs in which the sequence would have been lower-scoring).
    -filterNonFrag
     Pass --filter-non-frag to cmfinder
    -columnOnlyBasePairProbs
     Pass --column-only-base-pair-probs to cmfinder
    -saveTimer <file>
     create tab-delimited <file> containing timing stats on various sub-processes of this script.  the first tab-delimited field is the description of the sub-process, the second field is the total CPU time (user+sys) and the third field is the wall-clock time.  Sub-processes can occur in multiple lines if they are run multiple timers, so the caller should add them.  Due to my laziness, the time of the clustalw program (if used) is not counted.
    -cpu <num>
     use <num> CPUs for functionality that can use multi-CPUs (currently only the internal cmsearch commands in cmfinder04)
    -allCpus
     equivalent to -cpu X , where X is the number of available processors.
    -candsParallel
     run the two cands jobs in parallel, even if -cpu 1
    -outFileSuffix <string>
     add <string> to the output file names.  this is useful if you want to run the script in multiple ways in the same directory.
    -h               
     Show this list
    -version
     Show package version
HELP
}
