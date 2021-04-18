#!/usr/bin/perl -w

use Getopt::Long qw(:config no_ignore_case); 
$path= $ENV{CMfinder};

#default parameters

$CAND=40;
$MAXSPAN1=100;
$MINSPAN1=30;
$MAXSPAN2=100;
$MINSPAN2=40;
$CLUSTER = 3;
$FRACTION=0.8;
$SINGLE = 5;
$DOUBLE = 5;
$verbose = 0;
$help = 0;
$COMBINE = 0;
$ANCHOR = 0;
$HMM = 0;
$BANDED = 0;
$DEFAULT = 0;
$weight_option="";

if (!GetOptions(
	   "h" => \$help,
	   "v" => \$verbose,
	   "w=s" => \$weight_option,
	   "c=i" => \$CAND,	
	   "m1=i" => \$MINSPAN1,
	   "M1=i" => \$MAXSPAN1,
	   "m2=i" => \$MINSPAN2,
	   "M2=i" => \$MAXSPAN2,
	   "cl=i" => \$CLUSTER,	
	   "f=f" => \$FRACTION,
	   "s1=i" => \$SINGLE,
	   "s2=i" => \$DOUBLE,
	   "hmm" => \$HMM,
	   "hbanded" => \$BANDED,
	   "anchor" => \$ANCHOR,	
	   "combine" =>\$COMBINE,
	   "def" =>\$DEFAULT, 	
	   )
    ){
    print STDERR "Invalid options\n";
    print_help();
    exit(1);
}

if ($help) {
    print_help();
    exit(0);
}

if (!defined($path)) {
    die "you must define the environment variable \"CMfinder\" to be the path of the CMfinder directory in this distribution.  (This directory will have a child 'bin' containing this and other perl scripts.)";
}
$bin_path = "$path/bin";
$data_path = "$path/data";

if (scalar @ARGV==0){
    print STDERR "No sequence file is specfied\n";
    print_help();
    exit(1);
}

$SEQ= shift @ARGV;

$verbose_flag=$verbose ? "-v" : "";

if ($DEFAULT){
    #Set cluster numbers based on sequence length
    $tmp = `$bin_path/count_seq $SEQ`;
    if ($tmp =~ /(\d+)\s+(\d+)/){
	$seq_num =$1;
	$seq_len =$2;
	if ($seq_len < 500){
	    $CLUSTER = 1;
	}
	elsif ($seq_len < 2000){
	    $CLUSTER  = 3;
	}
	elsif ($seq_len < 10000){
	    $CLUSTER = 4;
	}
	else{
	    $CLUSTER = 5;
	}
	if ($seq_num > 4){
	    $FRACTION = 5/ $seq_num;
	    if ($FRACTION < 0.1) {
	       $FRACTION  = 0.1;
	    }
	}
    }
    $ANCHOR = 1;
}


if ($ANCHOR){
    RunCmd("$bin_path/align -m $data_path/align_matrix $SEQ > $SEQ.match");
    if($CLUSTER > 0){
	RunCmd("$bin_path/cluster $SEQ.match $CLUSTER $SEQ.cluster");
    }    
}
@clusters = glob("$SEQ.cluster*");

if (scalar @clusters > 0){
    foreach $c (@clusters){
	if ($SINGLE){
	    $cand = "$c.cand.h1";
	    $cand =~ s/\.cluster//;
	    $cmd1 = "$bin_path/candf -c $CAND -o $cand -r $c -M $MAXSPAN1 -m $MINSPAN1 -s 1 -S 1 $SEQ";
	    $cmd2 = "$bin_path/cands  -r $c -n $SINGLE -f $FRACTION -m $SEQ.match $SEQ $cand";
	    RunCmd($cmd1);
	    RunCmd($cmd2);
	}
	if ($DOUBLE){
	    $cand = "$c.cand.h2";
	    $cand =~ s/\.cluster//;
	    $cmd1 = "$bin_path/candf -c $CAND -o $cand -r $c -M $MAXSPAN2 -m $MINSPAN2 -s 2 -S 2 $SEQ";
	    $cmd2 = "$bin_path/cands  -r $c -n $DOUBLE -f $FRACTION -m $SEQ.match $SEQ $cand";
	    RunCmd($cmd1);
	    RunCmd($cmd2);
	}
    }
}
else{
    if ($SINGLE){
	$cmd1="$bin_path/candf -c $CAND -o $SEQ.cand.h1 -M $MAXSPAN1 -m $MINSPAN1 -s 1 -S 1 $SEQ";
	$cmd2="$bin_path/cands -n $SINGLE -f $FRACTION $SEQ $SEQ.cand.h1";
	if ($verbose) {
	    print "$cmd1\n";
	}
	RunCmd($cmd1);
	RunCmd($cmd2);
    }
    if ($DOUBLE){
	$cmd1="$bin_path/candf -c $CAND -o $SEQ.cand.h2 -M $MAXSPAN2 -m $MINSPAN2 -s 2 -S 2 $SEQ";
	$cmd2="$bin_path/cands -n $DOUBLE -f $FRACTION $SEQ $SEQ.cand.h2";    
	RunCmd($cmd1);
	RunCmd($cmd2);
    }
}

$hmm_option="";
if ($HMM){
    $hmm_option = "--hmm";
}
$banded_option="";
if ($BANDED){
    $banded_option= "--hbanded";
}
if ($weight_option ne ""){
    $weight_option = "-w $weight_option";
}

@cands = glob("$SEQ.*cand.h*");
foreach $cand (@cands){   
    if ($cand =~ /$SEQ(\.\d)?\.cand.*\_\d+$/){
	$pref = $1;
	$align =$cand;
	$align =~ s/cand/align/;
	RunCmd("$bin_path/canda $cand  $SEQ $align");
	$motif = $cand;
	$motif =~ s/cand/motif/;
	$cm = $cand;
	$cm =~ s/cand/cm/;
	$cluster_option="";
	if ($ANCHOR){
	    $cluster = "$SEQ.cluster";
	    if (defined $pref) {
		$cluster = "$cluster"."$pref";
	    }
	    $cluster_option = "-r $cluster";
	}
	$cmd = "$bin_path/cmfinder $weight_option $hmm_option $banded_option $cluster_option -o $motif -a $align $SEQ $cm"; 
	if ($verbose){
	    print $cmd, "\n";
	}
	$status = system($cmd);
	if ($status != 0){
	    system("rm $motif");
	    system("rm $cm");
	    print "note: cmd failed: \"$cmd\"\n";
	}
    }
}

if ($DEFAULT){
    $COMBINE = 1;
}

if ($COMBINE){
    if ($ANCHOR){
	foreach $c (@clusters){
	    $prefix = $c;
	    $prefix =~ s/\.cluster//;
	    RunCmd("perl $bin_path/comb_motif.pl $verbose_flag $hmm_option $banded_option $weight_option -r $c $SEQ $prefix.motif");
	}	
    }
    else{
	RunCmd("perl $bin_path/comb_motif.pl $verbose_flag $hmm_option $banded_option $weight_option $SEQ $SEQ.motif");
    }
}

if (! $verbose){
    system("rm -f $SEQ.*cand*");
    system("rm -f $SEQ.*align*");
    system("rm -f $SEQ.*-*.cm");
    if ($ANCHOR) {
	system("rm -f $SEQ.lalign");
	system("rm -f $SEQ.cluster*");
	system("rm -f $SEQ.match");	
    }
}

sub RunCmd {
    my ($cmd)=@_;
    if ($verbose){
	print("$cmd\n");
    }
    my $status=system($cmd);
    if ($status!=0) {
	die "problem with \"$cmd\": $status";
    }
}

sub print_help{
    print STDERR <<HELP;
CMFINDER [options] SEQ
Options:
    -b               Do not use BLAST search to locate anchors
    -c <number>      The maximum number of candidates in each sequence. Default 40. No bigger than 100.
    -m <number>      The minimum length of candidates. Default 30
    -M <number>      The maximum length of candidates. Default 100
    -f <number>      The fraction of the sequences expected to contain the motif. Default 0.80
    -s1 <number>     The max number of output single stem-loop motifs
    -s2 <number>     The max number of output double stem-loop motifs    
    -combine         Combine the output motifs. Default False
    -hmm             Apply HMM filter for speed up. Default false.
    -anchor [FASTA|BLAST|NONE]      Methods to compute the anchors    
    -h               Show this list
HELP
}
