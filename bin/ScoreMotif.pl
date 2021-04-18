use strict;
use Getopt::Long;
use FindBin qw($Bin);
use DB_File;
use lib $Bin;
use Stockholm; # Stockholm.pm should be in the same directory as this script

my $dataDir="$Bin/../data";
my $rnaphyloDataDir=undef;
my ($hmmpairDbFile,$rnaphyloDbFile,$pfoldRateFile);
my $maxNonCanonicalPairFreq=0.05;
my $usePartitionFunction=1;
my $numFlankingPolyN=200;
my $positionsToIgnoreFileName;
my $hmmpairExe="$Bin/hmmpair";
my $rnaphyloExe="$Bin/RNAPhylo";
my $disableHmmpair=0;
my $disablePscore=0;
my $disableCombScore=0;
my $tempBase="/tmp/$$";
my $pfoldPath="$Bin";
my $help;

if (!GetOptions(
	 "dataDir=s" => \$dataDir,
	 "rnaphyloDataDir=s" => \$rnaphyloDataDir,
	 "hmmpairDb=s" => \$hmmpairDbFile,
	 "rnaphyloDb=s" => \$rnaphyloDbFile,
	 "terminators=s" => \$positionsToIgnoreFileName,
	 "pfoldRateFile=s" => \$pfoldRateFile,
	 "tempBase=s" => \$tempBase,
	 "pfoldPath=s" => \$pfoldPath,
	 "hmmpairExe=s" => \$hmmpairExe,
	 "rnaphyloExe=s" => \$rnaphyloExe,
	 "disableHmmpair" => \$disableHmmpair,
	 "disablernaphylo" => \$disablePscore,
	 "disablePscore" => \$disablePscore,
	 "h" => \$help
    )) {
    die "problem with command line";
}

if ($help) {
    print "usage\nperl $0 [flags] <sto-file>\n";
    print "flags are:\n";
    print "-dataDir <dir> : directory for data files\n";
    print "-rnaphyloDataDir <dir> : directory for RNAphylo data files, which just needs to contain pscore.rate\n";
    print "-hmmpairDb <file> : the file containing a large number of hmmpair scores in which ScoreMotif.pl should calculate the rank.  The file is made in Perl using the DB_FILE module with the \$DB_TREE method, with \$DB_TREE->\{'compare'\} set to the <=> function.\n";
    print "-rnaphyloDb <file> : as for -hmmpairDb, the file containing scores for RNAphylo\n";
    print "-terminators <file> : comma-separated file containing positions of predicted transcription terminators, whose pairs hmmpair should ignore.  See user-manual.pdf for details.\n";
    print "-pfoldRateFile <file> : the rate file for inferring phylogenetic code.\n";
    print "-tempBase <path> : a prefix of paths to use to creating temporary files. By default it is /tmp/PID , where PID is the process ID of this script.\n";
    print "-pfoldPath <dir> : directory containing executables that are part of a version of Pfold that inferred phylogenetic trees.\n";
    print "-hmmpairExe <file> : the executable hmmpair file\n";
    print "-rnaphyloExe <file> : the executable RNAPhylo file\n";
    print "did help\n";
    exit(0);
}

my $stoFileName=shift;
if (!$stoFileName) {
    die "you must give a Stockholm alignment file on the command line";
}

if ($disableHmmpair || $disablePscore) {
    $disableCombScore=1;
}

if (!defined($tempBase)) {
    die "you need to use -tempBase <base-path-for-files> , since temporary files are required for the use of RNAphylo";
}

if (!defined($hmmpairDbFile)) {
    $hmmpairDbFile="$dataDir/runid550-scores.score-btree.hmmpair";
}
if (!defined($rnaphyloDbFile)) {
    $rnaphyloDbFile="$dataDir/runid550-scores.score-btree.pscore";
}
if (!defined($pfoldRateFile)) {
    $pfoldRateFile="$dataDir/pscore.rate";
}
if (!defined($rnaphyloDataDir)) {
    $rnaphyloDataDir=$dataDir;
}

my $hmmpairScore=undef;
if (!$disableHmmpair) {
    my $fragmentaryPolicy="f"; # just hardcode this
    my $uniformDistributionOfProfileHmmStartsAndEnds=0; # default, old thingy
    my $positionsToIgnoreFileName_cmdline=defined($positionsToIgnoreFileName) ? $positionsToIgnoreFileName : "NULL";
    my $cmd="$hmmpairExe $stoFileName $maxNonCanonicalPairFreq $fragmentaryPolicy $numFlankingPolyN $uniformDistributionOfProfileHmmStartsAndEnds $positionsToIgnoreFileName_cmdline";
    print "Running $cmd\n";
    if (!open(OUT,"$cmd |")) {
	die "problem opening $cmd: $? $!";
    }
    while (<OUT>) {
	s/[\r\n]//g;
	if (/^SCORE ([^ ]+) score=([-+.e0-9]+)/) {
	    my $statName=$1;
	    my $score=$2;
	    
	    my %statNameToAbbrev=("bestPair","b","bestPairPf","bpf","pairSum","s","pairSumPf","spf");
	    my $statAbbrev=$statNameToAbbrev{$statName};
	    if (!$statAbbrev) {
		die "cannot find abbreviation of statistic name $statName from line $_";
	    }
	    if ($statAbbrev eq "spf") {
		$hmmpairScore=$score;
	    }
	    print "hmmpair score $statName ($statAbbrev): $score\n";
	}
    }
    if (!close(OUT)) {
	die "problem reading from $cmd: $? $!";
    }
}
print "hmmpair score: $hmmpairScore\n";

my $rnaphyloScore=undef;
if (!$disablePscore) {
    my $quietlyRejectMotif=0;

    # convert Stockholm alignment to Fasta alignment (with dashes), so that we can use the fasta2col program.  Also check if there's any predicted pairs
    my $sto=new Stockholm;
    $sto->Load($stoFileName);

    # check for colons in hitIds, which creates problems with the phylo tree
    for my $hit ($sto->ListHits()) {
	if ($hit =~ /[():,]/) {
	    die "the sequence name \"$hit\" in the input alignment ($stoFileName) has a character that will create problems with the phylogenetic tree.  These problematic characters are: parentheses, colon and comma.  The characters create problems with the Newick format.  Please remove the problematic characters from the sequence names.";
	}
    }
    
    if ($sto->{gc}{SS_cons} =~ /^[-._:;~]+$/) {
	# there are no pairs at all -- can't be good evidence of pairing
	$rnaphyloScore=-100; # very bad score
    }
    else {
	my $fastaFileName="$tempBase.fasta";
	$sto->SaveAlignedFasta($fastaFileName);
	    
	my $pfoldOutFileName="$tempBase.pfoldout";
	unlink($pfoldOutFileName);
	my $colFileName2="$tempBase.col2";
	$sto->SavePfoldCol($colFileName2);
	my $findphylFileName="$tempBase.findphyl";
	RunCmd("$pfoldPath/findphyl $pfoldRateFile < $colFileName2 > $findphylFileName");
	RunCmd("$pfoldPath/mltree $pfoldRateFile < $findphylFileName > $pfoldOutFileName");
# was:	    RunCmd("cat $fastaFileName | $pfoldPath/fasta2col | sed 's/arbitrary/RNA/g' | $pfoldPath/findphyl $pfoldRateFile | $pfoldPath/mltree $pfoldRateFile > $pfoldOutFileName");

	my $treeFileName="$tempBase.tree";
	PfoldOutToPscoreTree($pfoldOutFileName,$treeFileName);

	my $rnaphyloFlags="--fragmentary --partition --partition-close-to-fuzzy-limit 3 --suspicious-degen-nucs 2";
	my ($score,$rnaScore)=RunPscore($treeFileName,$stoFileName,$rnaphyloFlags);
	$rnaphyloScore=$rnaScore;

	unlink($fastaFileName,$pfoldOutFileName,$colFileName2,$findphylFileName,$treeFileName);
    }
}
print "RNAphylo score: $rnaphyloScore\n";

my $combinedScore=undef;
my $pscoreRank=undef;
my $hmmpairRank=undef;
if (!$disableCombScore) {
    my $pscoreDb=OpenCombinedScoreDbToRead($rnaphyloDbFile);
    my $hmmpairDb=OpenCombinedScoreDbToRead($hmmpairDbFile);
    $pscoreRank=GetCombinedScoreRank($pscoreDb,$rnaphyloScore);
    $hmmpairRank=GetCombinedScoreRank($hmmpairDb,$hmmpairScore);
    $combinedScore=$pscoreRank<$hmmpairRank ? -$pscoreRank : -$hmmpairRank;
    CloseCombinedDb($pscoreDb);
    CloseCombinedDb($hmmpairDb);
}
print "combined score: $combinedScore (pscoreRank=$pscoreRank,hmmpairRank=$hmmpairRank)\n";

sub PfoldOutToPscoreTree {
    my ($pfoldOutFileName,$treeFileName)=@_;
    
    if (!open(PFOLD,$pfoldOutFileName)) {
	die "cannot open $pfoldOutFileName";
    }
    if (!open(TREE,">$treeFileName")) {
	die "cannot open $treeFileName";
    }

    my %tree=();
    my $state="start";
    while (<PFOLD>) {
	s/[\r\n]//g;
	if ($state eq "start") {
	    if (/^; +TYPE +TREE/) {
		$state="tree header";
	    }
	}
	if ($state eq "tree") {
	    if (/^; [*]/) {
		$state="rest of file";
	    }
	    else {
		if (/^ +N +([0-9]+) +([^ ]+) +([-+.e0-9]+) +([.0-9]+) +([.0-9]+)$/) {
		    my ($nodeNum,$nodeName,$branchLen,$child,$sibling)=($1,$2,$3,$4,$5);
		    $tree{$nodeNum}={name=>$nodeName,branchLen=>$branchLen,leftChild=>$child,nextSibling=>$sibling};
		}
		else {
		    die "in state \"tree\" could not parse \"$_\"";
		}
	    }
	}
	if ($state eq "tree header") {
	    if (/^; ----------/) {
		$state="tree";
	    }
	}
    }

    my $root="1";

    print TREE DumpPscoreTree_Recurse(\%tree,$root,0);
    print TREE ";\n";

    close(TREE);
    close(PFOLD);
}

sub RunPscore {
    my ($treeFileName,$stoFileName,$flagListStr)=@_;

    my $cmd="$rnaphyloExe --ignore-all-gap $flagListStr -t $treeFileName -s $rnaphyloDataDir/single_model -p $rnaphyloDataDir/pair_model -g $rnaphyloDataDir/scfg $stoFileName";
    print "Running $cmd\n";
    if (!open(PSCORE,"$cmd |")) {
	die "cannot open $cmd: $? $!";
    }
    
    # parse score
    my $score=undef;
    my $rnaScore=undef;
    while (<PSCORE>) {
	s/[\r\n]//g;
	if (/^Total pair posterior ([-+e.0-9]+)$/) {
	    $score=$1;
	}
	if (/^Total pair posterior [-]*nan$/) {
	    $score=-1e10;
	}
	if (/^Total RNA posterior ([-+e.0-9]+)$/) {
	    $rnaScore=$1;
	}
	if (/^Total RNA posterior inf$/) {
	    $rnaScore=1e10;
	}
	if (/^Total RNA posterior [-]*nan$/) {
	    $rnaScore=-1e10;
	}
	if (/BEGINPSCORE/) {
	    $score=undef;
	}
    }
    if (!close(PSCORE)) {
	die "problem reading from $cmd: $! $?";
    }
    if (!defined($score)) {
	die "could not parse score from $cmd";
    }
    if (!defined($rnaScore)) {
	if (defined($score) && $score==0.00) {
	    $rnaScore=0;
	}
	else {
	    die "could not parse RNA posterior score from $cmd";
	}
    }
    return ($score,$rnaScore);
}

sub NumericalCompare {
    my ($x,$y)=@_;
    return $x <=> $y;
}

sub DumpPscoreTree_Recurse {
    my ($tree,$node,$depth)=@_;
    my $x=$$tree{$node};
    my $linebreaks=1;
    my $s="";
    if ($x->{leftChild} eq ".") {
	# leaf
	if ($linebreaks) { $s .= " " x $depth; }
	$s .= "$x->{name}:$x->{branchLen}";
    }
    else {
	# internal node
	if ($linebreaks) { $s .= " " x $depth; }
	$s .= "(";
	my $child=$x->{leftChild};
	while (1) {
	    if ($child ne $x->{leftChild}) {
		$s .= ",";
	    }
	    if ($linebreaks) { $s .= "\n"; }
	    $s .= DumpPscoreTree_Recurse($tree,$child,$depth+1);
	    $child=$$tree{$child}->{nextSibling};
	    if ($child eq ".") {
		last;
	    }
	}
	$s .= ")";
	$s .= ":$x->{branchLen}";
    }
    return $s;
}

sub OpenCombinedScoreDbToRead {
    my ($dbFileName)=@_;
    $DB_BTREE->{'compare'}=\&NumericalCompare;
  
    my %tree;
    my $db=tie %tree, "DB_File", $dbFileName, O_RDONLY, 0666, $DB_BTREE
      or die "cannot open db file \"$dbFileName\"";

    return {db=>$db,tree=>\%tree};
}
sub CloseCombinedDb {
    my ($db)=@_;
    undef $db->{db};
    untie %{$db->{tree}};
}
sub GetCombinedScoreRank {
    my ($db,$score)=@_;

    my $orig_score=$score;
    my $value=-1;
    $db->{db}->seq($score,$value,R_CURSOR);
    if ($value==-1) {
	my ($highest_score,$highest_value);
	$db->{db}->seq($highest_score,$highest_value,R_LAST);
	if ($orig_score>$highest_score) {
	    $value=$highest_value;
	}
	else {
	    die "couldn't look up score=$score in db";
	}
    }
    else {
	if ($score!=$orig_score) {
	    $value++; # this is for compatibility with earlier scripts
	}
    }
    return $value;
}

sub RunCmd {
    my ($cmd)=@_;
    print "Running $cmd\n";
    my $status=system $cmd;
    if ($status != 0) {
	die "command failed: $? $! ($cmd)";
    }
}
