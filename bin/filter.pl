#! /usr/local/bin/perl -w

use Getopt::Long;
do "/projects/bio6/yzizhen/rna/bin/io.pl";

$lt = -1000;
$ut = 1000;

GetOptions(
	   "w" => \$do_weight,
           "s" => \$do_score,	   
	   "lt=f" => \$lt,
	   "ut=f" => \$ut,
       );

$file_name = shift @ARGV;
$out_file = shift @ARGV;

sub filter{
    my $lt = shift @_;
    my $ut = shift @_;
    my $alignment = shift @_;    
    
    my $temp = ($alignment->seqs);
    my %seqs = %$temp;
    $temp = ($alignment->flags);
    my %flags  = %$temp;
    if ($do_weight && ! exists $flags{"WGT"}) {
	die "No weight info is available in the alignment ";
    }
    if ($do_score && ! exists $flags{"SC"}) {
	die "No score info is available in the alignment ";
    }
    foreach $acc (keys %seqs) {
	if ($do_weight) {
	    if ($seqs{$acc}->weight < $lt || $seqs{$acc}->weight > $ut) {
		delete $seqs{$acc};
	    }
	}
	if ($do_score) {
	    if ($seqs{$acc}->score < $lt || $seqs{$acc}->score > $ut) {
		delete $seqs{$acc};
	    }
	}
	
    }
    $alignment->seqs(\%seqs);
    return ($alignment);
}

$alignment = read_stockholm($file_name);
$alignment = filter($lt, $ut, $alignment);
write_stockholm($alignment, $out_file);
#`$bin_path/sreformat --mingap stockholm $out_file  > $out_file.tmp`;
#`mv $out_file.tmp $out_file`;
