#! /usr/local/bin/perl -w

$file = shift @ARGV;

open IN, $file or die "Can't open file $file\n";
$query = "";
$subject = "";
while(<IN>) {
    chomp;
    if (/^Query=\s*(\S+)\s*/) {
	$query = $1;	
	next;
    }
    elsif (/^>(\S+)\s*$/) {
	$subject=$1;
    }
    elsif (/^Query:\s*(\d+)\s+-\s+(\d+)$/){
	$qstart = $1;
	$qend = $2;
    }
    elsif (/Expect = (\S+),/){
	$e_val = $1;
    }
    elsif (/^Sbjct:\s*(\d+)\s+-\s+(\d+)$/){
	$sstart = $1;
	$send = $2;
	if ($query ne $subject && $e_val < 0.1) {
	    print "$query :\t $qstart - $qend\t $e_val\n";
	    print "$subject :\t $sstart - $send\n\n";
	}
    }  
}
