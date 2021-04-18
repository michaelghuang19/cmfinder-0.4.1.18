#!/usr/bin/perl -w


$path= $ENV{CMfinder};
$bin_path = "$path/bin";
do "$bin_path/io.pl";
use Getopt::Long;
$no_header=0;
$sort_id = "Score";
$dump_file="temp";
GetOptions(
	   "h" => \$no_header,
	   "s=s" => \$sort_id,
	   "o=s"=>\$dump_file,
	   );

sub match_file{  
    my @files = ();
    my $f = shift @_;
    while(<$f*>) {
	push @files, $_;
    }
    return @files;
}

%sort_fields=();
%all_stat=();
%max=();
$max_file=30;


sub collect_summary{
    open OUT, ">$dump_file" or die "Can't write to $dump_file\n";    
    foreach $f (keys %all_files){
	if (length($f) > $max_file){
	    $max_file = length($f);
	}
	$result = `$bin_path/summarize $f`;
	if ($result eq ""){
	    print STDERR "Error with $f";
	}
	else{
	    print OUT "file=$f\t $result";
	}
    }
    close OUT;
}

if (! -e $dump_file){
    @files = @ARGV;
    %all_files=();
    foreach $f (@files){
	@match = match_file($f);
	foreach $m (@match){
	    $all_files{$m} = 1;
	}
    }
    collect_summary();
}
open IN, $dump_file or die "Can't open $dump_file for reading";
while(<IN>){
    next if /^\s*$/;
    $result = $_;    
    @fields = $result=~ /\s*(\S+)=\S+\s*/g;
    my %stat = $result=~ /\s*(\S+)=(\S+)\s*/g;
    $f = $stat{"file"};
    $all_stat{$f}=\%stat;
    $sort_fields{$f}= $stat{$sort_id};
    foreach $field (@fields){
	$l = length($stat{$field}); 
	if (!exists $max{$field} || $l > $max{$field}) {
	    $max{$field} = $l;
	}
    }    
}

foreach $field (@fields){
    if (length($field) > $max{$field}){
	$max{$field} = length($field);
    }
}


if (!$no_header){
    printf "%".$max_file."s,", "file";
    foreach $field (@fields){
	next if ($field =~ /^file/);
	printf "%".$max{$field}."s,", $field;
    }
    print "\n";
}

foreach $f (sort {$sort_fields{$a} <=> $sort_fields{$b}} keys %sort_fields){
    printf "%".$max_file."s,", $f;
    $temp = $all_stat{$f};
    foreach $field (@fields){
	next if ($field =~ /^file/);
	printf "%".$max{$field}."s,", $temp->{$field};
    }   
    print "\n";
}


