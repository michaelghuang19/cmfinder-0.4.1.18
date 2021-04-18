#!/usr/bin/perl -w

use Getopt::Long;

$rfam_dir= "";
$do_ncbi= 0;
$do_ranking=0;
$do_RNAz=0;
$do_weight = 0;
$match_dir="";
$dir = ".";
GetOptions(
	   "dir=s" => \$dir,
	   "weight" => \$do_weight,
           "rfam=s" => \$rfam_dir,
           "ncbi" => \$do_ncbi,
	   "rank" => \$do_ranking,
	   "match=s"=> \$match_dir,
	   "RNAz" => \$do_RNAz,	  
);
$file_pattern = shift @ARGV;
$out_file = shift @ARGV;

opendir(DIR, $dir) || die "can't open dir $dir";
@files =grep /$file_pattern/, readdir (DIR);
closedir(DIR);

if ($rfam_dir ne ""){
    mkdir($rfam_dir);
}

sub find_species{
    my $file = shift @_;
    my $result = `grep DE $file`;
    my @acc = ($result =~ /\#=GS\s+(\S+)\//g);
    my @sc = ($result =~ /(\d+\.\d+)$/g);
    my %temp=();
    my $total = 0;
    foreach $a (@acc){
	if (!exists $temp{$a}){	    
	    $temp{$a}=1;
	}
	else{
	    $temp{$a}++;
	}
	$total ++;	
    }
    $num_species = scalar keys %temp;
    return ($num_species);
}


sub rank{
    ($infile, $outfile)= @_;
    open R, ">$infile.r";
    print R<<EOF;
test.data <-  read.table("$infile", header=T, sep=",", na.strings="",fill=T);
test.data[is.na(test.data[,"Seq_id"]), "Seq_id"]= 1;
test.data[, "Rank.index"] = (test.data[, "Conserved_pos"] + 0.2) * (test.data[, "BP"]/test.data[,"Seq_id"] - test.data[,"Len"] * 0.1);
test.data[test.data[,"Rank.index"] < 0 ,"Rank.index"] = 0;
test.data[,"Rank.index"]= sqrt(test.data[,"Rank.index"]) *(test.data[, "Species"] - 2) *(1+ log(test.data[,"Num"]/test.data[,"Species"]));
ord <- order(test.data[,"Rank.index"], decreasing=T);
format(test.data, digits=3);
write.table(test.data[ord,], file = "$outfile", row.names=F, sep=",", quote=F,na="")
q();
EOF
    close R;
    `/homes/gws/yzizhen/R/R-2.5.0/bin/R --no-save < $infile.r > /dev/null`;
};

open OUT, ">$out_file" or die "Can't open $out_file for writing";
if ($do_ncbi){
    $rna_file = "/projects/bio5/yzizhen/NCBI_rnt/NCBI_rna";
    %rna_location=();
    open IN, $rna_file or die;
    while(<IN>){
	chomp;
	my @fields=split /\t/;
	if (!exists $rna_location{$fields[0]} ){
	    my @temp=(\@fields);
	    $rna_location{$fields[0]}= \@temp;
	}
	else{
	    $temp = $rna_location{$fields[0]};
	    push @$temp, \@fields;
	}
    }
}

$is_header = 1;
foreach $f (@files){   
    $f = "$dir/$f";
    print STDERR "$f\n";
    @paths = split "/", $f;
    $motif = $paths[$#paths];    
    if ($motif  =~ /(\d+)\.\S*motif/){
	$prefix = $1;
    }    
    #Get the features.    
    if ($do_weight){	
	$result = `summarize -w $f`;
    }
    else{
	$result = `summarize $f`;
    }
    #Invalid motif.
    next if ($result =~ /^\s*$/); 
    my %stat = $result=~ /\s*(\S+)=(\S+)\s*/g;
    my $n_species = find_species($f);
    if ($n_species == 0){
         $n_species = $stat{"Num"};
    }

    $stat{"Species"}= $n_species;
    @fields = sort keys %stat;
    if ($is_header){
	print OUT "motif,", (join ",",  @fields);	    
	if ($do_RNAz){
	    print OUT ",RNAz";
	}
	if ($rfam_dir ne ""){
	    print OUT ",Rfam";
	}
	if ($do_ncbi){
	    print OUT ",NCBI";
	}
	print OUT "\n";
	$is_header=0;
    }
    print OUT "$f";
    foreach $k (@fields){
	print OUT ",",$stat{$k};
    }	
    if ($do_RNAz){
	$RNAz_score  = `RNAz_score_cmfinder.pl $f`;
        $RNAz_score =~ s/\s//g;
	print ",$RNAz_score";
    }
    if ($rfam_dir ne ""){
	#Compare with Rfam	
	$rfam_hits="";
	$check_rfam=1;
	if ($match_dir){
	    $match_file="$match_dir/$prefix.match";            
	    if (! -e $match_file){
                $match_file="$match_dir/$f.match";            
	        if (! -e $match_file){
		     $check_rfam=0;
                }
	    }
	}
	if ($check_rfam){	    
	    $rfam_hits = `motif.rfam.overlap.pl $f $rfam_dir/$f.match`;
	    unlink("$f.blast");
	    $rfam_hits =~ s/\n//g;
	}
	print OUT ",$rfam_hits";
    }
    if ($do_ncbi){
	#Compare with NCBI	
	$rna_anno="";
	open IN, $f or die "Can't open file $f";
	
	while(<IN>){
	    if ( /#=GS\s+(\S+)\/(\d+)-(\d+) DE\s+(\d+)\.\.(\d+)/  ) {
		 $acc = $1;
		 $seq_start =$2;
		 $seq_end = $3;
		 $mstart = $4;
		 $mend = $5;		 
		 $strand = $seq_start < $seq_end;
		 if ($strand) {
		     $mstart = $mstart + $seq_start - 1;
		     $mend = $mend + $seq_start - 1;
		 }
		 else{
		     ($mend, $mstart) =  ($seq_start - $mstart + 1,$seq_start - $mend + 1);
		 }
		 $rna= $rna_location{$acc};
		 foreach $t (@$rna){
		     next if ($t->[2] < $mstart);
		     last if ($t->[1] > $mend);
		     $overlap_start = $mstart < $t->[1] ? $t->[1]: $mstart;
		     $overlap_end = $mend > $t->[2] ? $t->[2]: $mend;
		     $overlap = $overlap_end - $overlap_start;
		     $rna_anno=$rna_anno."\t $acc:".$t->[4].":".$overlap;
		 }
	     }
	}
	$rna_anno =~ s/(,|:)/ /g;
	print OUT ",", $rna_anno;
    }    
    print OUT "\n";
}
close OUT;

if ($do_ranking){
    rank($out_file, $out_file);
}
