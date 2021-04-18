#!/usr/bin/perl -w

use File::Copy;
use Cwd;

$out_motif_dir = shift @ARGV;
$rank_file = shift @ARGV;

$old_match_dir="";
if (scalar @ARGV >0 ){
    $old_match_dir = shift @ARGV;
    $new_match_dir = shift @ARGV;    
    if (! -e $new_match_dir){
	mkdir("$new_match_dir");
    }
}

if (! -e $out_motif_dir){
    mkdir("$out_motif_dir");
}

$threshold = 3;

open IN, $rank_file or die "Can't open $rank_file";
$is_header=1;
%motif_features=();
%header=();
while(<IN>){    
   s/\s+$//g;    
   next if /^\s*$/;
   @fields=split /,/;
   if ($is_header){
       for($i=0; $i <= $#fields;$i++){
	   $header{$fields[$i]} = $i;
       }
       #print $_,"\n";
       $is_header=0;
   }
   else{
       $motif = $fields[$header{"motif"}];           
       next if (! -e "$motif");
       my %features =();
       foreach $f (keys %header){
	   $features{$f} = $fields[$header{$f}];	   
	   if ($header{$f} > $#fields){
	       $features{$f} = 0;
	   }
       }
       next if ($features{"Seq_id"} > 0.95 && $features{"Weight"}<= 4);       
       $motif_features{$motif}=\%features;
   }
}

close IN;

@all_motifs = reverse sort {$motif_features{$a}->{"Rank.index"} <=> $motif_features{$b}->{"Rank.index"}} keys %motif_features;
%group_motifs=();
foreach $motif (@all_motifs){    
    $motif =~ /^(\S+)\.motif/;
    $pref = $1;
    if (!exists $group_motifs{$pref}){
	$group_motifs{$pref} = [$motif];
    }       
    else{
	$temp = $group_motifs{$pref};
	push @$temp,$motif;
    }
}

%processed=();
foreach $motif ( @all_motifs){
    $motif =~ /^(\S+)\.motif/;
    $pref = $1;
    next if (exists $processed{$pref});
    $processed{$pref}= 1;
    $temp = $group_motifs{$pref};        
    
    @m = reverse sort {$motif_features{$a}->{"Rank.index"} <=> $motif_features{$b}->{"Rank.index"}} @$temp;
    print join " ", @m;
    print "\n";
    @selected = ();
    foreach $n (@m){
	last if ($motif_features{$n}->{"Rank.index"} < $threshold || scalar @selected >= 5);
	next if ($motif_features{$n}->{"Weight"} < 2.5);
	$redundant=0;
	foreach $s (@selected){
	    $structure_diff = `diff_motif_structure.pl $n $s`;
	    %diff = ($structure_diff=~/(\S+)=(\S+)/g);
	    if( $diff{"nolap_seq1"} < 0.2 * $diff{"olap_seq"} &&
		$diff{"olap_len"} * 0.5 >  $diff{"nolap_len1"}){
		if ($diff{"FP"}+ $diff{"FM"} < 0.5 * $diff{"TP"}){
		    $redundant=1;
		    last;
		}
		if ($motif_features{$n}->{"Rank.index"}  < 0.5 * $motif_features{$s}->{"Rank.index"}){
		    $redundant = 1;
		    last;
		}
	    }
	}
	next if ($redundant);
	push @selected, $n;
    }
    if ($pref =~ /^\S*\//){
	$pref =~ s/^(\S*\/)/$out_motif_dir\//;	    
    }
    else{
	$pref = "$out_motif_dir/$pref";
    }	
    for($i=0; $i <= $#selected; $i++){
	$motif_features{$selected[$i]}->{"motif"}= "$pref.motif.$i";		
	print $selected[$i], " $pref.motif.$i\n";
	copy("$selected[$i]", "$pref.motif.$i");		
    }
}
    

