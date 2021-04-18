#!/usr/local/bin/perl -w

use Getopt::Long;
use Class::Struct;

$fam1_file=shift;
$fam2_file=shift;
$format = shift;
$FileName=shift;
$line_len = 210;

$maxHelixBulge=3; # >2 bases means we're not in a helix any more
$organismAssocFileName="";
$familyName = "";
$explanation = "";
%normalPairs=("AU",1,"UA",1,"CG",1,"GC",1,"GU",1,"UG",1);

struct 'MyMotif' => {seq_id=>'$', start=>'$', end=>'$', rev_flag=>'$', seq=>'$', align_seq=>'$', ss=>'$', align_ss=>'$', align_map=>'%', rev_map=>'%', helix_num=>'%'}; 

# first color is the one for bad base pairs
@colors=("#d0d0d0","#ff9999","#9999ff","#99ff99","#FF9900","#99ffff","#98C0C0","#ffff99","#ff33ff","#33ffff","#ffff33","#2ED70C","#F4BEF8","#ff9900","#B94F32","#FF0000","#ffcc99","#CCCCCC","#CC3366","#CCff66","#Ffcc66","#7e652f");


sub max{
    @a = @_;
    $m = $a[0];
    for($i =1; $i < length(@a); $i++) {
	if ($a[$i] > $m){
	    $m = $a[$i];
	}
    }
    return $m;
}

sub max2{
    ($a, $b)=@_;
    return $a if ($a > $b);
    return $b;
}

sub match_ss{
    my ($seq, $ss_cons) = @_;
    my $ss= $ss_cons;
    return $ss;
}

sub remove_gap{
    my $seq=shift;
    return join "", (split /\./, $seq);
}

sub remove_gap_ss{
    my ($align_seq, $align_ss) = @_;
    my $ss="";
    my $i;
    for($i=0; $i < length($align_seq); $i++) {
	$ch = substr $align_seq, $i, 1;
	if ($ch ne '.') {
	    $ss=$ss.(substr $align_ss, $i, 1);
	}
    }
    return $ss;
}

sub make_string{
    my ($n,$ch) = @_;
    my $b = "";
    my $i;
    for($i =0; $i < $n; $i++){
	$b = $b.$ch;
    }
    return $b;
}

sub align_seq_map{
    my $motif=shift @_;
    my $seq = $motif->align_seq;
    my $j= $motif->start;
    my @letters = split "", $seq;
    for($i=0; $i < length($seq); $i++) {
	if (! ($letters[$i] =~ /[\.-]/)) {	
	    $motif->align_map($i, $j);
	    $motif->rev_map($j, $i);
	    if ($motif->rev_flag) {
		$j--;
	    }
	    else{
		$j++;
	    }
	}
    }
}

sub read_stockholm{
    my ($rfam_file) = shift;
    my ($start, $end, $rev_flag);
    my %rfam_motif=();
    my $ss_cons="";
    my %scores = ();
    open RFAM, $rfam_file or die "Can't open file $rfam_file!\n";

    while (<RFAM>) {
	$line= $_;
	next if ($line =~ /^\s*$/);	    
	if ($line =~ /\#=GC SS_cons\s+(\S+)/) {
	    $ss_cons = $ss_cons.$1;
	}	
	elsif ($line =~ /^(\S+)\/(\d+)-(\d+)\s+(\S+)/){
	    $seq_id = $1;
	    $start = $2;
	    $end = $3;
	    $seq = $4;
	    $rev_flag = $start > $end;
	    $index = "$seq_id";
	    if (exists $rfam_motif{$index} ) {
		$old_seq = $rfam_motif{$index}->align_seq;
		$rfam_motif{$index}->align_seq( $old_seq.$seq);
	    }
	    else {
		$rfam_motif{$index} = MyMotif->new(seq_id => $seq_id, start => $start, end => $end, rev_flag=> $rev_flag, align_seq=>$seq);		
	    }
	}
	elsif ($line =~ /^\#=GR (\S+)\/.*\s+SS\s+(\S+)/){
	    $seq_id =$1;
	    $ss = $2;
	    $index = "$seq_id";
	    if (defined $rfam_motif{$index}->align_ss ) {
		$old_ss = $rfam_motif{$index}->align_ss;
		$rfam_motif{$index}->align_ss( $old_ss.$ss);
	    }
	    else {
		$rfam_motif{$index}->align_ss( $ss);
	    }
	}
	elsif ($line =~ /^\#=GS\s+(\S+)\/.*\s+DE\s+(\S+)/) {
	    $seq_id = $1;
	    $score = $2;
	    $scores{$seq_id} = $score;
	}
	else{
	    #print "Can't parse line : $line\n";
	}
    
    }
    close RFAM;
    foreach $id (keys %rfam_motif) {
	$rfam_motif{$id}->seq(remove_gap( $rfam_motif{$id}->align_seq));
	if (! defined $rfam_motif{$id}->align_ss || length($rfam_motif{$id}->align_ss) <= 0 ) {
	    $rfam_motif{$id}->align_ss( match_ss($rfam_motif{$id}->align_seq,  $ss_cons));
	}
	$rfam_motif{$id}->ss(remove_gap_ss( $rfam_motif{$id}->align_seq,$rfam_motif{$id}->align_ss ));
	align_seq_map($rfam_motif{$id});
    }
    return ($ss_cons, \%rfam_motif, \%scores);
}


sub pair_table{
    my $ss_str = shift;
    my $i=$j=0;
    #find the pair_table of the align_seq
    #print "ss  ", @ss, "\n";
    my $len = length($ss_str);
    my @ss = split "", $ss_str;   
    my %pt=();
    while($i < $len) {
	while($ss[$i] eq '<'){
	    $stack[$j] = $i;
	    $j++;
	    $i++;
	}
	while($i < $len && $ss[$i] eq '>') {
	    $j--;
	    $pt{$i} = $stack[$j];
	    $pt{$stack[$j]} = $i;
	    $i++;
	}
	if ($i < $len && $ss[$i] ne '<' && $ss[$i] ne '>'){
	    $i++;
	}
    }
    return %pt;
}

sub pair_table_motif{
    my ($motif) = shift @_;
    my ($i, $j, $i1, $j1);
    my $len = length($motif->align_ss);
    my %pt=();
    my $m = $motif->align_map;

    my %pt_temp = pair_table($motif->align_ss);


    for($i=0; $i < $len; $i++) {
	next if !exists $m->{$i};
	$i1 =  $m->{$i};
	next if !exists $pt_temp{$i};
	$j = $pt_temp{$i};
	next if !exists $m->{$j};
	$j1 = $m->{$j};
	$pt{$i1} = $j1;
	$pt{$j1} = $i1;
    }
    return %pt;
}    



sub map_helix{
    my ($ss_str, $pt) =  @_;
    my %posToHelixNum =();
    my $i;
    my @ss = split "", $ss_str;   
    my %helix_map=();
    

    for($i=0; $i < scalar @ss; $i++) {
	next if (!exists $pt->{$i} || $pt->{$i} < $i); 
	$connect = 0;
	for ($j=$i-1,$gap1 = 0; $j>=0 && $gap1 < $maxHelixBulge; $j--) {	   
	    if ( exists $pt->{$j} && $pt->{$j} > $j) {
		if ($pt->{$j}  - $pt->{$i} -1  <= $maxHelixBulge ) {
		    $connect=1;
		    last;
		}	    
		$gap2 = 0; 
		for($k =  $pt->{$i} + 1; $k <= $pt->{$j} - 1; $k++){
		    $gap2 ++ if $ss[$k] ne '.';
		}
		#print (substr $ss_str, $pt->{$i} + 1, $pt->{$j} - $pt->{$i}), " ", $gap2, "\n";
		#print $j, " " , $i, " " , $pt->{$i}, " ", $pt->{$j}, "\n";
		if ($gap2 <= $maxHelixBulge){
		    $connect = 1 ;
		    last;
		}			    
	    }
	    next if ($ss[$j] eq '.');
	    $gap1++;
	}

	if ($connect) {
	    # same helix
	}
	else {
	    # new helix
	    $currHelixNum++;
	}
	#print " $i $currHelixNum\n";
	$posToHelixNum{$i}=$currHelixNum;
	$posToHelixNum{$pt->{$i}}=$currHelixNum;
    }
    return %posToHelixNum;
}


sub map_helix_motif{
    my ($motif, $cons_pt, $pos) = @_;    
    my $ss_str = $motif->align_ss;
    my %pt = pair_table($ss_str);    

    #foreach $i (sort keys %$pos){
    #	print $i, " ", $pos->{$i}, "\t";
    #}
    #print "\n";
    
    #print $motif->seq_id, "\n";
    #foreach $i (sort keys %pt){
    #	print $i, "_", $pt{$i}, " ";
    #}
    #print "\n";
    my @ss = split "", $ss_str;  
    my @seq = split "", $motif->align_seq;
    #print @ss, "\n";
    #print @seq, "\n";
    for($i=0; $i < scalar @ss; $i++) {
	if (exists $pt{$i} ) {
	    next if $pt{$i} < $i;	    
	    $pair = $seq[$i].$seq[$pt{$i}];
	    #print $pair, "\n";
	    if ($normalPairs{$pair}) {
		if (exists $pos->{$i} ) {
		    $motif->helix_num($i, $pos->{$i});
		    $motif->helix_num($pt{$i}, $pos->{$i});
		}
		else{
		    # find the closest helix
		    for($j1=$i+1; $j1 < $pt{$i}; $j1++) {
			last if exists $pos->{$j1};
		    }
		    if (exists $pos->{$j1}){
			$diff = max2($j - $i, $pt{$i}-  $cons_pt->{$j1} );
			$j = $j1;
		    }
		    for($j2=$i-1; $j2 > 0; $j2--) {
			last if exists $pos->{$j2};
		    }
		    if (exists $pos->{$j2}) {
			$diff2 =  max2($i - $j , $cons_pt->{$j} - $pt{$i});
			if (defined $diff && $diff2 < $diff || !defined $diff) {
			    $diff = $diff2;
			    $j = $j2;
			}
		    }
		    if (defined $diff) {
			$motif->helix_num($i, $pos->{$j});
			$motif->helix_num($pt{$i} , $pos->{$j});
		    }
		}
	    }
	    else{
		$motif->helix_num($i, -1);
		$motif->helix_num($pt{$i}, -1);
	    }
	}
	else{
	    if (exists $pos->{$i}) {
		$motif->helix_num($pt{$i}, -1);
	    }
	}
    }
        
}

#overlay the second motif on the first motif. 
sub map_motif{
    my ($motif1, $motif2) = @_;
    
    my ($pre_seq, $pre_ss, $suf_seq, $suf_ss, $map_ss);
    my ($i, $j);
    my ($start1, $start2, $end1, $end2) = ($motif1->start, $motif2->start, $motif1->end, $motif2->end);

    $rev_flag = 0;
    $rev_flag = 1 if ($motif1->start > $motif1->end);	
        
    ##############
    if ( !$rev_flag && ($end1 < $start2 || $start1 > $end2) ||
	 $rev_flag && ($end1 > $start2 || $start1 < $end2)){      #No overlap
	return ();
    }
	

    if ($start1 <= $start2 && !$rev_flag || $start1 >= $start2 && $rev_flag) {
	    $pre_seq="";
	    $pre_ss = "";
	}
    else{
	$len = abs($start1 - $start2);
	$pre_seq = substr $motif2->seq, 0, $len;
	$pre_ss = substr $motif2->ss, 0, $len;	    
    }


    if ($end1 >= $end2  && !$rev_flag || $end2 >= $end1 && $rev_flag ) {
	$suf_seq="";
	$suf_ss = "";
    }
    else{
	$len = abs($end2 - $end1);
	$suf_seq = substr $motif2->seq, -$len;
	$suf_ss = substr $motif2->ss, -$len;	    
    }		
    

    my $map = $motif1->align_map;
    my %pt = pair_table_motif($motif2);
    my @ss=();
    
    
    for($i=0; $i < length($motif1->align_ss); $i++) {
	if (exists $map->{$i}) {
	    $j = $map->{$i};
	    if ($j < $motif2->start && $j <$motif2->end ||
		$j > $motif2->start && $j >$motif2->end) {
		$ss[$i] = '_';
		next;
	    }
	    if (exists $pt{$j}) {
		if ($pt{$j} > $j && $motif1->start < $motif1->end ||
		    $pt{$j} < $j && $motif1->start > $motif1->end ) {
		    $ss[$i] = '<';
		}
		else{
		    $ss[$i] = '>';
		}
		next;
	    }
	}
	if ($j < $motif2->start && $j < $motif2->end ||
	    $j > $motif2->start && $j > $motif2->end) {
	    $ss[$i] = '_';
	}
	else{
	    $ss[$i] = '.';
	}
    }

    $map_ss = join "", @ss;
    
    return ($pre_seq, $pre_ss, $suf_seq, $suf_ss, $map_ss);
}

sub parse_org{
    %emblIdToSpecies=();
    if (length($organismAssocFileName)>0) {
	open(ORGS,"$organismAssocFileName");
	while (<ORGS>) {
	    s/[\r\n]//g;
	    
	    ($emblId,$os,$oc)=split /\t/,$_;
	    
	    # just take first 2 words of species
	    $species=$os;
	    $_=$os;
	    if (/^([^ ]* [^ ]*)/) {
	    $species=$1;
	}
	    $emblIdToSpecies{$emblId}=$species;
	}
	close(ORGS)
    }    
}

sub plot_html{
    my $html_file = shift @_;
    open(HTML,">$html_file");
    
    print HTML "<html><head>";
    if ($familyName eq "") {
    }
    else {
	print HTML "<title>$familyName: multiple alignment of Rfam alignment and MFIND alignment</title>";
    }
    if (1) {
	#inline style
	print HTML "<style>b, .b { font-weight: normal;}\n";
	print HTML "#reallyBold { font-weight:bold; }\n";
	$i=0;
	for $color (@colors) {
	    print HTML "#b$i {background-color: $color;}\n";
	    $i++;
	}
	print HTML "#bad {background-color: $colors[0];}\n";
	print HTML "#shade {background-color: #DDDDDD;}\n";
	print HTML "</style>\n";
    }
    else {
	# .css file
	print HTML "<link rel=stylesheet href=\"./FancifyStockholm.css\" type=\"text/css\">";
    }
    print HTML "</head><body>\n";
        
    print HTML "<table>\n";
    print HTML "<tr>";
  
    #specfies header
    print HTML "<th><b ID=\"reallyBold\">Seq ID</b></th>";
    if (length($organismAssocFileName)>0) {
	print HTML "<th><b ID=\"reallyBold\">Species</b></th>";
    }
    #score header
    if (scalar(keys %$scores2)>0) {
	print HTML "<th> &#032 &#032 &#032 &#032 &#032 <b ID=\"reallyBold\">Score (bits)</b></th>";
    }
    else {
	print HTML "<th></th>";
    }

    print HTML "<th align=right></th> <th><b ID=\"reallyBold\">Sequence</b></th>  <th align=left></th>";
    print HTML "</tr>\n";
    foreach $id (sort keys %$fam1_motif) {
	if (exists $fam2_motif->{$id}) {
	    next if !exists $motif_overlap{$id};       #No overlap
	    $seq = $fam1_motif->{$id}->align_seq;
	    $helix_num = $fam1_motif->{$id}->helix_num;
	    $pre_seq = $motif_overlap{$id}{'pre_seq'};
	    $pre_ss = $motif_overlap{$id}{'pre_ss'};
	    $suf_seq = $motif_overlap{$id}{'suf_seq'};
	    $suf_ss = $motif_overlap{$id}{'suf_ss'};
	    $map_ss =  $motif_overlap{$id}{'map_ss'};
	
	    #print the id
	    print HTML  "<td><nobr>$id</nobr></td>";
	    #print the species
	    if (length($organismAssocFileName)>0) {
	    	print HTML "<td><nobr>$emblIdToSpecies{$id}</nobr></td>";
	    }
	    #print the score
	    if (exists $scores2->{$id}) {
		printf HTML "<td align=center> %3.2f </td>",$scores2->{$id};
	    }
	    else{
		print HTML "<td><nobr></nobr></td>";	
	    }

	    $onMouseOver="JavaScript:window.status='$id'";
	    #print sequence prefix
	    print HTML "<td align=right> <nobr><tt> <b>";
	    if (length($pre_seq) > 0) {
		print HTML "<b  id=\"shade\">$pre_seq </b>";
	    }

	    print HTML "</b> </tt></nobr> </td> ";
	    #print the colored seq
	    print HTML "<td onMouseOver=\"$onMouseOver\" onMouseOut=\"JavaScript:window.status=''\"><tt><nobr>";
	    for ($i=0; $i<length($seq); $i++) {
		$ch=substr $seq,$i,1;
		if (!exists $helix_num->{$i}) {
		    print HTML $ch;
		}
		else{
		    $h = $helix_num->{$i};
		    if ($h > 0) {
			print HTML "<b id=\"b$h\">$ch</b>";
		    }
		    else{
			print HTML "<b id=\"bad\">$ch</b>";
		    }
		}
	    }
	    print HTML "</nobr></tt></td>";
	    #print the sequence suffix
	    print HTML "<td><tt>";
	    if (length($suf_seq) > 0 ) {
		print HTML "<b  id=\"shade\">$suf_seq </b>";
	    }
	    print HTML "</tt> </td>";
	    print HTML "</tr>\n";		

	    

	    #Now the the secondary structure annotion of motif2
	    print HTML "</td><td>";

	    if (length($organismAssocFileName)>0) {
		print HTML "</td><td>";
	    }
	    
	    if (scalar(keys %$scores2)>0) {		
		print HTML "</td><td>";
	    }
	    print HTML "<td align=right> <tt>";
	    if (length($pre_ss) > 0 ) {
		print HTML "<b  id=\"shade\>"$pre_ss </b>";
	    }	    
	    print HTML "</tt></td>  <td><tt>", $map_ss,"</tt></td>";
	    print HTML "<td ><tt> ";
	    if (length($suf_ss) > 0 ) {
		print HTML "<b  id=\"shade\">$suf_ss </b>";
	    }	  
	    print HTML "</tt> </td>  </tr>\n";	    
	}
    }
    
    print HTML "</table></body></html>\n";
    close HTML;
}

sub plot_latex{
    my $latex_file = shift @_;
    open(LATEX,">$latex_file");
    print LATEX "\\documentclass{article}\n\\usepackage{color}\n\\setlength{\\fboxsep}{0pt}\n";
    print LATEX "\\usepackage{longtable}\n";
    print LATEX "\\usepackage{amsmath}\n";

    
    $i=0;
    for $color (@colors) {
	@component=();
	push @component,lc(substr $color,1,2);
	push @component,lc(substr $color,3,2);
	push @component,lc(substr $color,5,2);
	# print "$color,@component\n";

	print LATEX "\\definecolor{b$i}{rgb}{";
	for ($j=0; $j<3; $j++) {
	    if ($j>0) {
		print LATEX ",";
	    }
	    $level=hex($component[$j])/255.0;
	    print LATEX "$level";
	}
	print LATEX "}\n";
	$i++;
    }
    
    #print LATEX "%pdflatex params & extraction from hyperref package\n";
    #print LATEX "\\pdfhorigin0pt\n\\pdfvorigin0pt\n\\pdfpagewidth0pt\n";
    #print LATEX "\\pdfpageheight8.5in\n";
    #print LATEX "\\newcommand{\\myhref}[2]{\\pdfstartlink attr{/C [0 0 0.9] /Border [0 0 1]} user{/Subtype/Link/A<</Type/Action/S/URI/URI(#1)>>} #2 \\pdfendlink }\n";

    #print LATEX "\\addtolength\\textwidth{900pt}\n\\addtolength\\textwidth{72pt}\n";
    #print LATEX "\\setlength\\textheight{6.5in}\n";
    
    print LATEX "\\topmargin 0.0in\n";
    print LATEX "\\oddsidemargin 0.0in\n";
    print LATEX "\\headheight 0.0in\n";
    print LATEX "\\headsep 0.0in \n";
    print LATEX "\%\\textheight 7.0in \n";    
    print LATEX "\\pagestyle{empty}\n";

    print LATEX "\\begin{document}\n";
    $_=$explanation;
    s/\"([^\"]*)\"/``$1''/g;
    s/underlined/boxed/g;
    print LATEX "\\parbox{7in}{$_}\\\\\n";
    print LATEX "\\vspace*{3ex}\n";

    print LATEX "\\begin{longtable}{r";
    if (length($organismAssocFileName)>0) {
	print LATEX "c";
    }
    if (scalar(keys %$scores2)>0) {    
	print LATEX "c";
    }
    print LATEX "ccc}\n";

    # secondary structure string goes at the bottom of each page

    #$ss_cons1 =~ s/</&lt;/g;
    #$ss_cons1 =~ s/>/&gt;/g;
    $ss_cons1 =~ s/[-_]/\./g;

    # headers
    #ID
    print LATEX "{\\bf Seq ID}& ";
   
    #species
    if (length($organismAssocFileName)>0) {
	print LATEX "{\\bf Species}&";
    }

    #score
    if (scalar(keys %$scores2)>0) {    
	print LATEX "{\\bf Score (bits)}& ";
    }
    #sequence
    print LATEX "&{\\bf Sequence}\\\\\n";
  
    $align_len = length($ss_cons1);

    for($l = 0; $l < $align_len; $l += $line_len){
	
	foreach $id (sort keys %$fam1_motif) {    
	    if (exists $fam2_motif->{$id}) {
		next if !exists $motif_overlap{$id};       #No overlap
		$seq = $fam1_motif->{$id}->align_seq;
		$helix_num = $fam1_motif->{$id}->helix_num;
		$pre_seq = $motif_overlap{$id}{'pre_seq'};
		$pre_ss = $motif_overlap{$id}{'pre_ss'};
		$suf_seq = $motif_overlap{$id}{'suf_seq'};
		$suf_ss = $motif_overlap{$id}{'suf_ss'};
		$map_ss =  $motif_overlap{$id}{'map_ss'};
		$pre_ss =~ s/_/-/g;
		$map_ss =~ s/_/-/g;
		$suf_ss =~ s/_/-/g;
		

		print LATEX $fam1_motif->{$id}->seq_id, "&";
		if (scalar(keys %$scores2)>0) {    
		    print LATEX $scores2->{$id}, "&";
		}
		if (length($organismAssocFileName)>0) {
		    print LATEX "$emblIdToSpecies{$id}&";	
		}
		
		if (length($pre_seq) > 0) {
		    if ( $l == 0) {
			print LATEX "$pre_seq";
		    }
		}
		print LATEX "&";		    
		print LATEX "{\\tt ";
	    
	    
		$seq = $fam1_motif->{$id}->align_seq;
		$helix_num = $fam1_motif->{$id}->helix_num;
		
		for ($i=$l; $i<length($seq) && $i < $l + $line_len; $i++) {
		    $ch=substr $seq,$i,1;
		    if (!exists $helix_num->{$i}) {
			print LATEX $ch;
		    }
		    else{
			$h = $helix_num->{$i};
			if ($h > 0) {
			    print LATEX "\\colorbox{b$h}{$ch}";
			}
			else{
			    print LATEX "\\colorbox{b0}{$ch}";
			}
		    }
		}
		print LATEX "}";
		#print the sequence suffix
		print LATEX "&";
		if (length($suf_seq) > 0 ) {
		    if ($l == 0) {
			print LATEX "$suf_seq";
		    }
		}
		print LATEX "\\\\\n";			    

		#Now the the secondary structure annotion of motif2
		print LATEX "&";

		if (length($organismAssocFileName)>0) {
		    print LATEX "&";
		}
		
		if (scalar(keys %$scores2)>0) {		
		    print LATEX "&";
		}
		if (length($pre_ss) > 0 && $l == 0) {
		    print LATEX "{\\tt $pre_ss}";
		}
		
		print LATEX "&";
		

		$len = length($map_ss) - $l;
		$len = $line_len if ($len > $line_len) ;	
		$ss = substr $map_ss, $l, $len;
		print LATEX "{\\tt", $ss, "}";

		if (length($suf_ss) > 0 && $l == 0) {
		    print LATEX "& {\\tt $suf_ss}";
		}	  
		print LATEX "\\\\\n";			    
	    }
	    
	}
	print LATEX "Conserved SS&";
	
	if (length($organismAssocFileName)>0) {
	    print LATEX "&";
	}
	
	if (scalar(keys %$scores2)>0) {		
	    print LATEX "&";
	}
	$len = length($ss_cons1) - $l;
	$len = $line_len if ($len > $line_len) ;	
	$ss = substr $ss_cons1, $l, $len;
	print LATEX "&{\\tt $ss}\\\\\n";
	
    }
    print LATEX "\\end{longtable}\n";
    print LATEX "\\end{document}\n";
}


($ss_cons1, $fam1_motif, $scores1) = read_stockholm($fam1_file);
($ss_cons2, $fam2_motif, $scores2) = read_stockholm($fam2_file);


$ss_cons1 =~ s/\(/\</g;
$ss_cons1 =~ s/\)/\>/g;
$ss_cons1 =~ s/\]/\>/g;
$ss_cons1 =~ s/\[/\</g;

%pt = pair_table($ss_cons1);
foreach $id (sort {$a <=> $b} keys %pt) {
    #print $id, "\t", $pt{$id}, "\n";
}

%helix_map = map_helix($ss_cons1, \%pt);
%motif_overlap = ();
$max_pre_length=0;
foreach $id(sort keys %$fam1_motif) {
    if (exists $fam2_motif->{$id}) {	
	@result = map_motif($fam1_motif->{$id}, $fam2_motif->{$id});
	if ( scalar @result == 5) {
	    ($pre_seq, $pre_ss, $suf_seq, $suf_ss, $map_ss) = @result;
	    if (length($pre_seq) > $max_pre_length){
		$max_pre_length = length($pre_seq);
	    }	
	    $motif_overlap{$id} = {pre_seq=>$pre_seq, pre_ss => $pre_ss, suf_seq=>$suf_seq, suf_ss=>$suf_ss, map_ss=>$map_ss};     
	    
	    map_helix_motif($fam1_motif->{$id}, \%pt, \%helix_map);
	}

    }

}

if ($format eq "-h") {
    plot_html($FileName);
}
else {
    plot_latex($FileName);
}
#plot_html($htmlFileName);
#$format "%".$max_pre_length,"s";
#foreach $id(sort keys %$fam1_motif) {
    #sprintf $pre_seq $format   $motif_overlap{$id}->{'pre_seq'};
    #sprintf $pre_ss $format   $motif_overlap{$id}->{'pre_ss'};    
#}


