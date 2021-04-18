use Class::Struct;

## General String operation

sub remove_gap{
    my $seq=shift;     
    $seq =~ s/\.|-//g;
    return $seq;
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

sub pad_string{
    my ($seq, $n,$ch, $dir) = @_;
    my $l = $n  - length($seq) ;
    if ($l > 0) {
	# pad right;
	if ($dir == 1) {
	    $seq = $seq.make_string($l, $ch);
	}
	# pad left;
	if ($dir == 0) {
	    $seq = make_string($l, $ch).$seq;
	}	
    }
    return $seq;
}

#################
# FASTA FORMAT  #
#################
struct 'Seq' => {     		      
     acc =>'$', 
     id  =>'$',
     desc => '$',		      
     seq =>'$'
}; 		

     
sub read_fasta{
    my $file_name = shift @_;   
    open FASTA, $file_name or die "Can't open file $file_name!\n";    
    my %seqs = ();
    my $seq="";
    my $old = $/;
    my $acc="";
    my $desc ="";
    $/="\n>";
    my $i=1;
    while (<FASTA>){
	chomp;
	my @lines = split "\n", $_;
	if ($lines[0] =~ /^>{0,1}(\S+)\s*$/){
	    $acc=$1;
	    $desc="";
	}
	elsif ($lines[0] =~ /^>{0,1}(\S+)\s+(\S+.*)/){
	    $acc=  $1;
	    $desc= $2;
	}
	$seq = join "", @lines[1..$#lines];
	$new_seq = Seq->new(acc=>$acc, id=>$i, desc=>$desc, seq=>$seq);
	$seqs{$acc} = $new_seq;
	$i++;
    }
    close FASTA; 
    $/=$old;
    return \%seqs;
}

sub write_fasta{
    my $out_file = "";
    my $temp = shift @_;
    if (scalar @_ > 0){
	$out_file = shift @_;
    }
    my %seqs = %$temp;
    my $elem;

    if ($out_file ne "") {
	open OUT, ">$out_file" or die "Can't open file $out_file for writing";
	foreach $elem (sort {$a->id <=> $b->id} values %seqs){
	    next if (length($elem->seq) == 0);	    
	    print OUT ">", $elem->acc;
	    if ($elem->desc ne ""){
		print OUT "\t", $elem->desc;
	    }
	    print OUT "\n";
	    print OUT $elem->seq, "\n";
	}
	close OUT;
    }
    else{
	foreach $elem (sort {$a->id <=> $b->id} values %seqs){
	    next if (length($elem->seq) == 0);
	    print ">", $elem->acc,"\t", $elem->desc, "\n";
	    print $elem->seq, "\n";
	}
    }
}

#####################
# STOCKHOLM FORMAT  #
#####################

struct 'AlignSeq' => {     		      
     acc =>'$', 
     id  =>'$',
     desc => '$',		      
     weight=>'$',
     seq =>'$', 		      
     align_seq=>'$', 
     ss=>'$', 
     align_ss=>'$', 
     align_map=>'%', 
     rev_map=>'%', 

     start=>'$', 
     end=>'$', 
     score => '$'
}; 

struct 'Alignment' => {
     seqs => '%',
     flags => '%',
     ss_cons => '$',
     rf => '$',		       
     score => '$',	 
     weight => '$',
     len =>'$'		       
};		      



sub align_seq_map{
    my $motif=shift @_;
    my $seq = $motif->align_seq;
    my $j= $motif->start;
    my $dir = $motif->start < $motif->end;
    my @letters = split "", $seq;
    my $i;
    for($i=0; $i < length($seq); $i++) {
	if (! ($letters[$i] =~ /\.|-/)) {	
	    $motif->align_map($i, $j);
	    $motif->rev_map($j, $i);
	    if ($dir){
		$j++;
	    }
	    else{
		$j--;
	    }
	}
    }    
}

sub pair_table{
    my ($ss_str) = @_;
    my $i=$j=0;

    if (!defined($ss_str)) {
	die "!defined(ss_str), which is odd";
    }

    my $len = length($ss_str);
    my @ss = split "", $ss_str;
    my %pt=();
    my @stack=();
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

sub read_stockholm{
    my ($infile) = shift @_;
    my ($id, $acc, $weight,$score, $desc);
    my ($start, $end);
    my $ss_cons="";
    my $rf = "";
    my %seqs =();
    my %flags =();
    my $line;
    my $new_seq;
    open RFAM, $infile or die "Can't open file $infile!\n";
    my $new_id = 0;
    my $last_count="";
    my $last_acc="";
    while (<RFAM>) {
	$line= $_;
	if ($line =~ /^\#\s+/ || $line =~ /^\s*$/ || $line =~ /^\#=GF/){
	    $last_acc="";
	    $last_count= 0;
	    next;
	}      
	if ($line =~ /\#=GS\s+(\S+)\s+WT\s+(\S+)/ ){                
	    $acc = $1;
	    $weight = $2;
	    if ($last_acc ne $acc){
		$last_acc = $acc;
		$last_count = 0;
		$id = $acc;
	    }
	    else{
		$last_count++;
		$id = "$acc.$last_count";
	    }
	    if (! exists $seqs{$id}) {
		$new_seq = AlignSeq->new(acc=> $acc, weight =>$weight, id => $new_id);
		$flags{"WGT"} = 1;
		$new_id ++;
		$seqs{$id} = $new_seq;
	    }
	    else{
		$seqs{$id}->weight($weight);
	    }
	}
	elsif ($line =~ /\#=GS\s+(\S+)\s+DE\s+(.+)$/ ){
	    $acc = $1;
	    $desc = $2;
	    $flags{"DE"} = 1;
	    if ($last_acc ne $acc){
		$last_acc = $acc;
		$last_count = 0;
		$id = $acc;
	    }
	    else{
		$last_count++;
		$id = "$acc.$last_count";
	    }
	    if (!(exists $seqs{$id})) {
		$new_seq = AlignSeq->new(acc=> $acc, des=>$desc, id => $new_id, weight=>1);
		$new_id ++;
		$seqs{$id} = $new_seq;
	    }
	    else{
		$seqs{$id}->desc($desc);             
	    }
	    if ($desc =~ /(\d+)\.\.\s*(\d+)\s*(\S*)$/) {
		$start = $1;
		$end = $2;
		$flags{"START"} = 1;
		$flags{"END"} = 1;
		$seqs{$id}->start($start);             
		$seqs{$id}->end($end);
		if ($3 ne ""){
		    $score = $3;
		    $flags{"SC"} = 1;
		    $seqs{$id}->score($score);                          
		}
	    }
        }
	elsif ($line =~ /\#=GC\s+SS_cons\s+(\S+)/) {
	    $flags{"SS_cons"} = 1;
	    $ss_cons = $ss_cons.$1;
	}	
	elsif ($line =~ /\#=GC\s+RF\s+(\S+)/) {
	    $flags{"RF"} = 1;
	    $rf = $rf.$1;
	}	
	elsif ($line =~ /^\#=GR\s+(\S+)\s+SS\s+(\S+)/){
	    $acc =$1;
	    $ss = $2;
	    $flags{"SS"} = 1;
            if (!(exists $seqs{$id})) {
		$new_seq = AlignSeq->new(acc=> $acc, align_ss=>$ss, id => $new_id, weight=>1);
		$new_id ++;
		$seqs{$id} = $new_seq;
	    }
	    else{
		if (defined $seqs{$id}->align_ss ) {
		    $old_ss = $seqs{$id}->align_ss;
		    $seqs{$id}->align_ss($old_ss.$ss);
		}
		else {
		    $seqs{$id}->align_ss( $ss);
		}
            }
	}	
	elsif ($line =~ /^\#/ ){
	    print STDERR "Unrecognized : $line \n";
	}
	elsif ($line =~ /^(\S+)\s+(\S+)/){
	    $acc = $1;
	    $seq = $2;
	    if ($last_acc ne $acc){
		$last_acc = $acc;
		$last_count = 0;
		$id = $acc;
	    }
	    else{
		$last_count++;
		$id = "$acc.$last_count";
	    }
	    if (!(exists $seqs{$id})) {
		$new_seq = AlignSeq->new(acc=> $acc, align_seq=>$seq, id => $new_id, weight=>1);
		$new_id ++;
		$seqs{$id} = $new_seq;
	    }
	    else{
		if (defined $seqs{$id}->align_seq ) {
		    $old_seq = $seqs{$id}->align_seq;
		    $seqs{$id}->align_seq($old_seq.$seq);
		}
		else {
		    $seqs{$id}->align_seq($seq);
		}
	    }
	}      
	else{
	    $last_acc="";
	    $last_count=0;
	}
    }
    close RFAM;
    $sum_score=0;
    $sum_weight=0;
    $sum_len = 0;
    foreach $id (keys %seqs) {
	$seqs{$id}->seq(remove_gap( $seqs{$id}->align_seq));	
	if (exists $flags{"SS"} && $flags{"SS"} == 1){
	    $seqs{$id}->ss(remove_gap( $seqs{$id}->align_ss));
	}
	if (exists $flags{"SC"}){
	    $sum_score += $seqs{$id}->score * $seqs{$id}->weight;
	}
	$sum_weight += $seqs{$id}->weight;
	$sum_len += length($seqs{$id}->seq) * $seqs{$id}->weight;	
    }    
    if ($sum_weight == 0){
	$sum_weight = scalar keys %seqs;
    }
    return (Alignment->new(seqs=> \%seqs, flags => \%flags, weight=>$sum_weight, 
			   ss_cons => $ss_cons,rf => $rf, 
			   score => $sum_score/$sum_weight,len=>$sum_len/$sum_weight)); 
}


sub read_rfam{
    my $rfam_file = shift @_;
    my $alignment = read_stockholm($rfam_file);
    my $temp = $alignment->seqs;
    my %seqs = %$temp;
    my $acc;
    foreach $acc (keys %seqs) {
	if ($acc =~ /(\S+)\/(\d+)-(\d+)/){
	    $seqs{$acc}->start($2);
	    $seqs{$acc}->end($3);
	    $seqs{$1}=$seqs{$acc};
	    delete($seqs{$acc});
	}
    } 
    $alignment->seqs(\%seqs);
    return $alignment;
}


sub write_stockholm
{
    my $alignment = shift @_;
    my $file = shift @_;    
    open OUT, ">$file" or die "Can't open $file for writing";
    print OUT "# STOCKHOLM 1.0\n\n";
    my $max_name_length = 0;
    my $line_len = 80;
    my $temp = ($alignment->seqs); 
    my %seqs = %$temp;
    if (scalar keys %seqs == 0) {
	die "Empty alignments";
    }
    $temp =  $alignment->flags; 
    my %flags = %$temp;
    my $ss_cons = $alignment->ss_cons;
    my $rf = $alignment->rf;
    my $id;
    foreach $id (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {
	$temp = length($seqs{$id}->acc);
	$max_name_length = $temp if $temp > $max_name_length;
    }
    $max_name_length++;  
    my $name_format = "%-".$max_name_length."s";
    my $acc;
    my $seq;
    if (exists $flags{"WGT"}) {
	foreach $id (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {
	    print OUT "#=GS ";
	    printf OUT $name_format, $seqs{$id}->acc;
	    print OUT "WT\t", $seqs{$id}->weight, "\n";
	}
	print OUT "\n\n";
    }
    if (exists $flags{"DE"}) {
	foreach $id (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {
	    print OUT "#=GS ";
	    printf OUT $name_format, $seqs{$id}->acc;
	    print OUT "DE\t", $seqs{$id}->desc, "\n";
	}
	print OUT "\n\n";
    }        
    
    $temp_id = (keys %seqs)[0];
    my $ss_len = length($seqs{$temp_id}->align_ss);
    my $seq_len = $ss_len;
    my $gap1 = "            ";
    my $gap2 = "     ";
    for(my $len=0; $len < $seq_len; $len += $line_len){
	my $l = $line_len;
	$l = $ss_len - $len if ($ss_len - $len < $l)  ;
	foreach $id (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {
	    $acc = $seqs{$id}->acc;
	    $seq = substr $seqs{$id}->align_seq, $len, $l;
	    printf OUT $name_format, $acc;
	    print OUT $gap1, $seq, "\n";	    
	    if (exists $flags{"SS"}){
		my $ss = substr $seqs{$id}->align_ss, $len, $l;
		print OUT "#=GR ";
		printf OUT $name_format, $acc;
		print OUT "SS", $gap2, $ss, "\n";
	    }
	}
	if (exists $flags{"SS_cons"}){
	    my $ss_l = substr $ss_cons, $len, $l;
	    printf OUT $name_format, "#=GC SS_cons";
	    print OUT $gap1, $ss_l, "\n";
	}
	if (exists $flags{"RF"}){
	    my $rf_l = substr $rf, $len, $l;	
	    printf OUT $name_format, "#=GC RF";
	    print OUT $gap1, $rf_l, "\n";
	}
	print OUT "\n";
    }
    print OUT "//\n";
    close(OUT);
}

sub write_selex
{
    my $alignment = shift @_;
    my $file = shift @_;    
    open OUT, ">$file" or die "Can not open $file for writing";

    my $max_name_length = 0;
    my $line_len = 80;
    my $temp = ($alignment->seqs); 
    my %seqs = %$temp;
    $temp =  $alignment->flags; 
    my %flags = %$temp;
    my $ss_cons = $alignment->ss_cons;
    
    foreach $acc (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {
	$temp = length($seqs{$acc}->id."_".$acc) ;
	$max_name_length = $temp if $temp > $max_name_length;
    }
    $max_name_length++;  
    $name_format = "%-".$max_name_length."s";
    print OUT "#=AU CM RNA automatic alignment\n";
    foreach $acc (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {    
	print OUT "#=SQ ";
	printf OUT $name_format, $seqs{$acc}->id."_".$acc, "\t";
	print OUT "1.0000 - - ", $seqs{$acc}->start, "..", $seqs{$acc}->end, "::0 ", $seqs{$acc}->score, "\n" ;
	$seq_len = length($seqs{$acc}->align_seq);
    }    

    for($len=0; $len < $seq_len; $len += $line_len){
	$ss_len = length($ss_cons);
	$l = $line_len;
	$l = $ss_len - $len if ($ss_len - $len < $l)  ;
	if (length($ss_cons) > 0){
	    $ss_c = substr $ss_cons, $len, $l;
	    printf OUT $name_format, "#=CS"; 
	    print  OUT $ss_c, "\n";
	}
	foreach $acc (sort {$seqs{$a}->id <=> $seqs{$b}->id} keys %seqs) {
	    $seq_len = length($seqs{$acc}->align_seq);
	    my $name = $seqs{$acc}->acc;
	    $l = $line_len;
	    $l = $seq_len - $len if ($seq_len - $len < $l)  ;	    
	    my $seq = substr $seqs{$acc}->align_seq, $len, $l;
	    my $ss  = substr $seqs{$acc}->align_ss,  $len, $l;
	    printf OUT $name_format, $seqs{$acc}->id."_".$name;
	    print  OUT $seq, "\n";
	    printf OUT $name_format, "#=SS";
	    print  OUT $ss, "\n";
	}
	print OUT "\n";
    }
}

#$file=shift @ARGV;
#$temp = read_stockholm($file);
#write_stockholm($temp, "temp");
