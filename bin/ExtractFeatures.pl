#!/usr/bin/perl
#
use strict;
use warnings;
use Cwd 'abs_path';
use FindBin qw($Bin);
use File::Basename;



my $input = shift or die "Usage: $0 input.bam|input.sam|input.bed\n";
my $bed_file;

if ($input =~ /\.bam$/) {
    my $sam_file = basename($input, ".bam") . ".sam";
    system("samtools view -h $input > $sam_file") == 0 or die "Failed to convert BAM to SAM: $!";
    print "Converted BAM to SAM: $sam_file\n";

    $bed_file = basename($sam_file, ".sam") . ".bed";
    system("sam2bed < $sam_file > $bed_file") == 0 or die "Failed to convert SAM to BED: $!";
    print "Converted SAM to BED: $bed_file\n";
} elsif ($input =~ /\.sam$/) {
    $bed_file = basename($input, ".sam") . ".bed";
    system("sam2bed < $input > $bed_file") == 0 or die "Failed to convert SAM to BED: $!";
    print "Converted SAM to BED: $bed_file\n";
} elsif ($input =~ /\.bed$/) {
    $bed_file = $input;  
    print "Using existing BED file: $bed_file\n";
} else {
    die "Invalid input format. Please provide a BAM, SAM, or BED file.\n";
}
my %hs_basecnt = ();

my @init_tags = ("AS","XS","XM","XO","XG","NM","YS");
my %hs_multipleid = ();

open(F, $bed_file) or die "Failed to open BED file: $bed_file\n";
while(<F>){
    chomp;
	#$pm -> start and next;
    if ($_ =~ /^@/){next;}
    my @ar_tmp = split(/\t/,$_);
    my $chr = shift(@ar_tmp); 
    my $start = shift(@ar_tmp); 
    my $end = shift(@ar_tmp);
    my $id = shift(@ar_tmp);
    my $mapq = shift(@ar_tmp);
    my $order = shift(@ar_tmp);
    my $flag = shift(@ar_tmp);
    my $cigar = shift(@ar_tmp);
    my $mrnm = shift(@ar_tmp);
    my $mpos = shift(@ar_tmp);
    my $isize = shift(@ar_tmp);
    my $read = shift(@ar_tmp);
    my $qual = shift(@ar_tmp);

    $hs_multipleid{$id}++;
    my $alnid = $id
    my %hs_tag = ();
    foreach my $thistag (@init_tags){
        $hs_tag{$thistag} = "NA";
    }
    foreach my $thistag (@ar_tmp){
        my ($tag, $fmt, $val) = split(/:/,$thistag);
        $hs_tag{$tag} = $val;
    }
    print "$alnid\t$mapq"; #1 
    foreach my $thistag (@init_tags){ #2~8
        print "\t$hs_tag{$thistag}";    
    }
    my $diff = "NA";
    if ($hs_tag{"XS"} ne "NA" && $hs_tag{"AS"} ne "NA"){
        $diff = abs($hs_tag{"XS"}-$hs_tag{"AS"});
    }
    print "\t$diff";  #9 
    print "\t$isize";  #10
    my $gc = GC($read);
    print "\t$gc"; #11
    my ($lowq,$avgq) = QUAL($qual);
    print "\t$lowq";#12
    print "\t$avgq"; #13
    print "\n";
}
close(F);


sub GC {
    my $str = shift;
    %hs_basecnt = ("A" => 0 , "T" => 0 , "G" => 0 , "C" => 0 );
    foreach my $this (split(//,$str)){
        $hs_basecnt{$this} ++;
    }
	my $gc = 0;
	if(!($hs_basecnt{"A"}+$hs_basecnt{"C"}+$hs_basecnt{"G"}+$hs_basecnt{"T"}) == 0){
    $gc = ($hs_basecnt{"G"}+$hs_basecnt{"C"})/($hs_basecnt{"A"}+$hs_basecnt{"C"}+$hs_basecnt{"G"}+$hs_basecnt{"T"});
	
	}
    $gc = sprintf("%.3f", $gc);
    return $gc;
}

sub QUAL {
    my $seq = shift;
    my $total = 0;
    my $len = 0;
    my $cnt_under20 = 0;
    foreach my $b (split(//,$seq)){
        $total += ord($b);
        $len ++;
        if (ord($b)< 20){
            $cnt_under20 ++;
        }
    }
    my $avrg = $total / $len;
    $avrg = sprintf("%.3f", $avrg);

    return ($cnt_under20, $avrg)
}
