#!/usr/bin/perl
#
use strict;
use warnings;
use Cwd 'abs_path';
use Parallel::ForkManager;
use FindBin qw($Bin);


my $bam = shift;
#my $workdir = shift;
#`mkdir -p $workdir`;

my %hs_basecnt = ();

my @init_tags = ("AS","XS","XM","XO","XG","NM","YS");
my %hs_multipleid = ();
#my $pm = new Parallel::ForkManager(10);
open(F,"$bam");
#open(F,"samtools view $bam|");
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
    my $alnid = $id.".".$hs_multipleid{$id};
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
    ### READ features 
	#my $rsize = `python3 $Bin/read_feature/get_complexity_score.stdout.py $flag $read $workdir test` ;
	#chomp($rsize);
	#print "\t$rsize";#11
    my $gc = GC($read);
    print "\t$gc"; #12
    my ($lowq,$avgq) = QUAL($qual);
    print "\t$lowq";#13
    print "\t$avgq"; #14
    ### REF features 
	#my $ref_seq = `perl $Bin/REF_feature/get_fa.pl $chr $start $end /mss6/nayoung/CNV_kibon/See_Mapping_pattern/PIPE/FeatureExtract/REF_feature/genome/  test $workdir ` ;
    #my ($ref_cpg, $ref_size) = split(/\t/,$test);
	#my $ref_cpg = GC($ref_seq);
	#my $ref_size = `python3 $Bin/read_feature/get_complexity_score.stdout.py $flag $ref_seq $workdir test`;
	#chomp($ref_size);
	#print "\t$ref_cpg";#17
	#print "\t$ref_size";#18
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
