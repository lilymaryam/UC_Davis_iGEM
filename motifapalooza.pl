#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my $SEQS = 10;
my $LEN = 100;
my $FREQ = 0.9;
my $MPS  = 1;
my $PA   = 0.25;
my $PC   = 0.25;
my $PG   = 0.25;
my $PT   = 0.25;

my $usage = "
usage: $0 [options] <motif file>
options:
  -n <int>   number of sequences to generate [$SEQS]
  -l <int>   length of sequences to generate [$LEN]
  -f <float> frequency of sequences with motifs [$FREQ]
  -m <int>   number of motifs per sequence [$MPS]
  -b         include motifs on both strands [off]
  -a <float> background probablility of A [$PA]
  -c <float> background probablility of C [$PC]
  -g <float> background probablility of G [$PG]
  -t <float> background probablility of T [$PT]
";

our ($opt_n, $opt_l, $opt_f, $opt_m, $opt_b,
	$opt_a, $opt_c, $opt_g, $opt_t);
getopts('n:l:f:m:ba:c:g:t:');
$SEQS = $opt_n if $opt_n;
$LEN  = $opt_l if $opt_l;
$FREQ = $opt_f if $opt_f;
$MPS  = $opt_m if $opt_m;
$PA   = $opt_a if $opt_a;
$PC   = $opt_c if $opt_c;
$PG   = $opt_g if $opt_g;
$PT   = $opt_t if $opt_t;
my $BOTH_STRANDS = $opt_b;

die $usage unless @ARGV == 1;

my $motif = read_jaspar($ARGV[0]);
my $mlen = @$motif;
for (my $i = 0; $i < $SEQS; $i++) {
	my $seq = generate_seq($LEN, $PA, $PC, $PG, $PT);
	my @info;
	if (rand(1) < $FREQ) {
		for (my $m = 0; $m < $MPS; $m++) {
			my $site = generate_site($motif);
			my $strand = '+';
			if ($BOTH_STRANDS and rand(1) < 0.5) {
				$site =~ tr/ACGTacgt/TGCAtgca/;
				$site = reverse $site;
				$strand = '-';
			}
			my $pos = int(rand($LEN - $mlen +1));
			push @info, "$pos$strand";
			substr($seq, $pos, $mlen) = $site;
		}
	}
	print ">seq-$i motifs: @info\n", "$seq\n";
}

sub generate_seq {
	my ($len, $pa, $pc, $pg, $pt) = @_;
	my $seq = "";
	for (my $i = 0; $i < $len; $i++) {
		my $r = rand(1);
		if    ($r < $pa)         {$seq .= 'a'}
		elsif ($r < $pa+$pc)     {$seq .= 'c'}
		elsif ($r < $pa+$pc+$pg) {$seq .= 'g'}
		else                     {$seq .= 't'}
	}
	return $seq;
}

sub generate_site {
	my ($motif) = @_;
	
	my $site = "";
	for (my $i = 0; $i < @$motif; $i++) {
		my $tA = $motif->[$i]{A};
		my $tC = $tA + $motif->[$i]{C};
		my $tG = $tC + $motif->[$i]{G};
		my $r = rand(1);
		if    ($r < $tA) {$site .= 'A'}
		elsif ($r < $tC) {$site .= 'C'}
		elsif ($r < $tG) {$site .= 'G'}
		else             {$site .= 'T'}
	}
	
	return $site;
}

sub read_jaspar {
	my ($file) = @_;
	
	open(my $fh, $file) or die;
	my $head = <$fh>;
	my ($id, $name) = $head =~ /^>(\S+)\s+(\S+)/;
	my @motif;
	for (my $i = 0; $i < 4; $i++) {
		my $line = <$fh>;
		$line =~ s/[\[\]]//g;
		my ($nt, @val) = split(/\s+/, $line);
		for (my $j = 0; $j < @val; $j++) {
			$motif[$j]{$nt} = $val[$j];
		}
	}
	
	for (my $i = 0; $i < @motif; $i++) {
		my $total = 0;
		foreach my $nt (keys %{$motif[$i]}) {$total += $motif[$i]{$nt}}
		foreach my $nt (keys %{$motif[$i]}) {$motif[$i]{$nt} /= $total}
	}
	
	return \@motif;
}
