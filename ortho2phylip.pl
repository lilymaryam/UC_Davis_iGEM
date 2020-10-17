#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

my %d;
my %taxa;
my $test = 0;
while (<>) {
	chomp;
	my ($c1, $c2, $hit, $s1, $s2, $score) = split;
	
	my ($d1) = $c1 =~ /:(\d+)/;
	my ($d2) = $c2 =~ /:(\d+)/;
	next if not defined $d1 or not defined $d2;
	
	$d{$d1}{$d2} = 1/$score;
	$taxa{$d1} = 1;
	$taxa{$d2} = 1;
}

my $ntaxa = keys %taxa;
my @taxa = sort keys %taxa;

print "$ntaxa\n";
foreach my $t1 (@taxa) {
	print $t1;
	foreach my $t2 (@taxa) {
		my $val = defined $d{$t1}{$t2} ? $d{$t1}{$t2} : 0;
		printf "\t%.3f", $val;
	}
	print "\n";
}
   
