#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;

my $flag = 1;
my $counter = 0;
open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
my $lastFit = 0;
my $fit;
my $maxDiff = 0;
while (<FILE>) {
	if (/fit ([^ ]+) /) {
		my $b = defined($fit) ? 1 : 0;
		$b = defined($lastFit) ? $b : 0;
		my $tmp = ($b) ? abs($lastFit - $fit) : 0;
		$maxDiff = $tmp if ($tmp > $maxDiff);

		$lastFit = $fit;
		$fit = $1;

	}

	if ($flag and /fit ([^ ]+) /) {
		die "$0: Value $1 not equal fit $fit\n" if (defined($fit) and $1 != $fit);

		my $value = 3-$1;
		print "$counter $value $maxDiff\n";
		++$counter;
		$flag = 0;
		$maxDiff = 0;
	}

	last if ($counter == 3000);
	if (/^\-\-\-\-/) {
		$flag = 1;
	}
}

close(FILE);

