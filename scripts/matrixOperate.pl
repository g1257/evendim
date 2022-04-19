#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $ax, $bx) = @ARGV;
defined($bx) or die "USAGE: $0 filename a b\n";

my @a = loadFile($file);

operate(\@a, $ax, $bx);

printFile(\@a);

sub operate
{
	my ($a, $ax, $bx) = @_;
	my $rows = scalar(@$a);
	my $ptr0 = $a->[0];
	my $cols = scalar(@$ptr0);
	for (my $i = 0; $i < $rows; ++$i) {
		my $ptr = $a->[$i];
		if ($cols != scalar(@$ptr)) {
			die "$0: printFile: wrong number of columns\n";
		}

		for (my $j = 0; $j < $cols; ++$j) {
			my $val = $ax*$ptr->[$j];
			$val += $bx if ($i == $j);
			$a->[$i]->[$j] = $val;
		}
	}
}

sub printFile
{
	my ($a) = @_;
	my $rows = scalar(@$a);
	my $ptr0 = $a->[0];
	my $cols = scalar(@$ptr0);
	print "$rows $cols\n";
	for (my $i = 0; $i < $rows; ++$i) {
		my $ptr = $a->[$i];
		if ($cols != scalar(@$ptr)) {
			die "$0: printFile: wrong number of columns\n";
		}

		for (my $j = 0; $j < $cols; ++$j) {
			my $val = $ptr->[$j];
			print "$val ";
		}

		print "\n";
	}
}

sub loadFile
{
	my ($file) = @_;
	open(my $fh, "<", $file) or die "$0: Cannot open $file : $!\n";
	$_ = <$fh>;
	chomp;
	my ($rows, $cols) = split;
	my @a;
	my $row = 0;
	while (<$fh>) {
		my @temp = split;
		scalar(@temp) == $cols or die "$0: Wrong number of columns for row $row\n";
		$a[$row++] = \@temp;
	}

	close($fh);
	$row == $rows or die "$0: Rows wrong $row != $rows\n";
	return @a;
}

