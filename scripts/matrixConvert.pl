#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my @a = loadFile($file);

printFile(\@a);

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
	my ($rows, $cols);
	my @a;
	while (<$fh>) {
		next if (/^#/);
		chomp;
		if (!defined($rows)) {
			my @temp = split;
			scalar(@temp) == 2 or die "$0: Error for line $_\n";
			($rows, $cols) = @temp;
			next;
		}

		#(4, 4) (-159.58499847820988+0j)
		if (/^ *\( *(\d+) *, *(\d+) *\) *\( *([^+]+)\+0j\) *$/) {
			my $row = $1;
			my $col = $2;
			my $val = $3;
			if ($row >= $rows) {
				die "$0: Error $row >= $rows\n";
			}

			if ($col >= $cols) {
				die "$0: Error $col >= $cols\n";
			}

			my $ptr = $a[$row];
			if (!defined($ptr)) {
				my @temp;
				$temp[$col] = $val;
				$a[$row] = \@temp;
			} else {
				my $x = $a[$row]->[$col];
				if (defined($x)) {
					die "$0: Error: Redefining at ($row, $col) with $val\n";
				}

				$a[$row]->[$col] = $val;
			}
		} else {
			die "$0: ERROR: FILE $file and Line $_\n";
		}
	}

	close($fh);
	fillZeroes(\@a, $rows, $cols);
	return @a;
}

sub fillZeroes
{
	my ($ptr, $rows, $cols) = @_;
	for (my $i = 0; $i < $rows; ++$i) {
		my $ptr2 = $ptr->[$i];
		for (my $j = 0; $j < $cols; ++$j) {
			my $val = $ptr2->[$j];
			next if (defined($val));
			$ptr->[$i]->[$j] = 0;
		}
	}
}

