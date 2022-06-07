#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $label, $bits, $electrons) = @ARGV;
defined($file) or die "USAGE: $0 filename [label] [bits] [electrons]\n";
defined($label) or $label = "FCI_list_fc";
defined($bits) or $bits = 8;
defined($electrons) or $electrons = int($bits/2);

my $basis = readBasis($file, $label);

#print "$basis\n";

my @basis2 = stringToBasis($basis);

print STDERR "$0: Basis2 has ".scalar(@basis2)." elements.\n";

#my @basis3 = basisToConversion(\@basis2);

#print STDERR "$0: Basis3 has ".scalar(@basis3)." elements.\n";

printBasis(\@basis2);

sub basisToConversion
{
	my ($array) = @_;
	my @basis;
	my $hilbert = (1 << $bits);
	for (my $i = 0; $i < $hilbert; ++$i) {
		$basis[$i] = -1;
	}

	for (my $i = 0; $i < $hilbert; ++$i) {
		my $index = $array->[$i];
		defined($index) and $basis[$index] = $i;
	}

	return @basis;
}

sub printBasis
{
	my ($array) = @_;
	print "Basis=[".join(', ', @$array)."]\n";
}

sub stringToBasis
{
	my ($str) = @_;
	$str =~ s/\[//g;
	$str =~ s/\]//g;
	my @temp = split/,/, $str;
	my $n = scalar(@temp);
	my $r = $n % $electrons;
	$r == 0 or die "$0: $n does not divide $electrons\n";
	my $d = int($n/$electrons);
	my @basis;
	for (my $i = 0; $i < $d; ++$i) {
		$basis[$i] = convertFromBits(\@temp, $electrons*$i, $electrons);
	}

	return @basis;
}

sub convertFromBits
{
	my ($array, $start, $length) = @_;

	my @unique;

	my $sum = 0;
	for (my $i = 0; $i < $length; ++$i) {
		my $bitNumber = $array->[$start + $i];
		defined($bitNumber) or die "$0: Undefined at start $start\n";
		$bitNumber < $bits or die "$0: $bitNumber > $bits\n";
		$bitNumber < 0 and die "$0: Negative bit number $bitNumber\n";
		defined($unique[$bitNumber]) and die "$0: Repeated bit number $bitNumber\n";
		my $val = (1 << $bitNumber);
		$unique[$bitNumber] = 1;
		$sum += $val;
	}

	return $sum;
}



sub readBasis
{
	my ($file, $label) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";

	my $buffer;
	my $flag = 0;
	while (<FILE>) {
		chomp;
		if (/$label *= *(.*$)/) {
			$flag = 1;
			$buffer = $1;
			next;
		}

		if ($flag) {
			$buffer .= $_;
		}

		if (/\] *\] *$/) {
			$flag = 0;
		}
	}

	close(FILE);
	return $buffer;
}
