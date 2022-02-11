#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($lx, $ly) = @ARGV;
my @links = (defined($ly)) ? heisenberg2D($lx, $ly) : heisenberg1D($lx);

printHamiltonian(\@links);

sub printHamiltonian
{
	my ($links) = @_;
	my $firstCall = 1;
	for my $link (@$links) {
		printOneHlink($link, $firstCall);
		$firstCall = 0;
	}

	print "\n";
}

sub printOneHlink
{
	my ($link, $firstCall) = @_;
	my ($a, $b) = @$link;
	my @letters = ('x', 'y', 'z');
	foreach my $letter (@letters) {
		print "+" unless ($firstCall);
		$firstCall = 0;
		print "S${letter}$a*S${letter}$b";
	}
}

sub printHeader
{
	print "$lx";
	print " x $ly" if defined($ly);
	print "\n";
}

sub printLinks
{
	my ($links) = @_;
	for my $link (@$links) {
		printLink($link);
		print " ";
	}

	print "\n";
}

sub printLink
{
	my ($link) = @_;
	my @temp = @$link;
	print "$temp[0]:$temp[1]";
}

sub heisenberg1D
{
	my ($lx) = @_;
	defined($lx) or die "USAGE: $0 lx [ly]\n";
	my $isPeriodic = getPeriodic(\$lx);

	my @links;
	for (my $i = 0; $i < $lx; ++$i) {
		my $j = $i + 1;
		if ($j == $lx) {
			last if (!$isPeriodic);
			$j = 0;
		}

		die "$0: INTERNAL ERROR $j\n" if ($j >= $lx);

		my @temp = getMinMax($i, $j);
		push @links, \@temp;
	}

	return @links;
}

sub heisenberg2D
{
	my ($lx, $ly) = @_;
	my $isPeriodicX = getPeriodic(\$lx);
	my $isPeriodicY = getPeriodic(\$ly);

	for (my $i = 0; $i < $lx; ++$i) {
		my $ii = $i + 1;
		my $flagX = 1;
		if ($ii == $lx) {
			$flagX = 0 if (!$isPeriodicX);
			$ii = 0;
		}

		for (my $j = 0; $j < $ly; ++$j) {
			my $site1 = getSite($i, $j, $lx, $ly);

			my $jj = $j + 1;
			my $flagY = 1;
			if ($jj == $ly) {
				$flagY = 0 if (!$isPeriodicY);
				$jj = 0;
			}

			createLink(\@links, $site1, getSite($ii, $j, $lx, $ly)) if ($flagX);
			createLink(\@links, $site1, getSite($i, $jj, $lx, $ly)) if ($flagY);
		}
	}

	return @links;
}

sub getSite
{
	my ($ind, $jnd, $lx, $ly) = @_;
	die "$0: INTERNAL ERRORx $ind\n" if ($ind >= $lx);
	die "$0: INTERNAL ERRORy $jnd\n" if ($jnd >= $ly);
	return $ind + $lx*$jnd;
}

sub createLink
{
	my ($links, $site1, $site2) = @_;
	my @temp = getMinMax($site1, $site2);
	push @$links, \@temp;
}

sub getPeriodic
{
	my ($lx) = @_;
	if ($$lx =~ s/p$//) {
		return 1;
	}

	return 0;
}

sub getMinMax
{
	my ($ind, $jnd) = @_;
	die "$0: getMinMax $ind == $jnd\n" if ($ind == $jnd);
	my @temp;
	$temp[0] = ($ind < $jnd) ? $ind : $jnd;
	$temp[1] = ($ind < $jnd) ? $jnd : $ind;
	return @temp;
}

