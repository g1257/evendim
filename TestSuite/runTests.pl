#!/usr/bin/perl

use strict;
use warnings;
use utf8;

use Term::ANSIColor;

my $testN = 100;
announce($testN);
my $output = captureStdout($testN, "../src/gep2 -i 1 -h 5 -p 20 -t 10 -e 1");
compareStdout($testN, $output);
endTest($testN);

$testN = 103;
announce($testN);
$output = captureStdout($testN, "../src/gep2 -i 6 -h 14 -p 60 -t 30 -e 3");
compareStdout($testN, $output);
endTest($testN);

sub announce
{
	my ($testNumber) = @_;
	print "************* Test Number $testNumber *********** \n";
}

sub endTest
{
	my ($testNumber) = @_;
	print "----------------------------------\n\n";
}

sub captureStdout
{
	my ($testN, $cmd) = @_;
	my $fout = "out$testN.txt";
	system("$cmd 1> $fout 2> /dev/null");
	return $fout;
}

sub compareStdout
{
	my ($testN, $testOutput) = @_;
	my $oracle = "oracles/$testN.txt";
	system("diff -q $oracle $testOutput");
	my $x = $?;
	if ($x == 0) {
		print "$0: " . color("green") . " $testN passed " . color("reset") . "\n";
	} else {
		print "$0: " . color("red") . " $testN failed" . color("reset") . "\n";
	}
}

=pod

Red: \u001b[31m
Reset: \u001b[0m

Black: \u001b[30m
Red: \u001b[31m
Green: \u001b[32m
Yellow: \u001b[33m
Blue: \u001b[34m
Magenta: \u001b[35m
Cyan: \u001b[36m
White: \u001b[37m
Reset: \u001b[0m
=cut
