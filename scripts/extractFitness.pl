#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;

my $flag = 1;
my $counter = 0;
open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";

while (<FILE>) {
	if ($flag and /fitness ([^ ]+) /) {
		my $value = -$1;
		print "$counter $value\n";
		++$counter;
		#$flag = 0;
	}

	last if ($counter == 3000);
	if (/^\-\-\-\-/) {
		$flag = 1;
	}
}

close(FILE);

