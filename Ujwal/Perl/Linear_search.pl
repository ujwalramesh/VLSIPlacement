#!/usr/bin/perl
use warnings;
use strict;

my $num = 100;
my $Max =5000;

srand();
my $i;
my @array; 
print "numbers generated :\n(";
for $i (1 .. $num ) {
	push @array, sprintf("%d",rand(1) * $Max);
	print $array[$i-1];
	print "," unless ($i == $num);
}
print ")\n\n";
my $tosearch;
chomp($tosearch = <STDIN>);
my $counter =0;
my $hit =0;
my $num1;
foreach $num1(@array) {
	$counter++;
	if($num1 == $tosearch ) {
	print "\"$tosearch \" found at subscript ",$counter -1, "\n";
	$hit =1;
	last;
	}
}

if($hit == 0 ) { print "not found \n"};


