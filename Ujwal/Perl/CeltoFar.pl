#!/usr/bin/perl
use warnings;
use strict;

# The following program converts Degree celcius to degree Farenhite 

print " Please enter temperature in degree celcius : ";

my $cel;
chomp($cel = <STDIN> );
my $fah; 
$fah =($cel *1.8 )+32;

print "the farenhite equivalent of $cel is $fah \n";


