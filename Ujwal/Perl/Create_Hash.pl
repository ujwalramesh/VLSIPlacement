#!/usr/bin/perl
use warnings;
use strict;

my %first_hash;

%first_hash = ('1' => 'one', '2' => 'two', '3' => 'three');

use Data::Dumper;
print Dumper %first_hash;

print "$first_hash{1}\n";

my @array1 =(1,2,3,4,5);
my @array2 = @array1;
$, = ",";
my $d = $array1[2] + $array1[3];
print " array values are @array1 \n ";
print " sum of elements is $d \n";
