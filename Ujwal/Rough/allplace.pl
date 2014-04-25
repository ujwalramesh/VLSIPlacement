#!/usr/bin/perl -w
#use strict;

#### The below list of variables contain the paths required to run the script
$placerpath = "/home/tirumanl/thesis/VLSIPlacement/scripts";

#### Read the benchmark to be run from the command line ###

my $benchmarkname = $ARGV[0];

print "Placement is run for the benchmark \"$benchmarkname\" \n";

$resultspath="/home/tirumanl/results/$benchmarkname/bookshelf";
$rptsuffix="pl_post_route_timing.rpt";
$rptsuffix1="pl_pre_route_timing.rpt";

#### Array which contains all the placers to be run. Append new placers to the end of array 

my @placerlist;
push @placerlist, 'NTUPlace';
push @placerlist, 'Dragon';
push @placerlist, 'mPL';

my $placers = join(", ",@placerlist);
print "The Following placers \($placers\) are run for the benchmark $benchmarkname\n";

# The following loop runs in the sequence of steps mentioned below for each placer in the placerlist array for the input benchmark
#Step1: Run the placer
#Step2: Extract the timing details from the generated reports
#Step3: Tabulate the extracted timing details
@header = ('Benchmarkname','Placer','PreWNS','PostWNS');
@finaltable=([@header]);
foreach $plindex(@placerlist){
    print "Run details for $plindex\n";   
    (system("$placerpath/run_placement.sh $benchmarkname $plindex") == 0) || die ("Cannot run placer\n");
    if($plindex eq "Dragon"){
	my @temparray = ($benchmarkname,$plindex);
	open PRERPT, "$resultspath/$benchmarkname\_Dragon.$rptsuffix1" or die $!;
	while ($line1 = <PRERPT>){
	    if($line1 =~ /^.*?slack[\s]+\(VI.*/){
		$line1 =~ m/([0-9]+\.[0-9]+)/;
		my $temp = $1;
		push @temparray, $temp;
		last;
	    }
	}

	open POSTRPT, "$resultspath/$benchmarkname\_Dragon.$rptsuffix" or die $!; 
	while ($line2 = <POSTRPT>){
	    if($line2 =~ /^.*?slack[\s]+\(VI.*/){
		$line2 =~ m/([0-9]+\.[0-9]+)/;
		my $temp1 = $1;
		push @temparray, $temp1;
		last;
	    }
	}
	push (@finaltable,[@temparray]);
#	use Data::Dumper;
#	print Dumper(@finaltable);

    }
    elsif ($plindex eq "NTUPlace"){
	my @temparray = ($benchmarkname,$plindex);
	open PRERPT, "$resultspath/$benchmarkname.ntup.$rptsuffix1" or die $!;
	while ($line1 = <PRERPT>){
	    if($line1 =~ /^.*?slack[\s]+\(VI.*/){
		$line1 =~ m/([0-9]+\.[0-9]+)/;
		my $temp = $1;
		push @temparray, $temp;
		last;
	    }
	}
	open POSTRPT, "$resultspath/$benchmarkname.ntup.$rptsuffix" or die $!;
	while ($line2 = <POSTRPT>){
	    if($line2 =~ /^.*?slack[\s]+\(VI.*/){
		$line2 =~ m/([0-9]+\.[0-9]+)/;
		my $temp1 = $1;
		push @temparray, $temp1;
		last;
	    }
	}
	push (@finaltable,[@temparray]);
#	use Data::Dumper;
#	print Dumper(@finaltable);
    } 
    elsif($plindex eq "mPL"){
	open PRERPT, "$resultspath/$benchmarkname-mPL.$rptsuffix1" or die $!;
	open POSTRPT, "$resultspath/$benchmarkname-mPL.$rptsuffix" or die $!;
	my @temparray = ($benchmarkname,$plindex);
	open PRERPT, "$resultspath/$benchmarkname-mPL.$rptsuffix1" or die $!;
	while ($line1 = <PRERPT>){
	    if($line1 =~ /^.*?slack[\s]+\(VI.*/){
		$line1 =~ m/([0-9]+\.[0-9]+)/;
		my $temp = $1;
		push @temparray, $temp;
		last;
	    }
	}
	open POSTRPT, "$resultspath/$benchmarkname-mPL.$rptsuffix" or die $!;
	while ($line2 = <POSTRPT>){
	    if($line2 =~ /^.*?slack[\s]+\(VI.*/){
		$line2 =~ m/([0-9]+\.[0-9]+)/;
		my $temp1 = $1;
		push @temparray, $temp1;
		last;
	    }
	}
	push (@finaltable,[@temparray]);
#	use Data::Dumper;
#	print Dumper(@finaltable);
    } 

    close PRERPT;
    close POSTRPT;
######### The below step creates an html file to view the results############\
###Final table contains the entire tabulated information ######
    open FILE1, ">myfir.html" or die $!;
    print  FILE1 "Content-type: text/html \n\n";

    print  FILE1 "<html> <head>\n";
    print  FILE1 "<title>Hello, world!</title>";
    print  FILE1 "</head>\n";
    print  FILE1 "<body>\n";
    print  FILE1 "<h1>Benchmark statistics</h1>\n";
    print  FILE1 "<table border='1'>";
    for $i(0..(scalar @finaltable)-1){
	print  FILE1 "<tr>";
	for $j(0..3){
	    print FILE1 "<td>$finaltable[$i][$j]</td>";
	}
	print FILE1 "</tr>";
    }
    print FILE1 "</table>";
    print  FILE1 "</body> </html>\n";
    close FILE1;
}
exit(0);

