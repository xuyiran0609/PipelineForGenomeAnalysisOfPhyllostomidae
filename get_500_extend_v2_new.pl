#! user/bin/perl -w
use strict;

open (DB, $ARGV[1])||die;

my %genome;
#my $scaffold = '';
my $flag = 0;
while (<DB>){
chomp;
$flag ++;
print "$flag\n";
	my $sca = $_;
	$sca =~ s/\s+$//;
	$sca =~ s/^>//g;
	$/ = ">";
	my $seq = <DB>;
	$/ = "\n";
	$seq =~ s/>$//;
	$seq =~ s/\s//g;
	$genome{$sca} = $seq;
}
close DB;
print "The fasta file has been read over!\n";

my @keys = keys %genome;
open (FH, $ARGV[0])||die;
open (OUT, ">$ARGV[2]")||die;

while (my $line = <FH>){
chomp $line;
my $title = '';
my $seq = '';
my $length = 0;
my $large = 0;
my $small = 0;
my $start = 0;
my $end = 0;
my $start1 = 0;
my $end1 = 0;
	my ($gi, $up, $down, $evalue) = split /\t/, $line;
	if ($up > $down){
		$large = $up;
		$small = $down;
	}
	else {
		$large = $down;
		$small = $up;
	}
	foreach my $key (@keys){
		if ($key =~ /^$gi\s+/ || $key =~ /^$gi$/){
			$title = $key;
			$seq = $genome{$key};
			$length = length ($seq);

			last;
		}
	}
	if ($small > 500){
		$start = $small - 500;
		$start1 = 501;
	}
	if ($small <= 500){
		$start = 1;
		$start1 = $small;
	}
	if ($length > $large + 500){
		$end = $large + 500;
	}
	if ($length <= $large + 500){
		$end = $length;
	}
	$end1 = $start1 + $large - $small;
	my $len = $end - $start + 1;
print $len."\n";
print OUT ">gi=$gi\t$start1\t$end1\t$len=$up\t$down\n", substr ($seq, ($start-1), $len), "\n";
}
close OUT;		
