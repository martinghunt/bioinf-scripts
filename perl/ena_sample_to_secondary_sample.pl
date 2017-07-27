#!/usr/bin/env perl
use strict;
use warnings;

if (scalar @ARGV != 1) {
    die qq/usage: $0 sample_id

Given an ENA sample accession, returns the secondary sample accession
/;
}


my $sample_id = $ARGV[0];

my $url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$sample_id&result=read_run&fields=secondary_sample_accession&download=txt";

system("wget -q -O tmp.$$ '$url'") and die $!;

open F, "tmp.$$" or die $!;
my $line = <F>;
die unless $line eq "secondary_sample_accession\n";
my $secondary_accession = <F>;
chomp $secondary_accession;
length($secondary_accession) or die $!;
close F;
print "$sample_id\t$secondary_accession\n";
unlink "tmp.$$";
