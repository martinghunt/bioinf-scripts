#!/usr/bin/env perl
use strict;
use warnings;
use LWP::Simple;

if ($#ARGV != 2) {
    die qq/usage: genbank_downloader.pl filetype id outfile

Examples:

genbank_downloader.pl fasta GQ983346 GQ983346.fa

genbank_downloader.pl gb GQ983346 GQ983346.gb
/;
}

my $filetype=$ARGV[0];
my $id=$ARGV[1];
my $out=$ARGV[2];

# file has an empty line at the end, so need to remove it
my $file_contents = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=$filetype&retmode=text&id=$id");
chomp $file_contents;
open F, ">$out" or die $!;
print F $file_contents;
close F or die $!;
