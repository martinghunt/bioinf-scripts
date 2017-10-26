#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 3 or die "Usage: $0 <in.bam> <out1.fq> <out2.fq>

Converts paired BAM to pair of FASTQ files. Checks
read counts before/afeter. If output files end with
.gz, will gzip them.
Assumes fqtools and samtools installed.
";


my $bam_in = $ARGV[0];
my $out1 = $ARGV[1];
my $out2 = $ARGV[2];

-e $bam_in or die "Input BAM $bam_in not found";
my $samtools_cmd = "samtools fastq -N -1 $out1 -2 $out2 $bam_in";
system($samtools_cmd) and die "Error running samtools: $samtools_cmd";

my $expected_reads = `samtools view $bam_in | wc -l`;
chomp $expected_reads;
die "Error counting reads in BAM" unless $expected_reads =~ /^\d+$/;

my $got_reads = `fqtools count $out1 $out2`;
chomp $got_reads;
die "Error counting reads in output files" unless $got_reads =~ /^\d+$/;
$got_reads *= 2; # it reports the number of reads in each file

$expected_reads == $got_reads or die "Mismatch in read counts";
gzip($out1);
gzip($out2);

sub gzip {
    my $file = shift;
    if ($file =~ /\.gz$/) {
        system("gzip -9 -c $file > $file.$$") and die "Error gzip $file";
        rename("$file.$$", "$file");
    }
}
