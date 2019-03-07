#!/usr/bin/env perl

use strict;
use warnings;


@ARGV == 2 or die "usage: $0 <in.bam> <out.bed>\nOutputs BED file of low quality regions (where no read coverage or most common nucleotide is <95\% of total depth at a position)";

my $bam = $ARGV[0];
my $outfile = $ARGV[1];

open my $f_in, "samtools mpileup -aa $bam |" or die $!;
open my $f_out, "> $outfile" or die $!;


my $current_contig = "";
my $current_start = -1;
my $current_end = -1;
my $i = 0;
print STDERR "Gathering depths from mpileup.\n";

while (<$f_in>) {
    my @a = split;
    print STDERR "$a[0] $a[1]\n" if ($a[1] == 1 or $a[1] % 100000 == 0);
    my $is_good = pileup_array_line_is_good(\@a);
    if ($is_good and $current_contig eq "") {
        next;
    }
    elsif ($is_good and $current_contig ne "") {
        # we've come to the end of the bad region, so print it
        print $f_out "$current_contig\t" . ($current_start - 1) . "\t$current_end\n";
        $current_contig = "";
        $current_start = -1;
        $current_end = -1;
    }
    # if we're here, the position at this line is not good.
    # What we do depends on if this bad position is immediately
    # after the previous bad position
    elsif ($current_contig eq $a[0] and $current_end + 1 == $a[1]) {
        $current_end++;
    }
    else {
        unless ($current_contig eq "") {
            print $f_out "$current_contig\t" . ($current_start - 1) . "\t$current_end\n";
        }
        $current_contig = $a[0];
        $current_start = $a[1];
        $current_end = $a[1];
    }
}

if ($current_contig ne "") {
    print $f_out "$current_contig\t" . ($current_start - 1) . "\t$current_end\n";
}

print STDERR "finished\n";
close $f_in or die $!;
close $f_out or die $!;


sub pileup_array_line_is_good {
    my $array = shift;
    return 0 if ($array->[3] == 0);
    my %nuc_counts = (A => 0, G => 0, C => 0, T => 0);
    my $total = 0;
    for my $n (split(//, $array->[4])) {
        if (exists $nuc_counts{uc $n}) {
            $nuc_counts{uc $n}++;
            $total++;
        }
    }
    return 0 if $total < 1;
    return 1 if 1 == scalar keys %nuc_counts;
    my @counts = reverse sort {$a <=> $b} values %nuc_counts;
    return ($counts[0] / $total) >= 0.95;
}

