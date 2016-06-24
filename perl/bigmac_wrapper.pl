#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use File::Spec;
use File::Which;
use File::Basename;
use Getopt::Long;

my %options = (
     code_root => '~mh12/bin/MetaFinisherSC-5.1/srcRefactor/',
);

my $ops_ok = GetOptions(
    \%options,
    'code_root|c=s',
    'noclean',
);

if ($#ARGV != 2 or !($ops_ok)) {
    print STDERR "usage: $0 <reads.fasta> <assemble.fasta> <output directory>

Wrapper script to run BIGMAC. Assumes nucmer is in your path.

Options:

-c, --code_root
    Path to the srcRefactor/ directory [$options{code_root}]

--noclean
    Do not clean intermediate files. By default,
    almost everthing is deleted
";
    exit(1);
}

$options{code_root} = File::Spec->rel2abs(glob($options{code_root}));
my $reads = File::Spec->rel2abs($ARGV[0]);
my $input_contigs = File::Spec->rel2abs($ARGV[1]);
my $outdir = File::Spec->rel2abs($ARGV[2]);

# find MUMmer directory
print "Finding MUMmer directory...\n";
my $nucmer = which('nucmer');
die "nucmer not found in path. Cannot continue" unless defined($nucmer);
$nucmer = readlink($nucmer);
print "nucmer after following symlinks: $nucmer\n";
my $mummer_dir = dirname($nucmer);
print "MUMmer dir: $mummer_dir\n";


# setup data directory and symlinks to run bigmac
mkdir $outdir or die $!;
chdir $outdir or die $!;
mkdir 'data' or die $!;
system_call(q~perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge'~ . " $reads > data/LR.fasta");
system_call(q~perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge'~ . " $input_contigs > data/LC.fasta");
symlink "$options{code_root}/", 'srcRefactor' or die $!;
print "\nMade symlinks in output directory $outdir\n\n";


# run bigmac
system_call("python -m srcRefactor.misassemblyFixerLib.mFixer data $mummer_dir");
system_call("python -m srcRefactor.repeatPhaserLib.aSplitter data $mummer_dir");


if ($options{noclean}) {
    print "--noclean used, so keeping all files. Symlinking to final contigs file\n";
    symlink 'data/abun.fasta', 'contigs_out.fasta' or die $!;
}
else {
    print "Tidying up temporary files\n";
    rename 'data/abun.fasta', 'contigs_out.fasta' or die $!;
    system_call("rm -fr data srcRefactor");
}

sub system_call {
    my $cmd  = shift;
    print "Running command: $cmd\n";
    if (system($cmd)) {
        print STDERR "Error in system call:\n$cmd\n";
        exit(1);
    }
    print "Command finished: $cmd\n";
    print "_____________________________________________________________________\n\n";
}

