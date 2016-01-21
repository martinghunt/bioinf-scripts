#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use File::Which;
use File::Basename;
use File::Spec;

my %options = (
    min_len_for_query => 500,
    memory => 8000,
    threads => 1,
);

my $options_ok = GetOptions(\%options,
    'min_len_for_query=i',
    'memory=i',
    'threads=i',
);

if (($#ARGV != 4) or !($options_ok)) {
    print STDERR "usage: $0 [options] <celera bin/> <sprai bin/> <genome size> <reads> <output dir>

Runs SPRAI in single-node mode.
Assumes blastn is in your path and that makeblastdb is in the same directory as blastn

Options:

--min_len_for_query
    'the subreads longer than or equal to this value will be corrected' [$options{min_len_for_query}]

--memory
    Memory in MB that is passed to some celera options in the spec file [$options{memory}]

--threads
    Number of threads to use [$options{threads}]
";
    exit 1;
}

my $celera_bin = File::Spec->rel2abs($ARGV[0]);
my $sprai_bin = File::Spec->rel2abs($ARGV[1]);
my $genome_size = $ARGV[2];
my $reads = File::Spec->rel2abs($ARGV[3]);
my $outdir = $ARGV[4];
my $ec_spec = "ec.spec";
my $pbasm_spec = "pbasm.spec";

my $blastn = which 'blastn';
my($filename, $blastn_dir, $suffix) = fileparse($blastn);

mkdir $outdir or die $!;
chdir $outdir or die $!;

open F, ">$ec_spec" or die $!;
print F "
#### common ####
# input_for_database: filtered subreads in fasta or fastq format
input_for_database $reads

# min_len_for_query: the subreads longer than or equal to this value will be corrected
min_len_for_query $options{min_len_for_query}

#if you don't know the estimated genome size, give a large number
estimated_genome_size $genome_size
#if you don't know the estimated depth of coverage, give 0
estimated_depth 0

# ca_path: where Celera Assembler exist in
ca_path $celera_bin

# the number of processes used by all vs. all alignment
# = 'partition' (in single node mode)
# = 'pre_partition' * 'partition' (in many node mode)
pre_partition 2
partition $options{threads}

# sprai prefer full paths
# if you use ezez4qsub*.pl. you MUST specify blast_path & sprai_path
# blast_path: where blastn and makeblastdb exist in
blast_path $blastn_dir
# sprai_path: where binaries of sprai (bfmt72s, nss2v_v3 and so on) exist in
sprai_path $sprai_bin

#### many node mode (advanced) ####

#sge: options for all the SGE jobs
#sge -soft -l ljob,lmem,sjob
#queue_req: additional options for all the SGE jobs
#queue_req -l s_vmem=4G -l mem_req=4
#longestXx_queue_req: if valid, displaces queue_req
#longestXx_queue_req -l s_vmem=64G -l mem_req=64
#BLAST_RREQ: additional options for SGE jobs of all vs. all alignment
#BLAST_RREQ -pe def_slot 4

#### common (advanced) ####

# used by blastn
word_size 18
evalue 1e-50
num_threads 1
max_target_seqs 100

#valid_voters 11

#trim: both ends of each alignment by blastn will be trimmed 'trim' bases to detect chimeric reads
trim 42
";

close F or die $!;





open F, ">$pbasm_spec" or die $!;
print F "# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
unitigger = bogart
#utgErrorRate = 0.015
#utgErrorLimit = 4.5

cnsErrorRate = 0.25
cgwErrorRate = 0.25
ovlErrorRate = 0.015

frgMinLen = 1000
ovlMinLen = 40

merSize=14

merylMemory = $options{memory}
merylThreads = $options{threads}

ovlStoreMemory = $options{memory}

# grid info
useGrid = 0
scriptOnGrid = 0
frgCorrOnGrid = 0
ovlCorrOnGrid = 0

sge = -S /bin/bash -V -q all.q
#sge = -S /bin/bash -sync y -V -q all.q
sgeScript = -pe threads 1
sgeConsensus = -pe threads 1
sgeOverlap = -pe threads 4
sgeFragmentCorrection = -pe threads 4
sgeOverlapCorrection = -pe threads 1

#ovlHashBits = 22
#ovlHashBlockLength = 46871347
#ovlRefBlockSize =  537

ovlHashBits = 25
ovlThreads = $options{threads}
ovlHashBlockLength = 50000000
ovlRefBlockSize =  100000000

ovlConcurrency = $options{threads}
frgCorrThreads = $options{threads}
frgCorrBatchSize = 100000
ovlCorrBatchSize = 100000

cnsMinFrags = 7500
cnsConcurrency = $options{threads}

# change sgeName every time if you do not want to wait for the jobs not necessary to wait
sgeName = iroha
";
close F;


my $script = 'run.sh';
open F, ">$script" or die $!;
print F "#!/usr/bin/env bash
set -e
export PATH=$sprai_bin:\$PATH
ezez_vx1.pl $ec_spec $pbasm_spec
cp -s result_*/CA/9-terminator/asm_*.scf.fasta scaffolds.fasta
";
close F;

chmod 0755, $script;
my $cmd = "./$script";
exec($cmd) or die "error running $cmd: $!";

