#!/usr/bin/env bash
set -e

if [ $# -ne 4 ]
then
    echo "usage: $0 assembly.fa reads1.fq reads2.fq outprefix

Makes a BED file of bad regions, by mapping reads and running bam_to_low_qual_mask.pl.
BWA and samtools faidx indexes the assembly fasta if needed"
    exit
fi

assembly=$1
reads1=$2
reads2=$3
outprefix=$4
bam=$outprefix.bam
bed=$outprefix.mask.bed

if [ ! -f $assembly.fai ]; then
    samtools faidx $assembly
fi

if [ ! -f $assembly.bwt ]; then
    bwa index $assembly
fi

bwa mem -x intractg $assembly $reads1 $reads2 | samtools sort -o $bam
bam_to_low_qual_mask.pl $bam $bed

