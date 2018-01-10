#!/usr/bin/env python3

import os
import filecmp
import unittest

import pyfastaq

import bam_to_rmdup_fastq

tests_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(tests_dir, 'data', 'bam_to_rmdup_fastq')

def reads_same_in_fastq_files(file1, file2):
    d1 = {}
    d2 = {}
    pyfastaq.tasks.file_to_dict(file1, d1)
    pyfastaq.tasks.file_to_dict(file2, d2)
    return d1 == d2


def reads_in_same_order_in_fastq_files(file1, file2):
    ids1 = [seq.id.rstrip('1') for seq in pyfastaq.sequences.file_reader(file1)]
    ids2 = [seq.id.rstrip('2') for seq in pyfastaq.sequences.file_reader(file2)]
    return ids1 == ids2


class TestScript(unittest.TestCase):
    def test_whole_script(self):
        '''Test running whole script'''
        bam_in = os.path.join(data_dir, 'in.bam')
        fq_out1 = 'tmp.bam_to_rmdup_fastq.out_1.fq'
        fq_out2 = 'tmp.bam_to_rmdup_fastq.out_2.fq'
        count_out = 'tmp.bam_to_rmdup_fastq.counts.tsv'
        expected1 = os.path.join(data_dir, 'out1.fq')
        expected2 = os.path.join(data_dir, 'out2.fq')
        expected_count = os.path.join(data_dir, 'out.tsv')
        options = bam_to_rmdup_fastq.parser.parse_args([bam_in, fq_out1, fq_out2, count_out])
        bam_to_rmdup_fastq.run(options)
        self.assertTrue(filecmp.cmp(expected_count, count_out, shallow=False))
        self.assertTrue(reads_same_in_fastq_files(expected1, fq_out1))
        self.assertTrue(reads_in_same_order_in_fastq_files(fq_out1, fq_out2))
        os.unlink(fq_out1)
        os.unlink(fq_out2)
        os.unlink(count_out)


    def test_fails_on_bad_input(self):
        '''Test throws error when input BAM is bad'''
        bam_in = os.path.join(data_dir, 'in.bad.bam')
        fq_out1 = 'tmp.bam_to_rmdup_fastq.out_1.fq'
        fq_out2 = 'tmp.bam_to_rmdup_fastq.out_2.fq'
        count_out = 'tmp.bam_to_rmdup_fastq.counts.tsv'
        options = bam_to_rmdup_fastq.parser.parse_args([bam_in, fq_out1, fq_out2, count_out])
        with self.assertRaises(bam_to_rmdup_fastq.Error):
            bam_to_rmdup_fastq.run(options)

