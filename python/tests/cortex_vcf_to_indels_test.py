#!/usr/bin/env python3

import os
import filecmp
import unittest
import cortex_vcf_to_indels

tests_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(tests_dir, 'data', 'cortex_vcf_to_indels')


class TestScript(unittest.TestCase):
    def test_whole_script(self):
        '''Test running whole script'''
        genes_in = os.path.join(data_dir, 'genes.txt')
        genbank_in = os.path.join(data_dir, 'ref.gb')
        vcf_in = os.path.join(data_dir, 'in.vcf')
        tmp_out = 'tmp.cortex_vcf_to_indels.out'
        options = cortex_vcf_to_indels.parser.parse_args([genes_in, genbank_in, vcf_in, tmp_out])
        options.add_upstream = 10
        cortex_vcf_to_indels.run(options)
        expected = os.path.join(data_dir, 'expected.tsv')
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)

