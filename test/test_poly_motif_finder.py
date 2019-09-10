import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from pymotiffinder.motif_finder import poly_motif_finder
import pandas as pd
import numpy as np

class testMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_pmf_order(self):
        ## the example from the paper
        pmf_df_1, rate_1 = poly_motif_finder("test/dale_test_partis_1.csv",
                                   "test/dale_test_reference.fasta", k=2,
                                   max_mutation_rate=1,
                                   use_indel_seqs=True,
                                   kmer_dict=None,
                                   dale_method=True)
        pmf_df_2, rate_2 = poly_motif_finder("test/dale_test_partis_2.csv",
                                   "test/dale_test_reference.fasta", k=2,
                                   max_mutation_rate=1,
                                   use_indel_seqs=True,
                                   kmer_dict=None,
                                   dale_method=True)
        ## if ATT processed first, rate should be 1/3, if CTT processed first, rate should be 2/3
        def get_theoretical_rate(pmf_df):
            if pmf_df["query_sequence"][0] == "ATT":
                return 1./3
            if pmf_df["query_sequence"][0] == "CTT":
                return 2./3
        self.assertEqual(rate_1, get_theoretical_rate(pmf_df_1))
        self.assertEqual(rate_2, get_theoretical_rate(pmf_df_2))

    def test_pmf_no_order(self):
        pmf_df_1, rate_1 = poly_motif_finder("test/dale_test_partis_1.csv",
                                   "test/dale_test_reference.fasta", k=2,
                                   max_mutation_rate=1,
                                   use_indel_seqs=True,
                                   kmer_dict=None,
                                   dale_method=False)
        pmf_df_2, rate_2 = poly_motif_finder("test/dale_test_partis_2.csv",
                                   "test/dale_test_reference.fasta", k=2,
                                   max_mutation_rate=1,
                                   use_indel_seqs=True,
                                   kmer_dict=None,
                                   dale_method=False)
        ## in each case, 2 templated mutations, 5 mutations total
        self.assertEqual(rate_1, 2./5)
        self.assertEqual(rate_2, 2./5)

if __name__ == '__main__':
    unittest.main()
