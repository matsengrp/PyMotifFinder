import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from pymotiffinder.motif_finder import poly_motif_finder, templated_number
import pandas as pd
import numpy as np

class testMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_templated_number(self):
        pmf_out, pmf_frac = poly_motif_finder("test/test_hit_number_partis.csv",
                                    "test/dale_test_reference.fasta", k=2,
                                    max_mutation_rate=1,
                                    use_indel_seqs=True,
                                    kmer_dict=None,
                                    dale_method=True)

        (hits, num_mutations) = templated_number(pmf_out, dale_method=True)
        self.assertEqual(num_mutations, 12)
        self.assertEqual(hits, 6)
        (hits, num_mutations) = templated_number(pmf_out, dale_method=False)
        self.assertEqual(hits, 6)
        self.assertEqual(num_mutations, 14)

        pmf_out_15 = pmf_out[(pmf_out.query_name == "s1") | (pmf_out.query_name == "s5")]
        (hits, num_mutations) = templated_number(pmf_out_15, dale_method=True)
        self.assertEqual(num_mutations, 4)
        self.assertEqual(hits, 0)
        (hits, num_mutations) = templated_number(pmf_out_15, dale_method=False)
        self.assertEqual(hits, 0)
        self.assertEqual(num_mutations, 4)

        pmf_out_12 = pmf_out[(pmf_out.query_name == "s1") | (pmf_out.query_name == "s2")]
        (hits, num_mutations) = templated_number(pmf_out_12, dale_method=True)
        self.assertEqual(num_mutations, 4)
        self.assertEqual(hits, 2)
        (hits, num_mutations) = templated_number(pmf_out_12, dale_method=False)
        self.assertEqual(hits, 2)
        self.assertEqual(num_mutations, 4)

        pmf_out_145 = pmf_out[(pmf_out.query_name == "s1") | (pmf_out.query_name == "s4") | (pmf_out.query_name == "s5")]
        (hits, num_mutations) = templated_number(pmf_out_145, dale_method=True)
        self.assertEqual(num_mutations, 6)
        self.assertEqual(hits, 0)
        (hits, num_mutations) = templated_number(pmf_out_145, dale_method=False)
        self.assertEqual(hits, 0)
        self.assertEqual(num_mutations, 6)

        pmf_out_267 = pmf_out[(pmf_out.query_name == "s2") | (pmf_out.query_name == "s6") | (pmf_out.query_name == "s7")]
        (hits, num_mutations) = templated_number(pmf_out_267, dale_method=True)
        self.assertEqual(num_mutations, 6)
        self.assertEqual(hits, 6)
        (hits, num_mutations) = templated_number(pmf_out_267, dale_method=False)
        self.assertEqual(hits, 6)
        self.assertEqual(num_mutations, 6)

        pmf_out_46 = pmf_out[(pmf_out.query_name == "s4") | (pmf_out.query_name == "s6")]
        (hits, num_mutations) = templated_number(pmf_out_46, dale_method=True)
        self.assertEqual(num_mutations, 4)
        self.assertEqual(hits, 2)
        (hits, num_mutations) = templated_number(pmf_out_46, dale_method=False)
        self.assertEqual(hits, 2)
        self.assertEqual(num_mutations, 4)


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
