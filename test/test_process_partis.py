import unittest
from pymotiffinder.process_partis import get_pairs, process_partis, process_partis_poly
from pymotiffinder.motif_finder import likelihood_given_gcv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class testPP(unittest.TestCase):

    def setUp(self):
        pass

    def test_process_partis(self):
        # run partis on the test data
        partis_file = "test/partis_test.csv"
        mut_df = process_partis(partis_file)
        # there is one mutation in the test at position 7 with gl base
        # G and mutated base A
        self.assertEqual(mut_df.shape[0], 1)
        self.assertEqual(mut_df["mutated_seq"][0], "AAAAAAAA")
        self.assertEqual(mut_df["naive_seq"][0], "AAAAAAAG")
        self.assertEqual(mut_df["mutated_seq_id"][0], "s1")
        self.assertEqual(mut_df["mutation_index"][0], 7)
        self.assertEqual(mut_df["gl_base"][0], "G")
        self.assertEqual(mut_df["mutated_base"][0], "A")

    def test_process_partis_poly(self):
        # run partis on the test data
        partis_file = "test/partis_poly_test.csv"
        mut_df = process_partis_poly(partis_file, max_spacing=2)
        self.assertEqual(mut_df.shape[0], 1)
        self.assertEqual(mut_df["mutated_seq"][0], "AAAAAAAA")
        self.assertEqual(mut_df["naive_seq"][0], "GAAAAGAG")
        self.assertEqual(mut_df["mutated_seq_id"][0], "s1")
        self.assertEqual(mut_df["mutation_index"][0], (5,7))
        self.assertEqual(mut_df["gl_base"][0], "GG")
        self.assertEqual(mut_df["mutated_base"][0], "AA")

    def test_get_pairs(self):
        # set up the test data
        mut_indices = [0, 1, 8, 9, 10, 15, 18, 21, 23]
        # run get_pairs
        pairs = get_pairs(mut_indices, max_spacing=2)
        # check that 
        self.assertEqual((0, 1) in pairs, True)
        self.assertEqual((8, 9) in pairs, True)
        self.assertEqual((8, 10) in pairs, True)
        self.assertEqual((9, 10) in pairs, True)
        self.assertEqual((21, 23) in pairs, True)
        self.assertEqual(len(pairs), 5)


if __name__ == '__main__':
    unittest.main()
