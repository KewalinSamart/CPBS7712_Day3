'''
Perform unit tests for essential methods
'''
import unittest
from src import binary_search_tree, compare_sequences

class unitities(unittest.TestCase):

    def test_getkmers1(self):
        input_seq = "IAMCUTE"
        k = 3
        kmer_list = list(get_kmers(input_seq, k).keys())
        expected = ["IAM","AMC","MCU","CUT","UTE"]
        self.assertEqual(kmer_list.sort(),expected.sort())

    def test_getkmers2(self):
        input_seq = "IAMCCUTE"
        k = 4
        kmer_list = list(get_kmers(input_seq, k).keys())
        expected = ["IAM","AMCC","MCCU","CCUT","CUTE"]
        self.assertEqual(kmer_list.sort(),expected.sort())
        
class alignment(unittest.TestCase):

    def test_smithwaterman_100iden(self):
        input_seq1 = "IAMCUTE"
        input_seq2 = "IAMCUTE"
        expected_iden = smith_waterman(input_seq1, input_seq2)[0]
        self.assertEqual(expected_iden, 100)

    def test_smithwaterman_less100iden(self):
        input_seq1 = "IAMCUTE"
        input_seq2 = "NOTCUTE"
        expected_iden = round(smith_waterman(input_seq1, input_seq2)[0],2)
        self.assertEqual(expected_iden, 57.14)


class assembly(unittest.TestCase):
    pass 

if __name__ == '__main__':
    unittest.main()