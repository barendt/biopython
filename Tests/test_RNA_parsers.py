#!/usr/bin/env python

from unittest import TestCase, main
from Secstruc import Secstruc
from SecstrucParsers import parse_vienna, parse_bpseq, parse_ct

VIENNA = open('test_data/sample.vienna').readlines()
BPSEQ = open('test_data/sample.bpseq').readlines()
CT = open('test_data/sample.ct').readlines()

class ViennaParserTests(TestCase):
    
    def test_parse_vienna(self):
        """Should return a sequence + secstruc."""
        seq, ss = parse_vienna(VIENNA)
        self.assertEqual(len(seq), 274)
        self.assertEqual(len(ss), 274)
        self.assertEqual(seq, 'GGGAUCCCCCACAAUCCUGUCGUUACCUGUCAUGUAUCCGUCUAGAGAUCUCCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGAAGCUUGAUCCGAGGAAAAUGCUGUCCAUUAGACUAUUGAGUGCAUCUAG')
        self.assertEqual(ss, '....(((.....................((((.((.(((((......((((((((..((.....(((((....(((.(((........))).)))))))).((((..(((.(((((.(((.(((((.......))))).)))........(((......)))..)))))...))))))).))))))).)))....(((.((........)).))).....))))).)).))))....)))..((((..((.((......)).))..))))....')
        

class BPSeqParserTests(TestCase):
    
    def test_parse_bpseq(self):
        """Should return a sequence + secstruc."""
        #TODO: parsing of pknots doesnt work yet!
        seq, ss = parse_bpseq(BPSEQ)
        self.assertEqual(len(seq), 274)
        self.assertEqual(len(ss), 274)
        self.assertEqual(seq, 'GGGAUCCCCCACAAUCCUGUCGUUACCUGUCAUGUAUCCGUCUAGAGAUCUCCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGAAGCUUGAUCCGAGGAAAAUGCUGUCCAUUAGACUAUUGAGUGCAUCUAG')
        self.assertEqual(ss, '(((.((()))..................((((.((.(((((......((((((((..((.....(((((....(((.(((........))).)))))))).((((..(((.(((((.(((.(((((.......))))).)))........(((......)))..)))))...))))))).))))))).)))....(((.((........)).))).....))))).)).))))....)))..((((..((.((......)).))..))))....')
    

class CTParserTests(TestCase):
    
    def test_parse_ct(self):
        """Should return a sequence + secstruc."""
        #TODO: parsing of pknots doesnt work yet!
        seq, ss = parse_ct(CT)
        self.assertEqual(len(seq), 274)
        self.assertEqual(len(ss), 274)
        self.assertEqual(seq, 'GGGAUCCCCCACAAUCCUGUCGUUACCUGUCAUGUAUCCGUCUAGAGAUCUCCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGAAGCUUGAUCCGAGGAAAAUGCUGUCCAUUAGACUAUUGAGUGCAUCUAG')
        self.assertEqual(ss, '(((.((()))..................((((.((.(((((......((((((((..((.....(((((....(((.(((........))).)))))))).((((..(((.(((((.(((.(((((.......))))).)))........(((......)))..)))))...))))))).))))))).)))....(((.((........)).))).....))))).)).))))....)))..((((..((.((......)).))..))))....')
    
    
    
    
if __name__ == '__main__':
    main()
