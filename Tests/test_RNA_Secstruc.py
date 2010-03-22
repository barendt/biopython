#!/usr/bin/env python

from unittest import TestCase, main
from RNA.Secstruc.Secstruc import Secstruc, PseudoknotSecstruc

class SecstrucTests(TestCase):

    def setUp(self):
        self.hairpin = Secstruc(HAIRPIN_SS)
        self.bulge = Secstruc(BULGE_SS)
        self.junction = Secstruc(JUNCTION_SS)
        self.helix = Secstruc("((()))",[0,1,2,6,7,8])
        self.helix_shifted = Secstruc("((()))",[2,3,4,8,9,10])
        self.helix_a = Secstruc("(((",[0,1,2])
        self.helix_b = Secstruc(")))",[6,7,8])

    def test_eq(self):
        self.assertEqual(self.helix,Secstruc("((()))",[0,1,2,6,7,8]))
        self.assertNotEqual(self.helix, self.helix_a)
        self.assertNotEqual(self.helix_a, self.helix_b)
        self.assertNotEqual(self.helix, self.helix_shifted)
        self.assertEqual(self.hairpin, "((((...))))")

    def test_add(self):
        added = self.helix_a + self.helix_b
        self.assertEqual(added, self.helix)
        self.assertNotEqual(added, self.helix_shifted)

    def test_slice(self):
        a = self.helix[:3]
        self.assertEqual(a, self.helix_a)
        b = self.helix[3:]
        self.assertEqual(b, self.helix_b)
        self.assertEqual(self.helix[2:4], Secstruc("()",[2,6]))

    def test_getitem(self):
        self.assertEqual(self.helix[2], Secstruc("(",[2]))
        self.assertEqual(self.helix[3], Secstruc(")",[6]))
        
    def test_get_type(self):
        """Should recognize substructures."""
        self.assertEqual(Secstruc('()').get_type(),None)
        self.assertEqual(Secstruc('((()))').get_type(),'helix')
        self.assertEqual(Secstruc('((())))').get_type(),None)
        self.assertEqual(Secstruc('(())').get_type(),'helix')        
        self.assertEqual(Secstruc('...(').get_type(),"5'-overhang")
        self.assertEqual(Secstruc(').').get_type(),"3'-overhang")
        self.assertEqual(Secstruc('(....)').get_type(),'loop')
        self.assertEqual(Secstruc('(().)').get_type(),'bulge')
        self.assertEqual(Secstruc('(.())').get_type(),'bulge')
        self.assertEqual(Secstruc('(..())').get_type(),'bulge')
        self.assertEqual(Secstruc('(..().())').get_type(),'junction')
        self.assertEqual(Secstruc('(()()())').get_type(),'junction')

    def test_find_overhang5(self):
        """Should part the object into an overhang Secstruc."""
        # 5'-overhangs
        overhang = Secstruc("..((()))").find_overhang5()
        self.assertEqual(overhang, Secstruc('..(',[0,1,2]))

    def test_find_overhang3(self):
        """Should part the object into an overhang Secstruc."""
        # 3'-overhangs
        overhang = Secstruc("(((..)))....").find_overhang3()
        self.assertEqual(overhang, Secstruc(')....',[7,8,9,10,11]))
        
    def test_find_helices(self):
        """Should return all helices from a Secstruc."""
        # simple helix
        helices = Secstruc("...(())..").find_helices()
        self.assertEqual(len(helices),1)
        self.assertEqual(helices[0], Secstruc('(())', [3,4,5,6]))
        
        # multiple helices
        helices = Secstruc(".(())...((..)).((()))()..").find_helices()
        self.assertEqual(len(helices),4)
        
        # full helix in one.
        helices = Secstruc("((()))").find_helices()
        self.assertEqual(helices[0], Secstruc('((()))',[0,1,2,3,4,5]))
        
        # two helices with bulge in between
        helices = Secstruc("((((((..))).)))").find_helices()
        self.assertEqual(len(helices),2)
        self.assertEqual(helices[0], Secstruc('((()))',[0,1,2,12,13,14]))
        self.assertEqual(helices[1], Secstruc('((()))',[3,4,5,8,9,10]))

    def test_find_loops(self):
        """Should return all loops from a Secstruc."""
        loops = Secstruc("(((..))).").find_loops()
        self.assertEqual(loops[0], Secstruc('(..)',[2,3,4,5]))

    def test_find_junctions(self):
        """Should return all junctions from a Secstruc."""
        junctions = Secstruc('(()()())').find_junctions()
        self.assertEqual(junctions[0], Secstruc('(()()())'))
        # helices should be no junctions
        junctions = Secstruc('((()))').find_junctions()
        self.assertEqual(len(junctions), 0)

    def test_find_bulges(self):
        """Should return all bulges from a Secstruc."""
        bulges = Secstruc(BULGE_SS).find_bulges()
        self.assertEqual(len(bulges),1)
        self.assertEqual(bulges[0], Secstruc('(().)',[3,4,11,12,13]))
        # also check bulge in other direction
        bulges = Secstruc("((.((..))))").find_bulges()
        self.assertEqual(bulges[0], Secstruc('(.())',[1,2,3,8,9]))
        
        
    def test_find_Secstruc_elements(self):
        """Should divide a secondary structure into substructures and their indices."""
        result = Secstruc(HAIRPIN_SS).find_Secstruc_elements()
        self.assertTrue(Secstruc('(...)',[3,4,5,6,7]) in result)
        self.assertTrue(Secstruc('(((())))',[0,1,2,3,7,8,9,10]) in result)
        self.assertEqual(len(result), 2)

    def test_find_Secstruc_elements_bulge(self):
        result = Secstruc(BULGE_SS).find_Secstruc_elements()
        self.assertTrue(Secstruc('(....)',[5,6,7,8,9,10]) in result)
        self.assertTrue(Secstruc('(().)',[3,4,11,12,13]) in result)
        self.assertTrue(Secstruc('(())',[4,5,10,11]) in result)
        self.assertTrue(Secstruc('(((())))',[0,1,2,3,13,14,15,16]) in result)
        self.assertEqual(len(result), 4)

    def test_find_Secstruc_elements_junction(self):
        result = Secstruc(JUNCTION_SS).find_Secstruc_elements()
        self.assertTrue(Secstruc('..(',[0,1,2]) in result)
        self.assertTrue(Secstruc('((()))',[2,3,4,27,28,29]) in result)
        self.assertTrue(Secstruc('(..().())',[4,5,6,7,17,18,19,26,27]) in result)
        self.assertTrue(Secstruc('(((())))',[7,8,9,10,14,15,16,17]) in result)
        self.assertTrue(Secstruc('(...)',[10,11,12,13,14]) in result)
        self.assertTrue(Secstruc('(())',[19,20,25,26]) in result)
        self.assertTrue(Secstruc('(....)',[20,21,22,23,24,25]) in result)
        self.assertEqual(len(result), 7)
        
    def test_get_bp_indices(self):
        ss = Secstruc('((()))')
        bplist = [bp for bp in ss.get_base_pair_indices()]
        bplist.sort()
        self.assertEqual(bplist, [(0, 5), (1, 4), (2, 3)])
        
    
class PseudoknotSecstrucTests(TestCase):
    
    def test_get_bp_indices(self):
        ss = PseudoknotSecstruc(PSEUDOKNOT_SS)
        bplist = [bp for bp in ss.get_base_pair_indices()]
        bplist.sort()
        self.assertEqual(bplist, [(3, 16), (5, 13), (6, 12), (7, 15)])
            
    def test_get_pseudoknot_indices(self):
        ss = PseudoknotSecstruc(PSEUDOKNOT_SS)
        bplist = [bp for bp in ss.get_pseudoknot_indices()]
        bplist.sort()
        self.assertEqual(bplist, [(5, 13), (6, 12), (7, 15)])
        
        
        
HAIRPIN_SEQ = "AAAACCCUUUU"
HAIRPIN_SS  = "((((...))))"

BULGE_SEQ = "GGCCGGAAAACCUGGCC"
BULGE_SS  = "((((((....)).))))"

JUNCTION_SEQ = "AAUCGAAUAUAGCCUAUAGCCAUGCGGCGA"
JUNCTION_SS  = "..(((..((((...)))).((....)))))"

PSEUDOKNOT_SS = '...(.(({....)).}).'

if  __name__ == "__main__":
    main()

