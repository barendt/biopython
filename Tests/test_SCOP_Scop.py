# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Scop"""

import unittest
from StringIO import *

from Bio.SCOP import *

import sys

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    return test_suite()


class ScopTests(unittest.TestCase):


    def testParse(self):
  
        f = open("./SCOP/dir.cla.scop.txt_test")
        try:
            cla = f.read()
            f.close()
            
            f = open("./SCOP/dir.des.scop.txt_test")
            des = f.read()
            f.close()

            f = open("./SCOP/dir.hie.scop.txt_test")
            hie = f.read()
        finally:
            f.close()

        scop = Scop(StringIO(cla), StringIO(des), StringIO(hie))

        cla_out = StringIO()
        scop.write_cla(cla_out)
        assert cla_out.getvalue() == cla, cla_out.getvalue()
        
        des_out = StringIO()
        scop.write_des(des_out)
        assert des_out.getvalue() == des, des_out.getvalue()

        hie_out = StringIO()
        scop.write_hie(hie_out)
        assert hie_out.getvalue() == hie, hie_out.getvalue()

        domain = scop.getDomainBySid("d1hbia_")
        assert domain.sunid == '14996'

        domains = scop.getDomains()
        assert len(domains)==14
        assert domains[4].sunid == '14988'



    def testSccsOrder(self) :
        assert cmp_sccs("a.1.1.1", "a.1.1.1") == 0
        assert cmp_sccs("a.1.1.2", "a.1.1.1") == 1
        assert cmp_sccs("a.1.1.2", "a.1.1.11") == -1
        assert cmp_sccs("a.1.2.2", "a.1.1.11") == 1
        assert cmp_sccs("a.1.2.2", "a.5.1.11") == -1         
        assert cmp_sccs("b.1.2.2", "a.5.1.11") == 1
        assert cmp_sccs("b.1.2.2", "b.1.2") == 1        

    def testParseDomain(self) :
        s=">d1tpt_1 a.46.2.1 (1-70) Thymidine phosphorylase {Escherichia coli}"
        dom = parse_domain(s)

        assert dom.sid == 'd1tpt_1'
        assert dom.sccs == 'a.46.2.1'
        assert dom.residues.pdbid == '1tpt'
        assert dom.description == 'Thymidine phosphorylase {Escherichia coli}'

        s2="d1tpt_1 a.46.2.1 (1tpt 1-70) Thymidine phosphorylase {E. coli}"
        assert s2 == str(parse_domain(s2)), str(parse_domain(s2))



        #Genetic domains (See Astral release notes)
        s3="g1cph.1 g.1.1.1 (1cph B:,A:) Insulin {Cow (Bos taurus)}"
        assert s3 == str(parse_domain(s3)), str(parse_domain(s3))

        s4="e1cph.1a g.1.1.1 (1cph A:) Insulin {Cow (Bos taurus)}"
        assert s4 == str(parse_domain(s4))

        #Raw Astral header
        s5=">e1cph.1a g.1.1.1 (A:) Insulin {Cow (Bos taurus)}"
        assert s4 ==  str(parse_domain(s5))

        try:
            dom = parse_domain("Totally wrong")
            assert 0, "Should never get here"
        except SyntaxError, e :
            pass
            

def test_suite():
    return unittest.makeSuite(ScopTests)

if __name__ == '__main__':
    unittest.main()




