"""Builds that differ only in compiler flags must behave the same.

build/phyml-ndebug is compiled with -DNDEBUG, as distribution packages often
are. The scalar build (phyml-generic) is compared in test_likelihood_oracle.
"""

import unittest

from phymltest import PHYML_NDEBUG, PhymlTestCase, known_bug


@unittest.skipIf(PHYML_NDEBUG is None, "NDEBUG build (phyml-ndebug) not available")
class NoAssertBuild(PhymlTestCase):

    def test_interleaved_input_same_result(self):
        self.stage("sim_jc_8.phy")
        args = ["-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0", "--r_seed", "1"]
        a = self.run_phyml(*args)
        b = self.run_phyml(*args, binary=PHYML_NDEBUG, timeout=60)
        self.assertSuccess(a)
        self.assertSuccess(b)
        self.assertEqual(a.final_lnl(), b.final_lnl())

    @known_bug("io.c:1504 the sequential reader is the argument of assert(), "
               "so -DNDEBUG never reads the sequences and loops forever")
    def test_sequential_input(self):
        self.stage("oracle_dna.phy")
        res = self.run_phyml("-i", "oracle_dna.phy", "-q", "-m", "JC69", "-b", "0",
                             binary=PHYML_NDEBUG, timeout=20)
        self.assertSuccess(res)


if __name__ == "__main__":
    unittest.main()
