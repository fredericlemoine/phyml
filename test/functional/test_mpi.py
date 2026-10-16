"""phyml-mpi (bootstrap replicates spread over MPI ranks).

Skipped unless PHYML_MPI names an MPI build of PhyML; `make -C test check-mpi`
builds one with mpicc and sets it. MPIRUN (default: mpirun) launches it.
"""

import os
import shlex
import shutil
import unittest

from newick import transfer_support
from phymltest import PhymlTestCase, binary, known_bug, run, stage_files

PHYML_MPI = binary("PHYML_MPI", "phyml-mpi")
MPIRUN = shlex.split(os.environ.get("MPIRUN", "mpirun"))
N_RANKS = 2
N_BOOT = 10


@unittest.skipIf(PHYML_MPI is None or shutil.which(MPIRUN[0]) is None,
                 "MPI build (PHYML_MPI) or mpirun not available")
class MpiBootstrap(PhymlTestCase):

    @classmethod
    def mpi(cls, *args, ranks=N_RANKS):
        return run(MPIRUN[1:] + ["-np", str(ranks), PHYML_MPI] + list(args), cls.wd,
                   exe=shutil.which(MPIRUN[0]), timeout=600)

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.wd = cls.class_workdir()
        stage_files(cls.wd, "dna_small.phy")
        mpi = cls.mpi
        cls.fbp = mpi("-i", "dna_small.phy", "-b", str(N_BOOT), "--r_seed", "4", "--run_id", "fbp")
        cls.tbe = mpi("-i", "dna_small.phy", "-b", str(N_BOOT), "--tbe", "--r_seed", "4",
                      "--run_id", "tbe")

    def test_runs(self):
        self.assertSuccess(self.fbp)
        self.assertSuccess(self.tbe)

    def test_all_replicates_are_collected(self):
        self.assertSuccess(self.fbp)
        self.assertEqual(len(self.fbp.trees("dna_small.phy_phyml_boot_trees_fbp.txt")), N_BOOT)

    def test_fbp_counts_match_replicate_trees(self):
        self.assertSuccess(self.fbp)
        main = self.fbp.tree("dna_small.phy")
        reps = self.fbp.trees("dna_small.phy_phyml_boot_trees_fbp.txt")
        for split, v in main.supports().items():
            self.assertEqual(int(v), sum(1 for t in reps if split in t.splits()))

    def test_tbe_matches_replicate_trees(self):
        self.assertSuccess(self.tbe)
        main = self.tbe.tree("dna_small.phy")
        reps = self.tbe.trees("dna_small.phy_phyml_boot_trees_tbe.txt")
        expected = transfer_support(main, reps)
        for split, value in main.supports().items():
            self.assertAlmostEqual(value, expected[split], places=5, msg=sorted(split))

    @known_bug("mpi_boot.c:76 the replicate count is silently rounded up to a multiple of the ranks")
    def test_replicate_count_is_kept(self):
        res = self.mpi("-i", "dna_small.phy", "-b", "10", "--run_id", "r3", ranks=3)
        self.assertSuccess(res)
        self.assertEqual(len(res.trees("dna_small.phy_phyml_boot_trees_r3.txt")), 10)


if __name__ == "__main__":
    unittest.main()
