"""Branch supports: aLRT family, aBayes, bootstrap, TBE, Bayesian bootstrap."""

import unittest

from newick import Newick, transfer_support
from phymltest import PhymlTestCase, data, stage_files

ALN, TRUE = "sim_jc_8.phy", "sim_jc_8.true.nwk"
FAST = ["-m", "JC69", "-c", "1"]
N_INTERNAL = 5  # 8 taxa


def internal_labels(tree):
    return [n.label for n in tree.internal_nodes() if n is not tree.root]


class FastSupports(PhymlTestCase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        wd = cls.class_workdir()
        stage_files(wd, ALN, "dna_small.phy")
        cls.runs = {}
        for b in ("0", "-1", "-2", "-3", "-4", "-5"):
            cls.runs[b] = cls.class_run(wd, "-i", ALN, *FAST, "-b", b, "--run_id", "b" + b)
        cls.default = cls.class_run(wd, "-i", ALN, *FAST, "--run_id", "default")
        cls.weak = {b: cls.class_run(wd, "-i", "dna_small.phy", "-b", b, "--run_id", "w" + b)
                    for b in ("-1", "-2", "-4", "-5")}

    def tree(self, res, aln=ALN):
        self.assertSuccess(res)
        return res.tree(aln)

    def supports(self, b):
        tree = self.tree(self.runs[b])
        sup = tree.supports()
        self.assertEqual(len(sup), N_INTERNAL, tree.to_string())
        return list(sup.values())

    def test_no_supports_with_b0(self):
        self.assertTrue(all(l is None for l in internal_labels(self.tree(self.runs["0"]))))

    def test_alrt_statistic_is_non_negative(self):
        self.assertTrue(all(v >= 0 for v in self.supports("-1")))

    def test_probability_like_supports_in_unit_interval(self):
        for b in ("-2", "-3", "-4", "-5"):
            with self.subTest(b=b):
                self.assertTrue(all(0.0 <= v <= 1.0 for v in self.supports(b)))

    def test_strong_signal_gives_strong_support(self):
        # 1,500 sites simulated on well separated branches.
        for b in ("-2", "-4", "-5"):
            with self.subTest(b=b):
                self.assertTrue(all(v > 0.95 for v in self.supports(b)), self.supports(b))

    def test_min_chi2_sh_is_the_minimum(self):
        chi2 = dict(self.tree(self.runs["-2"]).supports())
        sh = dict(self.tree(self.runs["-4"]).supports())
        both = dict(self.tree(self.runs["-3"]).supports())
        for split, v in both.items():
            if split in chi2 and split in sh:
                self.assertLessEqual(v, max(chi2[split], sh[split]) + 1e-6)

    def test_default_is_abayes(self):
        self.assertEqual(self.tree(self.default).supports(),
                         self.tree(self.runs["-5"]).supports())

    def test_supports_on_weakly_resolved_data(self):
        for b, res in self.weak.items():
            with self.subTest(b=b):
                tree = self.tree(res, "dna_small.phy")
                values = list(tree.supports().values())
                self.assertTrue(values)
                if b != "-1":
                    self.assertTrue(all(0.0 <= v <= 1.0 for v in values))

    def test_topology_identical_across_support_types(self):
        ref = set(self.tree(self.runs["0"]).splits())
        for b in ("-1", "-2", "-4", "-5"):
            with self.subTest(b=b):
                self.assertEqual(set(self.tree(self.runs[b]).splits()), ref)

    def test_supports_on_fixed_topology(self):
        self.stage(ALN, TRUE)
        res = self.run_phyml("-i", ALN, *FAST, "-u", TRUE, "-o", "lr", "-b", "-4")
        tree = self.tree(res)
        self.assertEqual(len(tree.supports()), N_INTERNAL)


class Bootstrap(PhymlTestCase):

    N = 10

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        wd = cls.class_workdir()
        stage_files(wd, ALN, "dna_small.phy")
        n = str(cls.N)
        cls.fbp = cls.class_run(wd, "-i", ALN, *FAST, "-b", n, "--r_seed", "3", "--run_id", "fbp")
        cls.fbp2 = cls.class_run(wd, "-i", ALN, *FAST, "-b", n, "--r_seed", "3", "--run_id", "fbp2")
        cls.tbe = cls.class_run(wd, "-i", ALN, *FAST, "-b", n, "--tbe", "--r_seed", "3",
                                "--run_id", "tbe")
        cls.bayes = cls.class_run(wd, "-i", ALN, *FAST, "-b", n, "--bayesian_bootstrap",
                                  "--r_seed", "3", "--run_id", "bb")
        # weakly resolved data, where supports are not all maximal
        cls.wfbp = cls.class_run(wd, "-i", "dna_small.phy", "-b", n, "--r_seed", "9",
                                 "--run_id", "wfbp")
        cls.wtbe = cls.class_run(wd, "-i", "dna_small.phy", "-b", n, "--tbe", "--r_seed", "9",
                                 "--run_id", "wtbe")

    def sup(self, res, aln=ALN):
        self.assertSuccess(res)
        return res.tree(aln).supports()

    def test_replicate_trees_are_written(self):
        trees = self.fbp.trees(ALN + "_phyml_boot_trees_fbp.txt")
        self.assertEqual(len(trees), self.N)
        for t in trees:
            self.assertEqual(sorted(t.leaf_names()), sorted(Newick.read(data(TRUE)).leaf_names()))
        self.assertTrue(self.fbp.exists(ALN + "_phyml_boot_stats_fbp.txt"))

    def test_fbp_values_are_replicate_counts(self):
        values = list(self.sup(self.fbp).values())
        self.assertEqual(len(values), N_INTERNAL)
        for v in values:
            self.assertEqual(v, int(v))
            self.assertTrue(0 <= v <= self.N)

    def test_fbp_counts_match_replicate_trees(self):
        main = self.fbp.tree(ALN)
        reps = self.fbp.trees(ALN + "_phyml_boot_trees_fbp.txt")
        for split, v in main.supports().items():
            count = sum(1 for t in reps if split in t.splits())
            self.assertEqual(int(v), count)

    def test_strong_signal_full_bootstrap(self):
        self.assertTrue(all(v == self.N for v in self.sup(self.fbp).values()))

    def test_same_seed_same_bootstrap(self):
        self.assertEqual(self.fbp.read(ALN + "_phyml_boot_trees_fbp.txt"),
                         self.fbp2.read(ALN + "_phyml_boot_trees_fbp2.txt"))
        self.assertEqual(self.sup(self.fbp), self.sup(self.fbp2))

    def test_tbe_in_unit_interval(self):
        values = list(self.sup(self.tbe).values())
        self.assertEqual(len(values), N_INTERNAL)
        self.assertTrue(all(0.0 <= v <= 1.0 for v in values))

    def test_tbe_is_never_below_fbp(self):
        # With the same seed both runs see the same replicates; a branch found
        # in a replicate has transfer distance 0, so TBE >= FBP proportion.
        fbp = self.sup(self.wfbp, "dna_small.phy")
        tbe = self.sup(self.wtbe, "dna_small.phy")
        self.assertEqual(set(fbp), set(tbe))
        for split in fbp:
            self.assertGreaterEqual(tbe[split] + 1e-6, fbp[split] / self.N)

    def test_tbe_matches_reference_computation(self):
        main = self.wtbe.tree("dna_small.phy")
        reps = self.wtbe.trees("dna_small.phy_phyml_boot_trees_wtbe.txt")
        self.assertEqual(len(reps), self.N)
        expected = transfer_support(main, reps)
        for split, value in main.supports().items():
            self.assertAlmostEqual(value, expected[split], places=5, msg=sorted(split))

    def test_bayesian_bootstrap(self):
        values = list(self.sup(self.bayes).values())
        self.assertEqual(len(values), N_INTERNAL)
        self.assertTrue(all(0 <= v <= self.N for v in values))

    def test_bootstrap_rejected_with_multiple_data_sets(self):
        self.stage("multi3.phy")
        res = self.run_phyml("-i", "multi3.phy", "-q", "-n", "3", "-b", "5")
        self.assertNotEqual(res.returncode, 0)


if __name__ == "__main__":
    unittest.main()
