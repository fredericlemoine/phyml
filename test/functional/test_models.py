"""Substitution models and parameter optimisation on a fixed topology."""

import math
import os
import random
import unittest

import oracle
from newick import Newick, rf_distance
from phymltest import PhymlTestCase, data, stage_files

SIM = "sim_hky_12.phy"          # HKY85, kappa 4, freqs .3/.2/.2/.3, Gamma alpha 0.8
TRUE = "sim_hky_12.true.nwk"

DNA_MODELS = ["JC69", "K80", "F81", "F84", "HKY85", "TN93", "GTR"]
AA_MODELS = ["LG", "WAG", "JTT", "MtREV", "Dayhoff", "DCMut", "RtREV", "CpREV",
             "VT", "AB", "Blosum62", "MtMam", "MtArt", "HIVw", "HIVb"]

# Pairs (simpler, richer): the richer model nests the simpler one, so its
# maximised likelihood on the same topology cannot be lower.
NESTED = [("JC69", "K80"), ("JC69", "F81"), ("K80", "HKY85"), ("F81", "HKY85"),
          ("F81", "F84"), ("HKY85", "TN93"), ("F84", "TN93"), ("TN93", "GTR")]


class DnaModels(PhymlTestCase):
    """Every DNA model, branch lengths and parameters optimised on the true tree."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        wd = cls.class_workdir()
        stage_files(wd, SIM, TRUE)
        cls.results, cls.lnl = {}, {}
        for model in DNA_MODELS:
            freqs = "0.25,0.25,0.25,0.25" if model in ("JC69", "K80") else "m"
            res = cls.class_run(wd, "-i", SIM, "-u", TRUE, "-m", model, "-o", "lr",
                                "-b", "0", "-f", freqs, "--run_id", model)
            cls.results[model] = res
            if res.returncode == 0:
                cls.lnl[model] = res.final_lnl()

    def test_all_models_run(self):
        for model, res in self.results.items():
            with self.subTest(model=model):
                self.assertSuccess(res)
                self.assertTrue(math.isfinite(res.final_lnl()))

    def test_nested_models_do_not_lose_likelihood(self):
        for simple, rich in NESTED:
            with self.subTest(simple=simple, rich=rich):
                self.assertGreaterEqual(self.lnl[rich], self.lnl[simple] - 0.05)

    def test_true_model_beats_jc69_clearly(self):
        # Data simulated with kappa = 4 and unequal frequencies.
        self.assertGreater(self.lnl["HKY85"] - self.lnl["JC69"], 50)

    def test_extra_parameters_gain_little_over_true_model(self):
        # GTR has 4 more free parameters than HKY85; on HKY data the gain is
        # chi-square-like (99.99% quantile of chi2 with 4 df: ~23.5, halved).
        self.assertLess(self.lnl["GTR"] - self.lnl["HKY85"], 12.0)


class AminoAcidModels(PhymlTestCase):

    def test_all_empirical_models_run(self):
        self.stage("aa_small.phy")
        values = {}
        for model in AA_MODELS:
            with self.subTest(model=model):
                res = self.run_phyml("-i", "aa_small.phy", "-d", "aa", "-m", model,
                                     "-o", "lr", "-b", "0", "-c", "1")
                self.assertSuccess(res)
                self.assertEqual(res.stats("aa_small.phy").model, model)
                values[model] = res.final_lnl()
        # Distinct matrices give distinct likelihoods.
        self.assertGreater(len(set(round(v, 4) for v in values.values())), len(AA_MODELS) - 2)

    def test_model_frequencies_versus_empirical(self):
        self.stage("aa_small.phy")
        a = self.run_phyml("-i", "aa_small.phy", "-d", "aa", "-m", "LG", "-o", "lr",
                           "-b", "0", "-f", "m")
        b = self.run_phyml("-i", "aa_small.phy", "-d", "aa", "-m", "LG", "-o", "lr",
                           "-b", "0", "-f", "e")
        self.assertSuccess(a)
        self.assertSuccess(b)
        self.assertNotAlmostEqual(a.final_lnl(), b.final_lnl(), places=2)


class ParameterEstimation(PhymlTestCase):
    """Estimates on 2,000 simulated sites should land near the true values."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        wd = cls.class_workdir()
        stage_files(wd, SIM, TRUE)
        cls.res = cls.class_run(wd, "-i", SIM, "-u", TRUE, "-m", "HKY85", "-t", "e",
                                "-a", "e", "-f", "m", "-c", "4", "-o", "lr", "-b", "0")

    def setUp(self):
        super().setUp()
        self.assertSuccess(self.res)
        self.stats = self.res.stats(SIM)
        self.tree = self.res.tree(SIM)

    def test_kappa(self):
        self.assertTrue(3.0 < self.stats.ts_tv < 5.5, self.stats.ts_tv)

    def test_gamma_shape(self):
        self.assertTrue(0.5 < self.stats.alpha < 1.3, self.stats.alpha)

    def test_base_frequencies(self):
        for got, true in zip(self.stats.nt_freqs(), [0.3, 0.2, 0.2, 0.3]):
            self.assertAlmostEqual(got, true, delta=0.03)

    def test_tree_length(self):
        true_len = Newick.read(data(TRUE)).total_length()
        self.assertAlmostEqual(self.tree.total_length(), true_len, delta=0.15 * true_len)

    def test_topology_untouched(self):
        self.assertEqual(rf_distance(self.tree, Newick.read(data(TRUE))), 0)


class OptimisationLevels(PhymlTestCase):
    """Optimising more never lowers the likelihood; fixed things stay fixed."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        wd = cls.class_workdir()
        stage_files(wd, SIM, TRUE)
        cls.runs = {}
        for level in ("n", "l", "r", "lr", "tl", "tlr"):
            cls.runs[level] = cls.class_run(wd, "-i", SIM, "-u", TRUE, "-m", "HKY85",
                                            "-t", "e", "-a", "e", "-o", level, "-b", "0",
                                            "--run_id", level)

    def tree(self, level):
        self.assertSuccess(self.runs[level])
        return Newick.parse(self.runs[level].read(
            "%s_phyml_tree_%s.txt" % (SIM, level)).strip())

    def stats(self, level):
        self.assertSuccess(self.runs[level])
        from phymltest import Stats
        return Stats(self.runs[level].read("%s_phyml_stats_%s.txt" % (SIM, level)))

    def lnl(self, level):
        self.assertSuccess(self.runs[level])
        return self.runs[level].final_lnl()

    def test_all_levels_run(self):
        for level, res in self.runs.items():
            with self.subTest(level=level):
                self.assertSuccess(res)

    def test_more_optimisation_is_never_worse(self):
        tol = 1e-3
        self.assertGreaterEqual(self.lnl("l"), self.lnl("n") - tol)
        self.assertGreaterEqual(self.lnl("r"), self.lnl("n") - tol)
        self.assertGreaterEqual(self.lnl("lr"), self.lnl("l") - tol)
        self.assertGreaterEqual(self.lnl("lr"), self.lnl("r") - tol)
        self.assertGreaterEqual(self.lnl("tlr"), self.lnl("lr") - tol)
        self.assertGreaterEqual(self.lnl("tl"), self.lnl("l") - tol)

    def test_branch_lengths_fixed_with_r_and_n(self):
        start = Newick.read(data(TRUE)).total_length()
        for level in ("n", "r"):
            with self.subTest(level=level):
                self.assertAlmostEqual(self.tree(level).total_length(), start, places=5)

    def test_branch_lengths_move_with_l(self):
        start = Newick.read(data(TRUE)).total_length()
        self.assertNotAlmostEqual(self.tree("l").total_length(), start, places=3)

    def test_model_parameters_fixed_with_l(self):
        self.assertEqual(self.stats("l").alpha, self.stats("n").alpha)
        self.assertEqual(self.stats("l").ts_tv, self.stats("n").ts_tv)

    def test_model_parameters_move_with_r(self):
        self.assertNotAlmostEqual(self.stats("r").alpha, self.stats("n").alpha, places=2)
        self.assertNotAlmostEqual(self.stats("r").ts_tv, self.stats("n").ts_tv, places=2)

    def test_user_fixed_values_survive_full_optimisation(self):
        self.stage(SIM, TRUE)
        res = self.run_phyml("-i", SIM, "-u", TRUE, "-m", "HKY85", "-t", "2.5",
                             "-a", "0.7", "-v", "0.1", "-o", "tlr", "-b", "0")
        self.assertSuccess(res)
        st = res.stats(SIM)
        self.assertAlmostEqual(st.alpha, 0.7, places=3)
        self.assertAlmostEqual(st.ts_tv, 2.5, places=3)
        self.assertAlmostEqual(st.pinv, 0.1, places=3)

    def test_topology_fixed_without_t(self):
        truth = Newick.read(data(TRUE))
        for level in ("n", "l", "r", "lr"):
            with self.subTest(level=level):
                self.assertEqual(rf_distance(self.tree(level), truth), 0)


class RateHeterogeneity(PhymlTestCase):

    def run_on_true_tree(self, *opts):
        self.stage(SIM, TRUE)
        res = self.run_phyml("-i", SIM, "-u", TRUE, "-m", "HKY85", "-o", "lr", "-b", "0", *opts)
        self.assertSuccess(res)
        return res

    def test_gamma_improves_on_uniform_rates(self):
        uni = self.run_on_true_tree("-c", "1").final_lnl()
        gam = self.run_on_true_tree("-c", "4", "-a", "e").final_lnl()
        self.assertGreater(gam - uni, 10)  # data simulated with alpha = 0.8

    def test_more_gamma_categories_run(self):
        for k in ("2", "8", "12"):
            with self.subTest(categories=k):
                res = self.run_on_true_tree("-c", k, "-a", "e")
                self.assertEqual(len(res.stats(SIM).gamma_rates()), int(k))

    def test_invariant_sites_estimated_in_range(self):
        res = self.run_on_true_tree("-c", "4", "-a", "e", "-v", "e")
        pinv = res.stats(SIM).pinv
        self.assertTrue(0.0 <= pinv < 0.5, pinv)

    def test_freerate_model(self):
        gam = self.run_on_true_tree("-c", "4", "-a", "e").final_lnl()
        free = self.run_on_true_tree("-c", "4", "--freerates").final_lnl()
        # FreeRate generalises discrete Gamma with the same number of classes.
        self.assertGreaterEqual(free, gam - 0.5)


class Invariance(PhymlTestCase):

    def test_taxon_order_does_not_matter(self):
        names, seqs = oracle.read_phylip(data(SIM))
        order = list(range(len(names)))
        random.Random(3).shuffle(order)
        with open(os.path.join(self.tmp, "shuffled.phy"), "w") as fh:
            fh.write("%d %d\n" % (len(names), len(seqs[0])))
            for i in order:
                fh.write("%s %s\n" % (names[i], seqs[i]))
        self.stage(SIM, TRUE)
        common = ["-u", TRUE, "-m", "HKY85", "-o", "lr", "-b", "0", "-q"]
        a = self.run_phyml("-i", "shuffled.phy", *common)
        b = self.run_phyml("-i", SIM, *common[:-1])
        self.assertSuccess(a)
        self.assertSuccess(b)
        self.assertLnlClose(a.final_lnl(), b.final_lnl(), 1e-6)

    def test_site_order_does_not_matter(self):
        names, seqs = oracle.read_phylip(data(SIM))
        perm = list(range(len(seqs[0])))
        random.Random(5).shuffle(perm)
        with open(os.path.join(self.tmp, "cols.phy"), "w") as fh:
            fh.write("%d %d\n" % (len(names), len(seqs[0])))
            for n, s in zip(names, seqs):
                fh.write("%s %s\n" % (n, "".join(s[j] for j in perm)))
        self.stage(SIM, TRUE)
        a = self.run_phyml("-i", "cols.phy", "-q", "-u", TRUE, "-m", "HKY85", "-o", "lr", "-b", "0")
        b = self.run_phyml("-i", SIM, "-u", TRUE, "-m", "HKY85", "-o", "lr", "-b", "0")
        self.assertSuccess(a)
        self.assertSuccess(b)
        self.assertLnlClose(a.final_lnl(), b.final_lnl(), 1e-6)


if __name__ == "__main__":
    unittest.main()
