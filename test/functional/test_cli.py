"""Command-line behaviour: options, output files, reproducibility, errors."""

import os
import re
import unittest

from phymltest import ROOT, PhymlTestCase, known_bug

ALN = "sim_jc_8.phy"
FAST = ["-m", "JC69", "-c", "1"]


class Basics(PhymlTestCase):

    def test_version_is_printed(self):
        res = self.run_phyml("--version")
        self.assertRegex(res.output, r"PhyML version 3\.3\.\d{8}")

    def test_help_lists_main_options(self):
        res = self.run_phyml("--help")
        for opt in ("-i (or --input)", "-d (or --datatype)", "-m (or --model)",
                    "-b (or --bootstrap)", "-o params", "--r_seed"):
            self.assertIn(opt, re.sub(r"\x1b\[[0-9;]*m", "", res.output))

    @known_bug("--version and --help exit with status 1")
    def test_version_and_help_exit_zero(self):
        self.assertEqual(self.run_phyml("--version").returncode, 0)
        self.assertEqual(self.run_phyml("--help").returncode, 0)

    @known_bug("config.h is committed and stale: it still says 3.3.20250429")
    def test_version_matches_configure_ac(self):
        with open(os.path.join(ROOT, "configure.ac")) as fh:
            m = re.search(r'echo "(\d{8})"', fh.read())
        res = self.run_phyml("--version")
        self.assertIn("3.3." + m.group(1), res.output)

    def test_default_run_writes_tree_and_stats(self):
        self.stage(ALN)
        res = self.run_phyml("-i", ALN, *FAST)
        self.assertSuccess(res)
        self.assertTrue(res.exists(ALN + "_phyml_tree.txt"))
        self.assertTrue(res.exists(ALN + "_phyml_stats.txt"))
        stats = res.stats(ALN)
        self.assertEqual(stats.n_taxa, 8)
        self.assertEqual(stats.model, "JC69")
        self.assertAlmostEqual(stats.lnl, res.final_lnl(), places=4)
        tree = res.tree(ALN)
        self.assertEqual(sorted(tree.leaf_names()), ["t%d" % i for i in range(1, 9)])

    def test_default_model_is_hky85_with_gamma4(self):
        self.stage(ALN)
        res = self.run_phyml("-i", ALN, "-b", "0")
        self.assertSuccess(res)
        stats = res.stats(ALN)
        self.assertEqual(stats.model, "HKY85")
        self.assertEqual(len(stats.gamma_rates()), 4)

    def test_run_id_is_appended_to_output_names(self):
        self.stage(ALN)
        res = self.run_phyml("-i", ALN, *FAST, "--run_id", "ABC")
        self.assertSuccess(res)
        self.assertTrue(res.exists(ALN + "_phyml_tree_ABC.txt"), res.files())
        self.assertTrue(res.exists(ALN + "_phyml_stats_ABC.txt"), res.files())

    def test_quiet_mode_is_shorter(self):
        self.stage(ALN)
        loud = self.run_phyml("-i", ALN, *FAST)
        quiet = self.run_phyml("-i", ALN, *FAST, "--quiet")
        self.assertSuccess(quiet)
        self.assertLess(len(quiet.stdout), len(loud.stdout))

    def test_seed_is_recorded(self):
        self.stage(ALN)
        res = self.run_phyml("-i", ALN, *FAST, "--r_seed", "4242")
        self.assertSuccess(res)
        self.assertEqual(res.stats(ALN).seed, 4242)

    def test_same_seed_same_result(self):
        self.stage("dna_small.phy")
        a = self.run_phyml("-i", "dna_small.phy", "--r_seed", "11", "-b", "0")
        tree_a = a.read("dna_small.phy_phyml_tree.txt")
        b = self.run_phyml("-i", "dna_small.phy", "--r_seed", "11", "-b", "0")
        self.assertSuccess(a)
        self.assertSuccess(b)
        self.assertEqual(a.final_lnl(), b.final_lnl())
        self.assertEqual(tree_a, b.read("dna_small.phy_phyml_tree.txt"))

    def test_existing_outputs_are_overwritten(self):
        self.stage(ALN)
        with open(os.path.join(self.tmp, ALN + "_phyml_tree.txt"), "w") as fh:
            fh.write("garbage\n")
        res = self.run_phyml("-i", ALN, *FAST)
        self.assertSuccess(res)
        self.assertNotIn("garbage", res.read(ALN + "_phyml_tree.txt"))

    def test_input_in_another_directory(self):
        os.mkdir(os.path.join(self.tmp, "sub"))
        self.stage(ALN)
        os.rename(os.path.join(self.tmp, ALN), os.path.join(self.tmp, "sub", ALN))
        res = self.run_phyml("-i", os.path.join("sub", ALN), *FAST)
        self.assertSuccess(res)
        self.assertTrue(res.exists(os.path.join("sub", ALN + "_phyml_tree.txt")))


class OptionErrors(PhymlTestCase):

    def setUp(self):
        super().setUp()
        self.stage(ALN)

    def test_missing_input_file(self):
        res = self.run_phyml("-i", "does_not_exist.phy")
        self.assertCleanFailure(res)
        self.assertIn("does_not_exist.phy", res.output)

    def test_unknown_option(self):
        res = self.run_phyml("-i", ALN, "--no_such_option")
        self.assertCleanFailure(res)
        self.assertIn("no_such_option", res.output)

    def test_unknown_model(self):
        res = self.run_phyml("-i", ALN, "-m", "NOTAMODEL")
        self.assertCleanFailure(res)
        self.assertIn("model name is incorrect", res.output)

    def test_bootstrap_value_below_minus_five(self):
        res = self.run_phyml("-i", ALN, "-b", "-6")
        self.assertCleanFailure(res)

    @known_bug("cl.c:1205 the -b error message omits -5 (aBayes, the default)")
    def test_bootstrap_error_message_mentions_abayes(self):
        res = self.run_phyml("-i", ALN, "-b", "-6")
        self.assertRegex(res.output, r"-1 and -5|-5")

    def test_amino_acid_model_on_dna_is_rejected_or_switches_datatype(self):
        res = self.run_phyml("-i", ALN, "-m", "LG", "-b", "0")
        # LG implies amino acids; reading DNA as protein must not crash.
        self.assertFalse(res.crashed, res.describe())
        self.assertFalse(res.timed_out, res.describe())

    def test_search_option_values_are_accepted(self):
        for move in ("NNI", "SPR", "BEST"):
            with self.subTest(move=move):
                res = self.run_phyml("-i", ALN, *FAST, "-b", "0", "-s", move)
                self.assertSuccess(res)

    def test_nni_and_best_are_reported_deprecated(self):
        for move in ("NNI", "BEST"):
            with self.subTest(move=move):
                res = self.run_phyml("-i", ALN, *FAST, "-b", "0", "-s", move)
                self.assertIn("deprecated", res.output)

    @known_bug("help.c:245 still documents NNI as the default and fast search")
    def test_help_does_not_advertise_deprecated_search(self):
        res = self.run_phyml("--help")
        self.assertNotIn("NNI (default, fast)", re.sub(r"\x1b\[[0-9;]*m", "", res.output))

    @known_bug("cl.c:1882 dereferences io->mod->r_mat (NULL until -m is parsed)")
    def test_no_optimisation_with_default_model(self):
        res = self.run_phyml("-i", ALN, "-o", "n", "-b", "0")
        self.assertSuccess(res)

    def test_no_optimisation_with_explicit_model(self):
        res = self.run_phyml("-i", ALN, "-m", "HKY85", "-o", "n", "-b", "0")
        self.assertSuccess(res)


class MultipleDataSets(PhymlTestCase):

    def run_multi(self):
        self.stage("multi3.phy")
        res = self.run_phyml("-i", "multi3.phy", "-q", "-n", "3", "-b", "0", *FAST)
        self.assertSuccess(res)
        return res

    def test_one_stats_block_per_data_set(self):
        res = self.run_multi()
        self.assertEqual(len(res.stats("multi3.phy").all_lnl()), 3)

    def test_data_sets_are_analysed_independently(self):
        res = self.run_multi()
        lnls = res.stats("multi3.phy").all_lnl()
        self.assertEqual(len(set(lnls)), 3)

    @known_bug("main.c:384-397 writes each data set's tree twice to *_phyml_tree.txt")
    def test_one_tree_per_data_set(self):
        res = self.run_multi()
        self.assertEqual(len(res.trees("multi3.phy_phyml_tree.txt")), 3)

    @known_bug("main.c:333 random starts are switched off after the first data set")
    def test_random_starts_apply_to_every_data_set(self):
        self.stage("multi3.phy")
        res = self.run_phyml("-i", "multi3.phy", "-q", "-n", "2", "-b", "0", *FAST,
                             "--rand_start", "--n_rand_starts", "2", "-s", "SPR")
        self.assertSuccess(res)
        blocks = res.stdout.split("Data set [#2]")
        self.assertEqual(len(blocks), 2, res.describe())
        self.assertIn("[Random start", blocks[1])


if __name__ == "__main__":
    unittest.main()
