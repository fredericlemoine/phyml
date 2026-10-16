"""Structure models of Le & Gascuel (2010): -m EX|EHO|EX_EHO with
--struct_model mix|part|conf_lg|conf_mix and --struct_file.

Fixed-parameter likelihoods are compared with structure_ref.py (built on the
independent oracle). Estimation is checked on struct_aa.phy, simulated under
the EX_EHO classes with noisy annotations (20% random, 5% unknown).
"""

import math
import os
import random
import re
import unittest

import oracle
import structure_ref as sr
from newick import Newick
from phymltest import PHYML_OMP, PhymlTestCase, data, stage_files

MODELS = ("mix", "part", "conf_lg", "conf_mix")
PARTITIONS = {"EX": [0.3, 0.7], "EHO": [0.2, 0.5, 0.3],
              "EX_EHO": [0.1, 0.2, 0.15, 0.25, 0.2, 0.1]}
SIM = "struct_aa.phy"
SIM_TREE = "struct_aa.true.nwk"
SIM_ANNOT = "struct_aa.annot"


def random_annotations(n, rng):
    return ("".join(rng.choice("BE?") for _ in range(n)),
            "".join(rng.choice("EHO?") for _ in range(n)))


def stats_value(res, align, label):
    m = re.search(re.escape(label) + r"\s*([-0-9.]+)", res.stats(align).text)
    return float(m.group(1)) if m else None


class ExactLikelihood(PhymlTestCase):
    """Fixed tree, alpha, chi and proportions: PhyML must match the reference."""

    def check_all(self, aln, tree_file, names, seqs, expo, sstr):
        with open(os.path.join(self.tmp, "annot.txt"), "w") as fh:
            fh.write("exposure: %s\nstructure: %s\n" % (expo, sstr))
        tree = Newick.read(os.path.join(self.tmp, tree_file))
        lg = sr.lg_matrix()
        for part, props in PARTITIONS.items():
            classes = sr.site_classes(part, expo, sstr)
            for model in MODELS:
                with self.subTest(partition=part, model=model):
                    expected = sr.log_likelihood(tree, names, seqs, part, model, classes,
                                                 0.7, 4, chi=0.6, props=props, lg=lg)
                    args = ["-i", aln, "-d", "aa", "-m", part, "--struct_model", model,
                            "-u", tree_file, "-o", "n", "-a", "0.7", "-c", "4", "-b", "0",
                            "--struct_chi", "0.6",
                            "--struct_props", ",".join(map(str, props)),
                            "--run_id", part + model]
                    if model != "mix":
                        args += ["--struct_file", "annot.txt"]
                    res = self.run_phyml(*args)
                    self.assertSuccess(res)
                    self.assertLnlClose(res.final_lnl(), expected, tol=1e-8)

    def test_all_models_and_partitions(self):
        self.stage("oracle_aa.phy", "oracle5.nwk")
        names, seqs = oracle.read_phylip(data("oracle_aa.phy"))
        expo, sstr = random_annotations(len(seqs[0]), random.Random(1))
        self.check_all("oracle_aa.phy", "oracle5.nwk", names, seqs, expo, sstr)

    def test_default_starting_values(self):
        """With -o n, chi = 0.8 and the published class proportions (the
        authors' PhyML-Structure binary gives the same likelihoods)."""
        self.stage("oracle_aa.phy", "oracle5.nwk")
        names, seqs = oracle.read_phylip(data("oracle_aa.phy"))
        expo, sstr = random_annotations(len(seqs[0]), random.Random(3))
        with open(os.path.join(self.tmp, "annot.txt"), "w") as fh:
            fh.write("exposure: %s\nstructure: %s\n" % (expo, sstr))
        tree = Newick.read(os.path.join(self.tmp, "oracle5.nwk"))
        for part, model in (("EX", "mix"), ("EHO", "part"), ("EX_EHO", "conf_mix")):
            with self.subTest(partition=part, model=model):
                classes = sr.site_classes(part, expo, sstr)
                expected = sr.log_likelihood(tree, names, seqs, part, model, classes, 0.7, 4)
                res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-m", part,
                                     "--struct_model", model, "--struct_file", "annot.txt",
                                     "-u", "oracle5.nwk", "-o", "n", "-a", "0.7", "-c", "4",
                                     "-b", "0", "--run_id", part + model)
                self.assertSuccess(res)
                self.assertLnlClose(res.final_lnl(), expected, tol=1e-8)

    def test_identical_columns_with_different_classes(self):
        """Such columns must not be merged into one site pattern."""
        self.stage("oracle5.nwk")
        names, seqs = oracle.read_phylip(data("oracle_aa.phy"))
        seqs = [s + s[:15] for s in seqs]  # 15 duplicated columns
        with open(os.path.join(self.tmp, "dup.phy"), "w") as fh:
            fh.write("%d %d\n" % (len(names), len(seqs[0])))
            fh.writelines("%-10s %s\n" % (n, s) for n, s in zip(names, seqs))
        rng = random.Random(2)
        expo = "".join(rng.choice("BE") for _ in range(len(seqs[0])))
        sstr = "".join(rng.choice("EHO") for _ in range(len(seqs[0])))
        self.check_all("dup.phy", "oracle5.nwk", names, seqs, expo, sstr)


class SimulatedData(PhymlTestCase):
    """Models fitted on the true tree of struct_aa.phy (branch lengths and
    parameters estimated)."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.wd = cls.class_workdir()
        stage_files(cls.wd, SIM, SIM_TREE, SIM_ANNOT)
        cls.res = {}
        runs = {"LG": ["-m", "LG"]}
        for m in MODELS:
            runs["EX_" + m] = ["-m", "EX", "--struct_model", m, "--struct_file", SIM_ANNOT]
        runs["EX_EHO_part"] = ["-m", "EX_EHO", "--struct_model", "part",
                               "--struct_file", SIM_ANNOT]
        for rid, opts in runs.items():
            cls.res[rid] = cls.class_run(cls.wd, "-i", SIM, "-d", "aa", "-u", SIM_TREE,
                                         "-o", "lr", "-b", "0", *opts, "--run_id", rid,
                                         timeout=600)

    def lnl(self, rid):
        res = self.res[rid]
        self.assertSuccess(res)
        return res.stats(SIM).lnl

    def test_models_are_nested(self):
        # CONF/MIX contains MIX (chi = 0) and PART (chi = 1); CONF/LG contains PART.
        tol = 0.05
        self.assertGreaterEqual(self.lnl("EX_conf_mix"), self.lnl("EX_mix") - tol)
        self.assertGreaterEqual(self.lnl("EX_conf_mix"), self.lnl("EX_part") - tol)
        self.assertGreaterEqual(self.lnl("EX_conf_lg"), self.lnl("EX_part") - tol)

    def test_structure_models_beat_lg(self):
        for rid in ("EX_part", "EX_conf_mix", "EX_EHO_part"):
            with self.subTest(model=rid):
                self.assertGreater(self.lnl(rid) - self.lnl("LG"), 50)
        # The data were simulated with the 6 EX_EHO classes
        self.assertGreater(self.lnl("EX_EHO_part"), self.lnl("EX_part"))

    def test_chi_is_estimated(self):
        """80% of the class annotations are right; a random one keeps the right
        exposure half the time, so about 90% of exposures are right."""
        chi = stats_value(self.res["EX_conf_mix"], SIM, "Confidence coefficient chi:")
        self.assertIsNotNone(chi)
        self.assertGreater(chi, 0.75)
        self.assertLess(chi, 0.99)

    def test_outputs(self):
        res = self.res["EX_conf_mix"]
        text = res.stats(SIM).text
        self.assertIn("EX_CONF/MIX", text)
        self.assertIn("exposed", text)
        rate = stats_value(res, SIM, "Branch lengths multiplied by:")
        self.assertGreater(rate, 0.8)
        self.assertLess(rate, 1.5)
        # The output tree is rescaled: its size is the one in the stats file
        self.assertAlmostEqual(res.tree(SIM).total_length(), res.stats(SIM).tree_size,
                               places=3)

    def test_fixed_proportions(self):
        """--struct_fix_props keeps the published proportions, as the authors'
        program does for MIX; estimating them cannot lower the likelihood."""
        stage_files(self.tmp, SIM, SIM_TREE)
        res = self.run_phyml("-i", SIM, "-d", "aa", "-u", SIM_TREE, "-o", "lr", "-b", "0",
                             "-m", "EX", "--struct_model", "mix", "--struct_fix_props",
                             timeout=600)
        self.assertSuccess(res)
        text = res.stats(SIM).text
        self.assertIn("proportion 0.4480", text)
        self.assertIn("proportion 0.5520", text)
        self.assertGreaterEqual(self.lnl("EX_mix"), res.stats(SIM).lnl - 0.05)

    def test_clean_annotations_give_chi_near_one(self):
        with open(data("struct_aa.classes")) as fh:
            classes = [int(c) for c in fh.read().strip()]
        with open(os.path.join(self.tmp, "clean.annot"), "w") as fh:
            fh.write("exposure: %s\nstructure: %s\n" % (
                "".join("BE"[k // 3] for k in classes), "".join("EHO"[k % 3] for k in classes)))
        stage_files(self.tmp, SIM, SIM_TREE)
        res = self.run_phyml("-i", SIM, "-d", "aa", "-u", SIM_TREE, "-o", "lr", "-b", "0",
                             "-m", "EX", "--struct_model", "conf_lg",
                             "--struct_file", "clean.annot", timeout=600)
        self.assertSuccess(res)
        self.assertGreater(stats_value(res, SIM, "Confidence coefficient chi:"), 0.95)


class Options(PhymlTestCase):

    def setUp(self):
        super().setUp()
        self.stage("oracle_aa.phy")
        names, seqs = oracle.read_phylip(data("oracle_aa.phy"))
        self.n_sites = len(seqs[0])

    def write_annot(self, text):
        with open(os.path.join(self.tmp, "a.txt"), "w") as fh:
            fh.write(text)

    def test_refused_combinations(self):
        self.write_annot("exposure: %s\n" % ("B" * self.n_sites))
        for opts, msg in ((["-m", "EX", "--struct_model", "part"], "need site annotations"),
                          (["-m", "LG", "--struct_file", "a.txt"], "need a structure model"),
                          (["-m", "EX", "-v", "e"], "cannot be combined"),
                          (["-m", "EX", "-f", "e"], "cannot be combined"),
                          (["-m", "EX", "--freerates"], "cannot be combined"),
                          (["-m", "EX", "--struct_props", "0.5"], "needs 2 proportions"),
                          (["-m", "EX", "--struct_chi", "2"], "between 0 and 1"),
                          (["-m", "EX", "--struct_model", "nope"], "must be mix"),
                          (["-m", "LG", "--struct_fix_props"], "need a structure model"),
                          (["-m", "EX", "--struct_model", "part", "--struct_file", "a.txt",
                            "--struct_fix_props"], "only applies to the mix and conf_mix")):
            with self.subTest(opts=opts):
                res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-b", "0", *opts)
                self.assertCleanFailure(res)
                self.assertIn(msg, res.output)

    def test_bad_annotation_files(self):
        for text, msg in (("exposure: BE\n", "has 2 characters"),
                          ("exposure: %s\n" % ("Z" * self.n_sites), "invalid exposure"),
                          ("accessibility: %s\n" % ("B" * self.n_sites), "unknown line"),
                          ("exposure: %s\n" % ("B" * self.n_sites), "needs a 'structure:' line")):
            with self.subTest(text=text[:20]):
                self.write_annot(text)
                res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-b", "0",
                                     "-m", "EX_EHO", "--struct_file", "a.txt")
                self.assertCleanFailure(res)
                self.assertIn(msg, res.output)

    def test_documented_in_help(self):
        res = self.run_phyml("--help")
        for opt in ("EX_EHO", "--struct_model", "--struct_file", "--struct_chi",
                    "--struct_fix_props"):
            self.assertIn(opt, res.output)


@unittest.skipIf(PHYML_OMP is None, "OpenMP build (phyml-omp) not available")
class Threads(PhymlTestCase):

    binary = PHYML_OMP

    def test_same_results_for_any_number_of_threads(self):
        self.stage(SIM, SIM_TREE, SIM_ANNOT)
        out = {}
        for t in (1, 3):
            res = self.run_phyml("-i", SIM, "-d", "aa", "-u", SIM_TREE, "-o", "lr", "-b", "0",
                                 "-m", "EX", "--struct_file", SIM_ANNOT, "--threads", t,
                                 "--run_id", "t%d" % t, timeout=600)
            self.assertSuccess(res)
            out[t] = (res.stats(SIM).lnl, res.read(SIM + "_phyml_tree_t%d.txt" % t))
        self.assertEqual(out[1], out[3])


if __name__ == "__main__":
    unittest.main()
