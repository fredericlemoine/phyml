"""Model selection with SMS (--sms, --sms_criterion, --sms_only).

SMS fits the candidate models on a fixed tree and writes them to
*_phyml_sms.csv. The tests check the table (criterion formulas, number of
parameters, order), the selection path, nested models, and that the
analysis then uses the selected model. The fits were also checked against
the SMS 2.1.1 scripts run with the same PhyML (same fits to 1e-5), which
needs R and is not part of this suite.
"""

import csv
import io
import math
import os
import re
import unittest

import oracle
from phymltest import PHYML_OMP, PhymlTestCase, data, stage_files

DNA = "dna_small.phy"  # 8 taxa, 300 sites
AA = "aa_small.phy"    # 8 taxa, 200 sites: 'small' alignment (< 5000 residues)

DNA_MATRICES = {"GTR": 8, "TN93": 5, "HKY85": 4, "K80": 1}  # parameters besides lengths
AA_MATRICES = {"LG", "WAG", "JTT", "MtREV", "Dayhoff", "DCMut", "RtREV", "CpREV", "VT",
               "Blosum62", "MtMam", "MtArt", "HIVw", "HIVb", "AB", "FLU", "MtZoa",
               "Q.bird", "Q.insect", "Q.mammal", "Q.pfam", "Q.plant", "Q.yeast"}
DNA_DECOS = {"", "+G", "+I", "+G+I", "+R"}
AA_DECOS = {"", "+F", "+G", "+G+F", "+G+I", "+G+I+F", "+R", "+R+F"}


def n_sites(align):
    with open(data(align)) as fh:
        return int(fh.readline().split()[1])


def n_taxa(align):
    with open(data(align)) as fh:
        return int(fh.readline().split()[0])


def read_fits(res, align, run_id=None):
    name = "%s_phyml_sms%s.csv" % (align, "_" + run_id if run_id else "")
    rows = list(csv.DictReader(io.StringIO(res.read(name))))
    for r in rows:
        for k in ("lnL", "AIC", "BIC"):
            r[k] = float(r[k])
        r["K"] = int(r["K"])
    return rows


def selected(res):
    m = re.search(r"Selected model: (\S+) \((AIC|BIC) = ([-0-9.]+), lnL = ([-0-9.]+)\)",
                  res.output)
    return m.group(1), float(m.group(3)), float(m.group(4))


def expected_k(row, datatype, n_branches, n_catg=4):
    k = n_branches
    deco = row["Decoration"]
    if datatype == "nt":
        k += DNA_MATRICES[row["Model"]]
    elif "+F" in deco:
        k += 19
    if "+G" in deco:
        k += 1
    if "+R" in deco:
        k += 2 * n_catg - 2
    if "+I" in deco:
        k += 1
    return k


class SmsTableChecks:
    """Checks shared by the DNA and protein runs."""

    def check_table(self, rows, align, datatype, criterion="BIC", n_catg=4):
        n = n_sites(align)
        n_branches = 2 * n_taxa(align) - 3
        keys = set()
        for r in rows:
            with self.subTest(model=r["Model"], deco=r["Decoration"]):
                self.assertEqual(r["K"], expected_k(r, datatype, n_branches, n_catg))
                self.assertAlmostEqual(r["AIC"], 2 * r["K"] - 2 * r["lnL"], places=4)
                self.assertAlmostEqual(r["BIC"], r["K"] * math.log(n) - 2 * r["lnL"],
                                       places=4)
                self.assertIn(r["FixedPinv"], ("yes", "no"))
                if r["FixedPinv"] == "yes":
                    self.assertIn("+I", r["Decoration"])
            keys.add((r["Model"], r["Decoration"], r["FixedPinv"]))
        self.assertEqual(len(keys), len(rows), "a fit is listed twice")
        values = [r[criterion] for r in rows]
        self.assertEqual(values, sorted(values), "fits not sorted by %s" % criterion)


class DnaSelection(PhymlTestCase, SmsTableChecks):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.wd = cls.class_workdir()
        stage_files(cls.wd, DNA)
        cls.res = cls.class_run(cls.wd, "-i", DNA, "--sms", "--r_seed", "1", "-b", "0")
        cls.rows = read_fits(cls.res, DNA) if cls.res.returncode == 0 else []

    def setUp(self):
        super().setUp()
        self.assertSuccess(self.res)

    def test_table(self):
        self.check_table(self.rows, DNA, "nt")
        for r in self.rows:
            self.assertIn(r["Model"], DNA_MATRICES)
            self.assertIn(r["Decoration"], DNA_DECOS)

    def test_heuristic_path(self):
        fits = {(r["Model"], r["Decoration"]) for r in self.rows}
        # Steps 1-2: every decoration of GTR
        for deco in DNA_DECOS:
            self.assertIn(("GTR", deco), fits)
        # Step 4: every decoration of the selected matrix
        model = self.rows[0]["Model"]
        for deco in DNA_DECOS:
            self.assertIn((model, deco), fits)
        # Step 3 goes down GTR > TN93 > HKY85 > K80 only while it improves
        if ("K80", self.rows[0]["Decoration"]) in fits:
            self.assertIn(("HKY85", self.rows[0]["Decoration"]), fits)

    def test_selected_model_is_the_best_fit(self):
        name, crit, lnl = selected(self.res)
        best = self.rows[0]
        self.assertEqual(name, best["Model"] + best["Decoration"])
        self.assertAlmostEqual(crit, best["BIC"], places=3)
        self.assertAlmostEqual(lnl, best["lnL"], places=3)

    def test_nested_models_on_the_fixed_tree(self):
        lnl = {(r["Model"], r["Decoration"]): r["lnL"] for r in self.rows}
        pairs = [(("TN93", d), ("GTR", d)) for d in DNA_DECOS]
        pairs += [(("HKY85", d), ("TN93", d)) for d in DNA_DECOS]
        pairs += [((m, ""), (m, d)) for m in DNA_MATRICES for d in ("+G", "+I", "+R")]
        pairs += [((m, "+G"), (m, "+G+I")) for m in DNA_MATRICES]
        for simple, rich in pairs:
            if simple in lnl and rich in lnl:
                with self.subTest(simple=simple, rich=rich):
                    self.assertGreaterEqual(lnl[rich], lnl[simple] - 0.05)

    def test_analysis_uses_the_selected_model(self):
        stats = self.res.stats(DNA)
        best = self.rows[0]
        self.assertEqual(stats.model, best["Model"])
        self.assertEqual(stats.alpha is not None, "+G" in best["Decoration"])
        self.assertIn("Model selection:", self.res.output)
        # The search starts from the same BioNJ tree and re-estimates the parameters
        self.assertGreaterEqual(self.res.final_lnl(), best["lnL"] - 0.01)

    def test_reproducible(self):
        res = self.run_phyml("-i", DNA, "--sms", "--sms_only", "--r_seed", "1",
                             cwd=self.wd)
        self.assertSuccess(res)
        self.assertEqual(read_fits(res, DNA), self.rows)


class ProteinSelection(PhymlTestCase, SmsTableChecks):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.wd = cls.class_workdir()
        stage_files(cls.wd, AA)
        cls.res = cls.class_run(cls.wd, "-i", AA, "-d", "aa", "--sms_only", "--r_seed", "1",
                                timeout=300)
        cls.rows = read_fits(cls.res, AA) if cls.res.returncode == 0 else []

    def setUp(self):
        super().setUp()
        self.assertSuccess(self.res)

    def test_table(self):
        self.check_table(self.rows, AA, "aa")
        for r in self.rows:
            self.assertIn(r["Model"], AA_MATRICES)
            self.assertIn(r["Decoration"], AA_DECOS)

    def test_heuristic_path(self):
        fits = {(r["Model"], r["Decoration"]) for r in self.rows}
        for deco in ("", "+G", "+G+I", "+R"):  # Steps 1-2 on LG
            self.assertIn(("LG", deco), fits)
        model = self.rows[0]["Model"]
        for deco in AA_DECOS:                   # Step 4
            self.assertIn((model, deco), fits)
        # Step 3 tries every matrix
        self.assertEqual({m for m, _ in fits}, AA_MATRICES)
        self.assertIn("residues: small alignment", self.res.output)

    def test_selected_model_is_the_best_fit(self):
        name, crit, _ = selected(self.res)
        self.assertEqual(name, self.rows[0]["Model"] + self.rows[0]["Decoration"])
        self.assertAlmostEqual(crit, self.rows[0]["BIC"], places=3)

    def test_only_writes_the_table(self):
        self.assertFalse(self.res.exists(AA + "_phyml_tree.txt"))
        self.assertFalse(self.res.exists(AA + "_phyml_stats.txt"))


class Options(PhymlTestCase, SmsTableChecks):

    def test_aic(self):
        self.stage(DNA)
        res = self.run_phyml("-i", DNA, "--sms_only", "--sms_criterion", "aic", "--r_seed", "1")
        self.assertSuccess(res)
        rows = read_fits(res, DNA)
        self.check_table(rows, DNA, "nt", criterion="AIC")
        self.assertIn("(AIC = ", res.output)

    def test_rate_categories(self):
        self.stage(DNA)
        res = self.run_phyml("-i", DNA, "--sms_only", "-c", "6", "--r_seed", "1",
                             "--run_id", "c6")
        self.assertSuccess(res)
        self.check_table(read_fits(res, DNA, "c6"), DNA, "nt", n_catg=6)

    def test_fixed_model_parameters(self):
        """With -o tl the parameters stay at the SMS estimates."""
        self.stage(DNA)
        res = self.run_phyml("-i", DNA, "--sms", "-o", "tl", "--r_seed", "1", "-b", "0")
        self.assertSuccess(res)
        best = read_fits(res, DNA)[0]
        if "+G" in best["Decoration"]:
            shown = re.search(r"Gamma distribution parameter:\s+([0-9.]+)", res.output)
            self.assertIsNotNone(shown)
            self.assertAlmostEqual(float(shown.group(1)), res.stats(DNA).alpha, places=2)

    def test_user_tree(self):
        self.stage(DNA, "dna_small.user.nwk")
        res = self.run_phyml("-i", DNA, "--sms", "-u", "dna_small.user.nwk", "--r_seed", "1",
                             "-b", "0")
        self.assertSuccess(res)
        self.assertIn("tree from the -u tree", res.output)

    def test_multiple_data_sets(self):
        self.stage("multi3.phy")
        res = self.run_phyml("-i", "multi3.phy", "-n", "3", "--sms", "--r_seed", "1", "-b", "0")
        self.assertSuccess(res)
        self.assertEqual(res.output.count("Selected model:"), 3)

    def test_model_options_are_rejected(self):
        self.stage(DNA)
        for opts in (["-m", "GTR"], ["-f", "e"], ["-v", "e"], ["-a", "0.5"],
                     ["--freerates"], ["-m", "HKY85", "-t", "4"]):
            with self.subTest(opts=opts):
                res = self.run_phyml("-i", DNA, "--sms", *opts)
                self.assertCleanFailure(res)
                self.assertIn("cannot be combined", res.output)

    def test_invalid_settings(self):
        self.stage(DNA)
        for opts, msg in ((["--sms_criterion", "xyz"], "criterion"),
                          (["--sms", "-c", "1"], "two rate categories"),
                          (["--sms", "-d", "generic"], "DNA (-d nt) or amino acids")):
            with self.subTest(opts=opts):
                res = self.run_phyml("-i", DNA, *opts)
                self.assertCleanFailure(res)
                self.assertIn(msg, res.output)

    def test_documented_in_help(self):
        res = self.run_phyml("--help")
        for opt in ("--sms", "--sms_criterion", "--sms_only", "Q.pfam", "MtZoa"):
            self.assertIn(opt, res.output)


@unittest.skipIf(PHYML_OMP is None, "OpenMP build (phyml-omp) not available")
class Threads(PhymlTestCase):

    binary = PHYML_OMP

    def test_same_fits_for_any_number_of_threads(self):
        # 10 taxa and 1000 sites of threads_dna: enough patterns for 3 threads
        names, seqs = oracle.read_phylip(data("threads_dna.phy"))
        with open(os.path.join(self.tmp, "sub.phy"), "w") as fh:
            fh.write("10 1000\n")
            for name, seq in list(zip(names, seqs))[:10]:
                fh.write("%-10s %s\n" % (name, seq[:1000]))
        fits = {}
        for t in (1, 3):
            res = self.run_phyml("-i", "sub.phy", "--sms_only", "--r_seed", "1",
                                 "--threads", t, "--run_id", "t%d" % t, timeout=300)
            self.assertSuccess(res)
            fits[t] = read_fits(res, "sub.phy", "t%d" % t)
        self.assertEqual(fits[1], fits[3])


if __name__ == "__main__":
    unittest.main()
