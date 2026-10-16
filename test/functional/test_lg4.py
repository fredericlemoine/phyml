"""LG4M and LG4X (Le, Dang & Gascuel 2012): -m LG4M, -m LG4X.

Fixed-parameter likelihoods and ancestral states are compared with
lg4_ref.py (built on the independent oracle, from the authors' matrix file).
Estimation is checked on lg4_aa.phy, simulated under LG4X.
"""

import math
import os
import re
import unittest

import lg4_ref
import oracle
from newick import Newick, rf_distance
from phymltest import PHYML_OMP, PhymlTestCase, data, stage_files

SIM = "lg4_aa.phy"
SIM_TREE = "lg4_aa.true.nwk"
TRUE_WEIGHTS = [0.3, 0.3, 0.25, 0.15]
_RAW = [0.15, 0.6, 1.4, 3.2]
TRUE_RATES = [r / sum(w * x for w, x in zip(TRUE_WEIGHTS, _RAW)) for r in _RAW]


def matrices(res, aln):
    """Rate and weight of each matrix, from the stats file."""
    return [(float(r), float(w)) for r, w in
            re.findall(r"Matrix [XM]\d:\s+rate ([0-9.]+), weight ([0-9.]+)", res.stats(aln).text)]


class ExactLikelihood(PhymlTestCase):
    """Fixed tree and rates: PhyML must match the reference."""

    def test_lg4m_and_lg4x(self):
        self.stage("oracle_aa.phy", "oracle5.nwk")
        names, seqs = oracle.read_phylip(data("oracle_aa.phy"))
        tree = Newick.read(data("oracle5.nwk"))
        gamma = oracle.discrete_gamma_rates(0.7, 4)
        for model in ("LG4M", "LG4X"):
            with self.subTest(model=model):
                # LG4X with -o n: equal weights and the Gamma rates of -a
                expected = lg4_ref.log_likelihood(tree, names, seqs, model, gamma, [0.25] * 4)
                res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-m", model,
                                     "-u", "oracle5.nwk", "-o", "n", "-a", "0.7", "-b", "0",
                                     "--run_id", model)
                self.assertSuccess(res)
                self.assertLnlClose(res.final_lnl(), expected, tol=1e-8)

    def test_ancestral_states(self):
        """Marginal probabilities mix the four categories, each with its own
        matrix and frequencies."""
        self.stage("oracle_aa.phy", "oracle5.nwk")
        res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-m", "LG4M", "-u", "oracle5.nwk",
                             "-o", "n", "-a", "0.7", "-b", "0", "--ancestral")
        self.assertSuccess(res)
        rows = {}
        for line in res.read("oracle_aa.phy_phyml_ancestral_seq.txt").splitlines():
            parts = line.split()
            if len(parts) == 23 and parts[0].isdigit():
                rows[(int(parts[1]), int(parts[0]))] = [float(x) for x in parts[2:22]]
        anc = Newick.parse(res.read("oracle_aa.phy_phyml_ancestral_tree.txt").strip())
        names, seqs = oracle.read_phylip(data("oracle_aa.phy"))
        ref_tree = Newick.read(data("oracle5.nwk"))
        ref_clades = {frozenset(n.name for n in Newick(c).leaves())
                      for c in ref_tree.internal_nodes()}
        mats = lg4_ref.read_matrices()
        gamma = oracle.discrete_gamma_rates(0.7, 4)
        checked = 0
        for node in anc.internal_nodes():
            clade = frozenset(n.name for n in Newick(node).leaves())
            if node.label is None or clade not in ref_clades or len(clade) == len(names):
                continue
            joint = None
            for k in range(4):
                exch, freqs = mats["M%d" % (k + 1)]
                q = [[x * gamma[k] for x in row] for row in oracle.custom_aa_q(exch, freqs)]
                part = oracle.ancestral_marginals(ref_tree, names, seqs, q, freqs, oracle.AA,
                                                  clade, normalise=False)
                joint = part if joint is None else [[a + b for a, b in zip(x, y)]
                                                    for x, y in zip(joint, part)]
            for site, p in enumerate(joint):
                tot = sum(p)
                got = rows[(int(node.label), site + 1)]
                for a, b in zip(got, p):
                    self.assertAlmostEqual(a, b / tot, delta=2e-5)
            checked += 1
        self.assertGreaterEqual(checked, 1)


class SimulatedData(PhymlTestCase):
    """lg4_aa.phy was simulated under LG4X: weights 0.3, 0.3, 0.25, 0.15."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.wd = cls.class_workdir()
        stage_files(cls.wd, SIM, SIM_TREE)
        cls.res = {m: cls.class_run(cls.wd, "-i", SIM, "-d", "aa", "-m", m, "-u", SIM_TREE,
                                    "-o", "lr", "-b", "0", "--run_id", m, timeout=600)
                   for m in ("LG", "LG4M", "LG4X")}

    def lnl(self, model):
        res = self.res[model]
        self.assertSuccess(res)
        return res.stats(SIM).lnl

    def test_lg4x_fits_best(self):
        self.assertGreater(self.lnl("LG4X") - self.lnl("LG4M"), 30)
        self.assertGreater(self.lnl("LG4M") - self.lnl("LG"), 30)

    def test_constraints_and_recovery(self):
        rw = matrices(self.res["LG4X"], SIM)
        self.assertEqual(len(rw), 4)
        self.assertAlmostEqual(sum(w for _, w in rw), 1.0, places=3)
        self.assertAlmostEqual(sum(r * w for r, w in rw), 1.0, places=3)
        for (r, w), tr, tw in zip(rw, TRUE_RATES, TRUE_WEIGHTS):
            self.assertLess(abs(r - tr), 0.2 * tr + 0.02)
            self.assertLess(abs(w - tw), 0.05)

    def test_lg4m_categories(self):
        rw = matrices(self.res["LG4M"], SIM)
        self.assertEqual([w for _, w in rw], [0.25] * 4)
        alpha = float(re.search(r"Gamma shape parameter:\s+([0-9.]+)",
                                self.res["LG4M"].stats(SIM).text).group(1))
        expected = oracle.discrete_gamma_rates(alpha, 4)
        for (r, _), e in zip(rw, expected):
            self.assertAlmostEqual(r, e, delta=2e-3)

    def test_reported_state_reproduces_likelihood(self):
        """The output tree and the printed rates and weights give back the
        printed likelihood."""
        names, seqs = oracle.read_phylip(data(SIM))
        res = self.res["LG4X"]
        rw = matrices(res, SIM)
        got = lg4_ref.log_likelihood(res.tree(SIM), names, seqs, "LG4X",
                                     [r for r, _ in rw], [w for _, w in rw])
        self.assertAlmostEqual(got, self.lnl("LG4X"), delta=0.05)

    def test_search_recovers_true_tree(self):
        stage_files(self.tmp, SIM, SIM_TREE)
        res = self.run_phyml("-i", SIM, "-d", "aa", "-m", "LG4X", "-b", "0", "--r_seed", "3",
                             timeout=600)
        self.assertSuccess(res)
        self.assertEqual(rf_distance(res.tree(SIM), Newick.read(data(SIM_TREE))), 0)
        # The search starts from the LG4M fit, then optimises LG4X
        self.assertIn("first fitting LG4M", res.output)
        self.assertGreaterEqual(res.stats(SIM).lnl, self.lnl("LG4X") - 0.05)


class Options(PhymlTestCase):

    def test_refused_combinations(self):
        self.stage("oracle_aa.phy")
        for opts in (["-f", "e"], ["-v", "e"], ["--freerates"], ["-c", "8"], ["--use_median"],
                     ["--sms"]):
            for model in ("LG4M", "LG4X"):
                with self.subTest(model=model, opts=opts):
                    res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-m", model,
                                         "-b", "0", *opts)
                    self.assertCleanFailure(res)

    def test_documented_in_help(self):
        res = self.run_phyml("--help")
        self.assertIn("LG4M", res.output)
        self.assertIn("LG4X", res.output)

    def test_stats_name_the_matrices(self):
        self.stage("oracle_aa.phy", "oracle5.nwk")
        res = self.run_phyml("-i", "oracle_aa.phy", "-d", "aa", "-m", "LG4X", "-u", "oracle5.nwk",
                             "-o", "n", "-b", "0")
        self.assertSuccess(res)
        text = res.stats("oracle_aa.phy").text
        self.assertIn("LG4X", text)
        for k in range(1, 5):
            self.assertIn("Matrix X%d" % k, text)
        self.assertNotIn("f(Ala)", text)


@unittest.skipIf(PHYML_OMP is None, "OpenMP build (phyml-omp) not available")
class Threads(PhymlTestCase):

    binary = PHYML_OMP

    def test_same_results_for_any_number_of_threads(self):
        self.stage(SIM, SIM_TREE)
        out = {}
        for t in (1, 3):
            res = self.run_phyml("-i", SIM, "-d", "aa", "-m", "LG4X", "-u", SIM_TREE, "-o", "lr",
                                 "-b", "0", "--threads", t, "--run_id", "t%d" % t, timeout=600)
            self.assertSuccess(res)
            out[t] = (res.stats(SIM).lnl, res.read(SIM + "_phyml_tree_t%d.txt" % t))
        self.assertEqual(out[1], out[3])


if __name__ == "__main__":
    unittest.main()
