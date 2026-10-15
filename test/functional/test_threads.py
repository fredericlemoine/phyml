"""Multithreaded likelihood: --threads, in the OpenMP build (build/phyml-omp).

Threads share the site patterns of each likelihood computation, and the
per-site values are summed in site order afterwards, so the results must not
depend on the number of threads: these tests compare outputs exactly. The
threads_* data sets have enough patterns (128 per thread at least) for the
threads to be used.
"""

import unittest

from phymltest import PHYML, PHYML_OMP, PhymlTestCase

THREADS = (1, 2, 4)


def output(res, align, kind, t):
    if kind == "ancestral_seq":  # the run ID goes before "seq" in this one
        text = res.read("%s_phyml_ancestral_t%d_seq.txt" % (align, t))
    else:
        text = res.read("%s_phyml_%s_t%d.txt" % (align, kind, t))
    return text.replace("_t%d_" % t, "_tN_")  # names of other output files


@unittest.skipIf(PHYML_OMP is None, "OpenMP build (phyml-omp) not available")
class ThreadInvariance(PhymlTestCase):

    binary = PHYML_OMP

    def run_threads(self, align, *args, threads=THREADS):
        """Run the same analysis with each number of threads."""
        self.stage(align)
        if "--r_seed" not in args:  # the default seed comes from the clock
            args += ("--r_seed", "1")
        results = {}
        for t in threads:
            res = self.run_phyml("-i", align, *args, "--threads", t, "--run_id", "t%d" % t,
                                 timeout=300)
            self.assertSuccess(res)
            self.assertRegex(res.stats(align).text, r"Threads:\s*%d\b" % t)
            results[t] = res
        return results

    def assertSameOutputs(self, results, align, kinds=("tree",)):
        first = min(results)
        ref = results[first]
        for t, res in results.items():
            with self.subTest(threads=t):
                self.assertEqual(res.stats(align).all_lnl(), ref.stats(align).all_lnl())
                for kind in kinds:
                    self.assertEqual(output(res, align, kind, t), output(ref, align, kind, first),
                                     "%s differs between %d and %d threads" % (kind, first, t))

    def test_fixed_tree_site_likelihoods(self):
        self.stage("threads_dna.true.nwk")
        res = self.run_threads("threads_dna.phy", "-u", "threads_dna.true.nwk", "-m", "GTR",
                               "-c", "4", "-a", "0.6", "-v", "0.1", "-o", "n", "-b", "0",
                               "--print_site_lnl")
        self.assertSameOutputs(res, "threads_dna.phy", ("tree", "lk"))

    def test_dna_search_gamma_invariant(self):
        res = self.run_threads("threads_dna.phy", "-m", "GTR", "-c", "4", "-a", "e", "-v", "e",
                               "-o", "tlr", "-b", "0", "--r_seed", "11", threads=(1, 4))
        self.assertSameOutputs(res, "threads_dna.phy")

    def test_dna_spr_search_free_rates(self):
        res = self.run_threads("threads_dna.phy", "-m", "HKY85", "--freerates", "-c", "3",
                               "-s", "SPR", "-b", "0", "--r_seed", "5", threads=(1, 4))
        self.assertSameOutputs(res, "threads_dna.phy")

    def test_amino_acid_search(self):
        res = self.run_threads("threads_aa.phy", "-d", "aa", "-m", "LG", "-c", "4", "-a", "e",
                               "-o", "tlr", "-b", "0", "--r_seed", "3", threads=(1, 4))
        self.assertSameOutputs(res, "threads_aa.phy")

    def test_branch_supports(self):
        res = self.run_threads("threads_aa.phy", "-d", "aa", "-m", "WAG", "-c", "4",
                               "-b", "-5", "--r_seed", "3", threads=(1, 4))
        self.assertSameOutputs(res, "threads_aa.phy")

    def test_bootstrap(self):
        res = self.run_threads("threads_aa.phy", "-d", "aa", "-m", "LG", "-c", "2",
                               "-b", "3", "--r_seed", "3", threads=(1, 4))
        self.assertSameOutputs(res, "threads_aa.phy", ("tree", "boot_trees"))

    def test_ancestral_states(self):
        self.stage("threads_dna.true.nwk")
        res = self.run_threads("threads_dna.phy", "-u", "threads_dna.true.nwk", "-m", "HKY85",
                               "-c", "4", "-o", "lr", "-b", "0", "--ancestral", threads=(1, 4))
        self.assertSameOutputs(res, "threads_dna.phy", ("tree", "ancestral_seq"))

    @unittest.skipIf(PHYML is None, "default build not available")
    def test_same_as_build_without_openmp(self):
        """Up to rounding, as the two builds come from different compiler output.
        Parameters are fixed: optimisation would magnify rounding differences."""
        self.stage("threads_dna.phy", "threads_dna.true.nwk")
        args = ["-i", "threads_dna.phy", "-u", "threads_dna.true.nwk", "-m", "GTR", "-c", "4",
                "-a", "0.6", "-o", "n", "-b", "0"]
        serial = self.run_phyml(*args, binary=PHYML, timeout=300)
        threaded = self.run_phyml(*args, "--threads", "4", "--run_id", "t4", timeout=300)
        self.assertSuccess(serial)
        self.assertSuccess(threaded)
        self.assertLnlClose(threaded.final_lnl(), serial.final_lnl(), tol=1e-9)

    def test_small_alignment_runs_serially(self):
        """Below 128 patterns per thread the loops stay serial; results unchanged."""
        res = self.run_threads("oracle_dna.phy", "-q", "-m", "HKY85", "-c", "4", "-o", "tlr",
                               "-b", "0", threads=(1, 8))
        self.assertSameOutputs(res, "oracle_dna.phy")


class ThreadsOption(PhymlTestCase):

    binary = PHYML_OMP or PHYML

    def test_invalid_values_are_rejected(self):
        self.stage("sim_jc_8.phy")
        for value in ("0", "-2", "abc", "4x", "100000"):
            with self.subTest(value=value):
                res = self.run_phyml("-i", "sim_jc_8.phy", "-m", "JC69", "-b", "0",
                                     "--threads", value)
                self.assertCleanFailure(res)
                self.assertIn("number of threads", res.output)

    @unittest.skipIf(PHYML is None, "default build not available")
    def test_ignored_without_openmp(self):
        self.stage("sim_jc_8.phy")
        res = self.run_phyml("-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0",
                             "--threads", "4", binary=PHYML)
        self.assertSuccess(res)
        self.assertIn("compiled without OpenMP", res.output)

    def test_documented_in_help(self):
        res = self.run_phyml("--help")
        self.assertIn("--threads", res.output)


if __name__ == "__main__":
    unittest.main()
