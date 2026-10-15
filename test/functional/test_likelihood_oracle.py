"""PhyML's likelihood on a fixed tree with fixed parameters (-o n) must equal
the value computed by the independent reference implementation in oracle.py.
"""

import math
import unittest

import oracle
from newick import Newick
from phymltest import PHYML_GENERIC, PhymlTestCase, data, known_bug

DNA_TREE = "oracle6.nwk"
DNA_ALN = "oracle_dna.phy"
AMBIG_ALN = "oracle_dna_ambig.phy"
AA_TREE = "oracle5.nwk"
AA_ALN = "oracle_aa.phy"
FREQS = [0.1, 0.2, 0.3, 0.4]
FREQ_ARG = "0.1,0.2,0.3,0.4"
EQUAL = [0.25] * 4


def dna_reference(aln=DNA_ALN, rates6=(1,) * 6, freqs=EQUAL, alpha=None, ncat=1,
                  pinv=0.0, median=False):
    names, seqs = oracle.read_phylip(data(aln))
    tree = Newick.read(data(DNA_TREE))
    q = oracle.gtr_q(list(rates6), freqs)
    rates = [1.0] if alpha is None else oracle.discrete_gamma_rates(alpha, ncat, median)
    return oracle.log_likelihood(tree, names, seqs, q, freqs, oracle.DNA,
                                 rates=rates, pinv=pinv)


class FixedParameterLikelihood(PhymlTestCase):
    """One PhyML run per case, all on oracle_dna.phy / oracle6.nwk."""

    def phyml_lnl(self, *opts, aln=DNA_ALN, tree=DNA_TREE, binary=None):
        self.stage(aln, tree)
        res = self.run_phyml("-i", aln, "-q", "-u", tree, "-o", "n", "-b", "0",
                             *opts, binary=binary)
        self.assertSuccess(res)
        return res.final_lnl()

    def check(self, opts, expected, **kw):
        got = self.phyml_lnl(*opts, **kw)
        self.assertLnlClose(got, expected, 1e-8)

    # --- substitution models ------------------------------------------
    def test_jc69(self):
        self.check(["-m", "JC69", "-c", "1"], dna_reference())

    def test_k80_t_is_kappa(self):
        # -t is used as the rate ratio kappa of K80 (not the ts/tv ratio
        # kappa/2 described in --help).
        self.check(["-m", "K80", "-t", "2.0", "-c", "1"],
                   dna_reference(rates6=(1, 2.0, 1, 1, 2.0, 1)))

    def test_k80_kappa_matters(self):
        a = self.phyml_lnl("-m", "K80", "-t", "1.0", "-c", "1")
        b = self.phyml_lnl("-m", "K80", "-t", "5.0", "-c", "1")
        self.assertNotAlmostEqual(a, b, places=3)
        self.assertLnlClose(a, dna_reference(), 1e-8)  # K80 with kappa 1 is JC69

    def test_f81_user_frequencies(self):
        self.check(["-m", "F81", "-f", FREQ_ARG, "-c", "1"], dna_reference(freqs=FREQS))

    def test_hky85(self):
        self.check(["-m", "HKY85", "-t", "3.0", "-f", FREQ_ARG, "-c", "1"],
                   dna_reference(rates6=(1, 3.0, 1, 1, 3.0, 1), freqs=FREQS))

    def test_gtr_with_default_rates_is_f81(self):
        # Without optimisation all GTR exchangeabilities stay at 1.
        self.check(["-m", "GTR", "-f", FREQ_ARG, "-c", "1"], dna_reference(freqs=FREQS))

    def test_empirical_frequencies(self):
        names, seqs = oracle.read_phylip(data(DNA_ALN))  # no gaps or ambiguities
        freqs = oracle.empirical_freqs(seqs, oracle.DNA)
        self.check(["-m", "F81", "-f", "e", "-c", "1"], dna_reference(freqs=freqs))

    # --- rate heterogeneity -------------------------------------------
    def test_gamma_categories(self):
        for ncat in (2, 4, 8):
            with self.subTest(ncat=ncat):
                self.check(["-m", "JC69", "-c", str(ncat), "-a", "0.5"],
                           dna_reference(alpha=0.5, ncat=ncat))

    def test_gamma_shape_values(self):
        for alpha in (0.2, 1.0, 5.0):
            with self.subTest(alpha=alpha):
                self.check(["-m", "HKY85", "-t", "2.0", "-f", FREQ_ARG, "-c", "4",
                            "-a", str(alpha)],
                           dna_reference(rates6=(1, 2, 1, 1, 2, 1), freqs=FREQS,
                                         alpha=alpha, ncat=4))

    def test_gamma_median_rates(self):
        self.check(["-m", "JC69", "-c", "4", "-a", "0.5", "--use_median"],
                   dna_reference(alpha=0.5, ncat=4, median=True))

    def test_invariant_sites(self):
        # Branch lengths keep their meaning over all sites: variable sites
        # evolve 1/(1-pinv) times faster.
        self.check(["-m", "JC69", "-c", "1", "-v", "0.2"], dna_reference(pinv=0.2))

    def test_invariant_plus_gamma(self):
        self.check(["-m", "HKY85", "-t", "2.0", "-f", FREQ_ARG, "-c", "4",
                    "-a", "0.6", "-v", "0.3"],
                   dna_reference(rates6=(1, 2, 1, 1, 2, 1), freqs=FREQS,
                                 alpha=0.6, ncat=4, pinv=0.3))

    # --- data handling ------------------------------------------------
    def test_iupac_ambiguities_gaps_and_missing(self):
        self.check(["-m", "JC69", "-c", "1"], dna_reference(aln=AMBIG_ALN), aln=AMBIG_ALN)

    def test_ambiguities_with_unequal_frequencies(self):
        self.check(["-m", "HKY85", "-t", "4.0", "-f", FREQ_ARG, "-c", "4", "-a", "1.0"],
                   dna_reference(aln=AMBIG_ALN, rates6=(1, 4, 1, 1, 4, 1),
                                 freqs=FREQS, alpha=1.0, ncat=4),
                   aln=AMBIG_ALN)

    def test_stats_file_agrees_with_stdout(self):
        self.stage(DNA_ALN, DNA_TREE)
        res = self.run_phyml("-i", DNA_ALN, "-q", "-u", DNA_TREE, "-o", "n", "-b", "0",
                             "-m", "JC69", "-c", "1")
        self.assertSuccess(res)
        self.assertAlmostEqual(res.stats(DNA_ALN).lnl, res.final_lnl(), places=4)

    # --- amino acids --------------------------------------------------
    def test_custom_amino_acid_matrix(self):
        names, seqs = oracle.read_phylip(data(AA_ALN))
        tree = Newick.read(data(AA_TREE))
        exch, freqs = oracle.read_paml_matrix(data("custom_aa.dat"))
        expected = oracle.log_likelihood(tree, names, seqs, oracle.custom_aa_q(exch, freqs),
                                         freqs, oracle.AA)
        self.stage("custom_aa.dat")
        self.check(["-d", "aa", "-m", "custom", "--aa_rate_file", "custom_aa.dat",
                    "-f", "m", "-c", "1"], expected, aln=AA_ALN, tree=AA_TREE)

    def test_custom_amino_acid_matrix_with_gamma(self):
        names, seqs = oracle.read_phylip(data(AA_ALN))
        tree = Newick.read(data(AA_TREE))
        exch, freqs = oracle.read_paml_matrix(data("custom_aa.dat"))
        expected = oracle.log_likelihood(tree, names, seqs, oracle.custom_aa_q(exch, freqs),
                                         freqs, oracle.AA,
                                         rates=oracle.discrete_gamma_rates(0.8, 4))
        self.stage("custom_aa.dat")
        self.check(["-d", "aa", "-m", "custom", "--aa_rate_file", "custom_aa.dat",
                    "-f", "m", "-c", "4", "-a", "0.8"], expected, aln=AA_ALN, tree=AA_TREE)

    def test_built_in_matrices_from_sms(self):
        """Q.* (QMaker), MtZoa, AB and FLU, tested by SMS, against their PAML files.
        AB.mt is the AB_model file of the paper's supplementary material."""
        names, seqs = oracle.read_phylip(data(AA_ALN))
        tree = Newick.read(data(AA_TREE))
        for model in ("Q.pfam", "Q.bird", "Q.insect", "Q.mammal", "Q.plant", "Q.yeast",
                      "MtZoa", "AB", "FLU"):
            with self.subTest(model=model):
                exch, freqs = oracle.read_paml_matrix(data("aa_matrices/%s.mt" % model))
                freqs = [f / sum(freqs) for f in freqs]  # MtZoa's do not sum to 1
                expected = oracle.log_likelihood(tree, names, seqs,
                                                 oracle.custom_aa_q(exch, freqs),
                                                 freqs, oracle.AA)
                self.check(["-d", "aa", "-m", model, "-f", "m", "-c", "1"], expected,
                           aln=AA_ALN, tree=AA_TREE)

    # --- custom DNA models (currently crash, see test_known_bugs too) -----
    @known_bug("cl.c:956/968 writes to io->mod->r_mat before it is allocated")
    def test_custom_model_string(self):
        self.check(["-m", "000000", "-f", FREQ_ARG, "-c", "1"], dna_reference(freqs=FREQS))

    @known_bug("cl.c:980 writes to io->mod->r_mat before it is allocated")
    def test_gtr_user_rates(self):
        rates = (1.5, 3.0, 0.7, 1.2, 4.0, 1.0)
        self.check(["-m", ",".join(map(str, rates)), "-f", FREQ_ARG, "-c", "4", "-a", "0.7"],
                   dna_reference(rates6=rates, freqs=FREQS, alpha=0.7, ncat=4))


class SiteLikelihoods(PhymlTestCase):

    def test_per_site_values_match_reference(self):
        self.stage(DNA_ALN, DNA_TREE)
        res = self.run_phyml("-i", DNA_ALN, "-q", "-u", DNA_TREE, "-o", "n", "-b", "0",
                             "-m", "JC69", "-c", "4", "-a", "0.5", "--print_site_lnl")
        self.assertSuccess(res)
        rows = [l.split() for l in res.read(DNA_ALN + "_phyml_lk.txt").splitlines()]
        rows = [r for r in rows if r and r[0].isdigit()]
        names, seqs = oracle.read_phylip(data(DNA_ALN))
        ref = oracle.site_log_likelihoods(Newick.read(data(DNA_TREE)), names, seqs,
                                          oracle.gtr_q([1] * 6, EQUAL), EQUAL, oracle.DNA,
                                          rates=oracle.discrete_gamma_rates(0.5, 4))
        self.assertEqual(len(rows), len(seqs[0]))
        for i, row in enumerate(rows):
            self.assertEqual(int(row[0]), i + 1)
            self.assertAlmostEqual(math.log(float(row[1])), ref[i], places=4,
                                   msg="site %d" % (i + 1))
        total = sum(math.log(float(r[1])) for r in rows)
        self.assertAlmostEqual(total, res.final_lnl(), places=3)


@unittest.skipIf(PHYML_GENERIC is None, "scalar build (phyml-generic) not available")
class VectorisedVersusScalar(PhymlTestCase):
    """The SSE/AVX/NEON kernels must give the same likelihood as the scalar code."""

    CASES = [
        (DNA_ALN, DNA_TREE, ["-m", "HKY85", "-t", "2.0", "-f", FREQ_ARG, "-c", "4", "-a", "0.5"]),
        (AMBIG_ALN, DNA_TREE, ["-m", "JC69", "-c", "4", "-a", "0.3", "-v", "0.1"]),
        ("aa_small.phy", None, ["-d", "aa", "-m", "LG", "-c", "4", "-a", "0.7"]),
    ]

    def test_same_likelihood(self):
        for aln, tree, opts in self.CASES:
            with self.subTest(aln=aln):
                self.stage(aln)
                common = ["-i", aln, "-b", "0", "-o", "n"] + opts
                if tree:
                    self.stage(tree)
                    common += ["-q", "-u", tree]
                a = self.run_phyml(*common)
                b = self.run_phyml(*common, binary=PHYML_GENERIC)
                self.assertSuccess(a)
                self.assertSuccess(b)
                self.assertLnlClose(a.final_lnl(), b.final_lnl(), 1e-9)


if __name__ == "__main__":
    unittest.main()
