"""Tree search: topology recovery, optimality, starting trees, constraints."""

import unittest

from newick import Newick, rf_distance
from phymltest import PhymlTestCase, data, known_bug, stage_files

JC8, JC8_TRUE = "sim_jc_8.phy", "sim_jc_8.true.nwk"
HKY12, HKY12_TRUE = "sim_hky_12.phy", "sim_hky_12.true.nwk"


class TopologyRecovery(PhymlTestCase):
    """Long simulated alignments with clear signal: the search must find the
    generating tree and a likelihood at least as good as the true tree's."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        wd = cls.class_workdir()
        stage_files(wd, JC8, JC8_TRUE, HKY12, HKY12_TRUE)
        cls.jc_search = cls.class_run(wd, "-i", JC8, "-m", "JC69", "-c", "1", "-b", "0",
                                      "--r_seed", "1", "--run_id", "search")
        cls.jc_true = cls.class_run(wd, "-i", JC8, "-m", "JC69", "-c", "1", "-b", "0",
                                    "-u", JC8_TRUE, "-o", "lr", "--run_id", "true")
        cls.jc_bionj = cls.class_run(wd, "-i", JC8, "-m", "JC69", "-c", "1", "-b", "0",
                                     "-o", "lr", "--run_id", "bionj")
        cls.hky_search = cls.class_run(wd, "-i", HKY12, "-b", "0", "--r_seed", "1",
                                       "--run_id", "search")
        cls.hky_true = cls.class_run(wd, "-i", HKY12, "-b", "0", "-u", HKY12_TRUE,
                                     "-o", "lr", "--run_id", "true")

    def test_jc69_true_topology_found(self):
        self.assertSuccess(self.jc_search)
        self.assertEqual(rf_distance(self.jc_search.tree(JC8), Newick.read(data(JC8_TRUE))), 0)

    def test_hky_gamma_true_topology_found(self):
        self.assertSuccess(self.hky_search)
        self.assertEqual(rf_distance(self.hky_search.tree(HKY12),
                                     Newick.read(data(HKY12_TRUE))), 0)

    def test_search_at_least_as_good_as_true_tree(self):
        for search, true in ((self.jc_search, self.jc_true), (self.hky_search, self.hky_true)):
            self.assertSuccess(search)
            self.assertSuccess(true)
            self.assertGreaterEqual(search.final_lnl(), true.final_lnl() - 1e-2)

    def test_search_at_least_as_good_as_bionj_start(self):
        self.assertSuccess(self.jc_bionj)
        self.assertGreaterEqual(self.jc_search.final_lnl(), self.jc_bionj.final_lnl() - 1e-3)

    def test_output_tree_is_well_formed(self):
        for res, aln, n in ((self.jc_search, JC8, 8), (self.hky_search, HKY12, 12)):
            tree = res.tree(aln)
            names = tree.leaf_names()
            self.assertEqual(len(names), n)
            self.assertEqual(len(set(names)), n)
            self.assertEqual(len(tree.root.children), 3)  # unrooted output
            self.assertEqual(len(tree.splits()), n - 3)   # fully resolved
            for node in tree.nodes():
                if node is not tree.root:
                    self.assertGreaterEqual(node.length, 0.0)
                    self.assertEqual(len(node.children) in (0, 2), True)


class StartingTrees(PhymlTestCase):

    def setUp(self):
        super().setUp()
        self.stage(JC8, JC8_TRUE)

    def test_parsimony_starting_tree(self):
        res = self.run_phyml("-i", JC8, "-m", "JC69", "-c", "1", "-b", "0", "-p")
        self.assertSuccess(res)
        self.assertEqual(rf_distance(res.tree(JC8), Newick.read(data(JC8_TRUE))), 0)

    def test_user_starting_tree_is_improved(self):
        wrong = ("((t1:0.1,t5:0.1):0.1,(t3:0.1,t8:0.1):0.1,((t2:0.1,t6:0.1):0.1,"
                 "(t7:0.1,t4:0.1):0.1):0.1);")
        with open(self.tmp + "/wrong.nwk", "w") as fh:
            fh.write(wrong + "\n")
        res = self.run_phyml("-i", JC8, "-m", "JC69", "-c", "1", "-b", "0",
                             "-u", "wrong.nwk", "-o", "tlr")
        self.assertSuccess(res)
        self.assertEqual(rf_distance(res.tree(JC8), Newick.read(data(JC8_TRUE))), 0)

    def test_random_starts(self):
        res = self.run_phyml("-i", JC8, "-m", "JC69", "-c", "1", "-b", "0", "-s", "SPR",
                             "--rand_start", "--n_rand_starts", "3", "--r_seed", "5")
        self.assertSuccess(res)
        self.assertGreaterEqual(res.stdout.count("[Random start"), 3)
        trees = res.trees(JC8 + "_phyml_rand_trees.txt")
        self.assertGreaterEqual(len(trees), 3)
        best = res.tree(JC8)
        self.assertEqual(rf_distance(best, Newick.read(data(JC8_TRUE))), 0)

    def test_all_search_flags_give_the_same_answer(self):
        # NNI and BEST are deprecated aliases of the default search.
        lnls = []
        for move in ("NNI", "SPR", "BEST"):
            res = self.run_phyml("-i", JC8, "-m", "JC69", "-c", "1", "-b", "0",
                                 "-s", move, "--r_seed", "2")
            self.assertSuccess(res)
            lnls.append(round(res.final_lnl(), 3))
        self.assertEqual(len(set(lnls)), 1, lnls)


class Constraints(PhymlTestCase):
    """--constraint_file: the search must keep the constrained clades.

    The constraint tree used here, ((t1,t3),t2,t4,t5,t6,t7,t8), contradicts
    the true tree, which groups t1 with t2.
    """

    def write(self, text):
        with open(self.tmp + "/c.nwk", "w") as fh:
            fh.write(text + "\n")
        return "c.nwk"

    def constrained(self, text):
        self.stage(JC8)
        return self.run_phyml("-i", JC8, "-m", "JC69", "-c", "1", "-b", "0",
                              "--constraint_file", self.write(text), timeout=60)

    @known_bug("constraint trees have no branch lengths, so Read_Node_Name (io.c:327) overflows")
    def test_constraint_without_branch_lengths(self):
        with open(data("sim_jc_8.constraint.nwk")) as fh:
            res = self.constrained(fh.read().strip())
        self.assertSuccess(res)
        self.assertTrue(res.tree(JC8).has_split(["t1", "t3"]))

    @known_bug("a constraint on a subset of the taxa crashes (bus error)")
    def test_constraint_on_some_taxa(self):
        res = self.constrained("((t1:0.1,t3:0.1):0.1,t2:0.1,t4:0.1);")
        self.assertSuccess(res)
        self.assertTrue(res.tree(JC8).has_split(["t1", "t3"]))

    def test_constraint_with_lengths_on_all_taxa_runs(self):
        res = self.constrained("((t1:0.1,t3:0.1):0.1,t2:0.1,t4:0.1,t5:0.1,t6:0.1,t7:0.1,t8:0.1);")
        self.assertSuccess(res)

    @known_bug("the SPR search moves away from the constraint: output groups t1 with t2")
    def test_constraint_is_respected(self):
        res = self.constrained("((t1:0.1,t3:0.1):0.1,t2:0.1,t4:0.1,t5:0.1,t6:0.1,t7:0.1,t8:0.1);")
        self.assertSuccess(res)
        tree = res.tree(JC8)
        self.assertTrue(tree.has_split(["t1", "t3"]), tree.to_string())


class AminoAcidSearch(PhymlTestCase):

    def test_search_on_protein_data(self):
        self.stage("aa_small.phy")
        bionj = self.run_phyml("-i", "aa_small.phy", "-d", "aa", "-b", "0", "-o", "lr",
                               "--run_id", "b")
        search = self.run_phyml("-i", "aa_small.phy", "-d", "aa", "-b", "0")
        self.assertSuccess(bionj)
        self.assertSuccess(search)
        self.assertGreaterEqual(search.final_lnl(), bionj.final_lnl() - 1e-3)


if __name__ == "__main__":
    unittest.main()
