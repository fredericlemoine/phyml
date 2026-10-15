"""Optional outputs: ancestral states, site weights, traces, statistics."""

import json
import os
import re
import unittest

import oracle
from newick import Newick
from phymltest import PhymlTestCase, data, known_bug

ALN = "oracle_dna.phy"
TREE = "oracle6.nwk"
FIXED_JC = ["-q", "-u", TREE, "-m", "JC69", "-c", "1", "-o", "n", "-b", "0"]


class AncestralStates(PhymlTestCase):

    def setUp(self):
        super().setUp()
        self.stage(ALN, TREE)
        self.res = self.run_phyml("-i", ALN, *FIXED_JC, "--ancestral")
        self.assertSuccess(self.res)

    def table(self):
        rows = {}
        for line in self.res.read(ALN + "_phyml_ancestral_seq.txt").splitlines():
            parts = line.split()
            if len(parts) == 7 and parts[0].isdigit():
                site, node = int(parts[0]), int(parts[1])
                rows[(node, site)] = ([float(x) for x in parts[2:6]], parts[6])
        return rows

    def labelled_clades(self):
        """Map node label -> set of leaves below it in the ancestral tree."""
        tree = Newick.parse(self.res.read(ALN + "_phyml_ancestral_tree.txt").strip())
        out = {}

        def rec(node):
            leaves = (frozenset([node.name]) if node.is_leaf()
                      else frozenset().union(*(rec(c) for c in node.children)))
            if not node.is_leaf() and node.label is not None:
                out[int(node.label)] = leaves
            return leaves
        rec(tree.root)
        return out, frozenset(tree.leaf_names())

    def test_files_exist(self):
        self.assertTrue(self.res.exists(ALN + "_phyml_ancestral_seq.txt"))
        self.assertTrue(self.res.exists(ALN + "_phyml_ancestral_tree.txt"))

    def test_one_row_per_site_and_internal_node(self):
        rows = self.table()
        nodes = {k[0] for k in rows}
        self.assertEqual(len(nodes), 4)          # 6 taxa: 4 internal nodes
        self.assertEqual(len(rows), 4 * 80)

    def test_probabilities_sum_to_one(self):
        for key, (probs, _) in self.table().items():
            self.assertAlmostEqual(sum(probs), 1.0, places=4, msg=str(key))
            self.assertTrue(all(0.0 <= p <= 1.0 for p in probs))

    def test_marginals_match_reference(self):
        rows = self.table()
        clades, all_taxa = self.labelled_clades()
        names, seqs = oracle.read_phylip(data(ALN))
        ref_tree = Newick.read(data(TREE))
        ref_clades = {frozenset(n.name for n in Newick(c).leaves())
                      for c in ref_tree.internal_nodes()}
        q = oracle.gtr_q([1] * 6, [0.25] * 4)
        checked = 0
        for label, clade in clades.items():
            # The ancestral tree may be rooted elsewhere: use whichever side of
            # the node matches a node of the reference tree as written.
            match = clade if clade in ref_clades else None
            if match is None or len(match) == len(all_taxa):
                continue
            ref = oracle.ancestral_marginals(ref_tree, names, seqs, q, [0.25] * 4,
                                             oracle.DNA, match)
            for site in range(len(seqs[0])):
                got = rows[(label, site + 1)][0]
                for a, b in zip(got, ref[site]):
                    self.assertAlmostEqual(a, b, delta=2e-5,
                                           msg="node %d site %d" % (label, site + 1))
            checked += 1
        self.assertGreaterEqual(checked, 1)


class SiteWeights(PhymlTestCase):

    def lnl(self, aln, *extra):
        res = self.run_phyml("-i", aln, *FIXED_JC, *extra)
        self.assertSuccess(res)
        return res.final_lnl()

    def write(self, name, text):
        with open(os.path.join(self.tmp, name), "w") as fh:
            fh.write(text)
        return name

    def test_unit_weights_change_nothing(self):
        self.stage(ALN, TREE)
        w = self.write("w.txt", ",".join(["1"] * 80) + "\n")
        self.assertEqual(self.lnl(ALN, "--weights", w), self.lnl(ALN))

    def test_weight_two_equals_duplicated_column(self):
        self.stage(ALN, TREE)
        w = self.write("w.txt", ",".join(["2"] + ["1"] * 79) + "\n")
        names, seqs = oracle.read_phylip(data(ALN))
        dup = self.write("dup.phy", "%d %d\n" % (len(names), 81) +
                         "".join("%s %s\n" % (n, s[0] + s) for n, s in zip(names, seqs)))
        self.assertLnlClose(self.lnl(ALN, "--weights", w), self.lnl(dup), 1e-10)

    def test_zero_weight_drops_a_column(self):
        self.stage(ALN, TREE)
        w = self.write("w.txt", ",".join(["0"] + ["1"] * 79) + "\n")
        names, seqs = oracle.read_phylip(data(ALN))
        cut = self.write("cut.phy", "%d %d\n" % (len(names), 79) +
                         "".join("%s %s\n" % (n, s[1:]) for n, s in zip(names, seqs)))
        self.assertLnlClose(self.lnl(ALN, "--weights", w), self.lnl(cut), 1e-10)


class Traces(PhymlTestCase):

    def setUp(self):
        super().setUp()
        self.stage("sim_jc_8.phy")

    def test_json_trace_is_valid_json(self):
        res = self.run_phyml("-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0",
                             "--json_trace")
        self.assertSuccess(res)
        states = json.loads(res.read("sim_jc_8.phy_phyml_trace.json"))
        self.assertGreater(len(states), 1)
        for st in states:
            Newick.parse(st["state"]["tree"])
        self.assertAlmostEqual(float(states[-1]["state"]["likelihood"]),
                               res.final_lnl(), places=1)

    @known_bug("--print_trace creates *_phyml_trace.txt but writes nothing to it")
    def test_print_trace_lists_explored_trees(self):
        res = self.run_phyml("-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0",
                             "--print_trace")
        self.assertSuccess(res)
        self.assertGreater(len(res.read("sim_jc_8.phy_phyml_trace.txt").strip()), 0)


class StatsFile(PhymlTestCase):

    def test_reports_requested_fixed_parameters(self):
        self.stage(ALN, TREE)
        res = self.run_phyml("-i", ALN, "-q", "-u", TREE, "-o", "n", "-b", "0",
                             "-m", "HKY85", "-t", "3.0", "-c", "4", "-a", "0.5", "-v", "0.25")
        self.assertSuccess(res)
        st = res.stats(ALN)
        self.assertAlmostEqual(st.alpha, 0.5, places=3)
        self.assertAlmostEqual(st.pinv, 0.25, places=3)
        self.assertAlmostEqual(st.ts_tv, 3.0, places=3)
        self.assertEqual(len(st.gamma_rates()), 4)

    def test_tree_size_matches_tree_file(self):
        self.stage("sim_jc_8.phy")
        res = self.run_phyml("-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0")
        self.assertSuccess(res)
        self.assertAlmostEqual(res.stats("sim_jc_8.phy").tree_size,
                               res.tree("sim_jc_8.phy").total_length(), places=3)

    def test_parsimony_score_is_reported(self):
        self.stage("sim_jc_8.phy")
        res = self.run_phyml("-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0")
        self.assertSuccess(res)
        self.assertGreater(res.stats("sim_jc_8.phy").parsimony, 0)


if __name__ == "__main__":
    unittest.main()
