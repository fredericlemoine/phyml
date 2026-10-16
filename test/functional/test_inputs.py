"""Reading alignments and trees: formats, character handling, bad input."""

import os
import unittest

import oracle

from phymltest import PhymlTestCase, data, known_bug

TREE = "oracle6.nwk"
ALN = "oracle_dna.phy"
FIXED = ["-m", "HKY85", "-t", "2.0", "-c", "4", "-a", "0.5", "-o", "n", "-b", "0"]


class AlignmentFormats(PhymlTestCase):
    """Equivalent encodings of the same data must give the same likelihood."""

    def lnl(self, aln, *extra, tree=TREE):
        self.stage(tree)
        res = self.run_phyml("-i", aln, "-u", tree, *FIXED, *extra)
        self.assertSuccess(res)
        return res.final_lnl()

    def write(self, name, text):
        with open(os.path.join(self.tmp, name), "w", newline="") as fh:
            fh.write(text)
        return name

    def original(self):
        with open(data(ALN)) as fh:
            return fh.read()

    def test_interleaved_sequential_and_nexus_agree(self):
        self.stage("dna_small.phy", "dna_small_seq.phy", "dna_small.nex")
        tree = "dna_small.user.nwk"
        a = self.lnl("dna_small.phy", tree=tree)
        b = self.lnl("dna_small_seq.phy", "-q", tree=tree)
        c = self.lnl("dna_small.nex", tree=tree)
        self.assertEqual(a, b)
        self.assertEqual(a, c)

    def test_windows_line_endings(self):
        self.stage(ALN)
        crlf = self.write("crlf.phy", self.original().replace("\n", "\r\n"))
        self.assertEqual(self.lnl(crlf, "-q"), self.lnl(ALN, "-q"))

    def test_lowercase_sequences(self):
        self.stage(ALN)
        lines = self.original().splitlines()
        low = [lines[0]] + ["%s %s" % (l.split()[0], l.split()[1].lower()) for l in lines[1:]]
        name = self.write("lower.phy", "\n".join(low) + "\n")
        self.assertEqual(self.lnl(name, "-q"), self.lnl(ALN, "-q"))

    def test_rna_u_is_read_as_t(self):
        self.stage(ALN)
        lines = self.original().splitlines()
        rna = [lines[0]] + ["%s %s" % (l.split()[0], l.split()[1].replace("T", "U"))
                            for l in lines[1:]]
        name = self.write("rna.phy", "\n".join(rna) + "\n")
        self.assertEqual(self.lnl(name, "-q"), self.lnl(ALN, "-q"))

    def test_long_taxon_names_are_kept(self):
        long_name = "a_rather_long_taxon_name_of_more_than_fifty_characters_o1"
        text = self.original().replace("o1 ", long_name + " ", 1)
        name = self.write("long.phy", text)
        with open(data(TREE)) as fh:
            self.write("long.nwk", fh.read().replace("o1:", long_name + ":"))
        res = self.run_phyml("-i", name, "-q", "-u", "long.nwk", *FIXED)
        self.assertSuccess(res)
        self.assertIn(long_name, res.tree(name).leaf_names())

    def test_amino_acid_alignment(self):
        self.stage("aa_small.phy")
        res = self.run_phyml("-i", "aa_small.phy", "-d", "aa", "-b", "0")
        self.assertSuccess(res)
        self.assertEqual(res.stats("aa_small.phy").model, "LG")
        self.assertEqual(len(res.tree("aa_small.phy").leaf_names()), 8)


class Duplicates(PhymlTestCase):

    def test_duplicates_are_removed_then_put_back(self):
        self.stage("duplicates.phy")
        res = self.run_phyml("-i", "duplicates.phy", "-q", "-b", "0", "-m", "JC69")
        self.assertSuccess(res)
        self.assertIn("duplicate", res.output)
        names = sorted(res.tree("duplicates.phy").leaf_names())
        self.assertEqual(names, sorted(["t1", "t2", "t3", "t4", "t5", "t6", "dup1", "dup2"]))

    def test_leave_duplicates(self):
        self.stage("duplicates.phy")
        res = self.run_phyml("-i", "duplicates.phy", "-q", "-b", "0", "-m", "JC69",
                             "--leave_duplicates")
        self.assertSuccess(res)
        self.assertEqual(len(res.tree("duplicates.phy").leaf_names()), 8)

    def test_duplicates_sit_together_on_zero_length_branches(self):
        self.stage("duplicates.phy")
        res = self.run_phyml("-i", "duplicates.phy", "-q", "-b", "0", "-m", "JC69")
        self.assertSuccess(res)
        tree = res.tree("duplicates.phy")
        self.assertTrue(tree.has_split(["t1", "dup1", "dup2"]), tree.to_string())
        for leaf in tree.leaves():
            if leaf.name in ("t1", "dup1", "dup2"):
                self.assertLess(leaf.length, 1e-6)


class BadAlignments(PhymlTestCase):
    """Malformed input must stop PhyML with an error message, not a crash."""

    def run_bad(self, name, *extra):
        self.stage(name)
        return self.run_phyml("-i", name, "-q", *extra)

    def test_invalid_character(self):
        res = self.run_bad("bad_char.phy")
        self.assertCleanFailure(res)
        self.assertIn("Unknown character state", res.output)

    def test_non_numeric_header(self):
        self.assertCleanFailure(self.run_bad("bad_header.phy"))

    def test_empty_file(self):
        res = self.run_bad("empty.phy")
        self.assertCleanFailure(res)
        self.assertIn("empty", res.output)

    def test_missing_input_tree_file(self):
        self.stage(ALN)
        self.assertCleanFailure(self.run_phyml("-i", ALN, "-q", "-u", "nope.nwk"))

    def test_tree_with_wrong_taxa(self):
        self.stage("sim_jc_8.phy", "bad_tree_taxa.nwk")
        res = self.run_phyml("-i", "sim_jc_8.phy", "-u", "bad_tree_taxa.nwk", "-m", "JC69")
        self.assertCleanFailure(res)

    def test_sequence_longer_than_header(self):
        """Read_One_Line_Seq used to write past the sequence buffer."""
        names, seqs = oracle.read_phylip(data("oracle_dna.phy"))
        long_first = [seqs[0] + "ACGTACGTAC"] + seqs[1:]
        for interleaved in (False, True):
            with self.subTest(interleaved=interleaved):
                path = os.path.join(self.tmp, "long.phy")
                with open(path, "w") as fh:
                    fh.write("%d %d\n" % (len(names), len(seqs[0])))
                    for n, sq in zip(names, long_first):
                        fh.write("%-10s %s\n" % (n, sq))
                opts = ("-i", "long.phy", "-m", "JC69", "-b", "0") + (() if interleaved else ("-q",))
                res = self.run_phyml(*opts)
                self.assertCleanFailure(res)
                self.assertIn("longer than", res.output)

    @known_bug("io.c:1504 a short sequence trips assert(Read_One_Line_Seq(...))")
    def test_sequence_shorter_than_header(self):
        self.assertCleanFailure(self.run_bad("bad_short_seq.phy"))

    @known_bug("io.c:1504 a truncated file trips assert(Read_One_Line_Seq(...))")
    def test_truncated_file(self):
        self.assertCleanFailure(self.run_bad("bad_truncated.phy"))

    @known_bug("io.c:1105 fewer than three taxa ends in assert(FALSE)")
    def test_fewer_than_three_taxa(self):
        res = self.run_bad("bad_too_few_taxa.phy")
        self.assertCleanFailure(res)

    @known_bug("duplicate taxon names are accepted silently")
    def test_duplicate_taxon_names(self):
        self.assertCleanFailure(self.run_bad("bad_dup_names.phy", "-b", "0"))


class UserTrees(PhymlTestCase):

    def setUp(self):
        super().setUp()
        self.stage("sim_jc_8.phy")

    def write_tree(self, text):
        with open(os.path.join(self.tmp, "user.nwk"), "w") as fh:
            fh.write(text + "\n")
        return "user.nwk"

    def run_tree(self, text, opt):
        return self.run_phyml("-i", "sim_jc_8.phy", "-u", self.write_tree(text),
                              "-m", "JC69", "-c", "1", "-o", opt, "-b", "0", timeout=60)

    def test_rooted_user_tree(self):
        res = self.run_tree("(((t1:0.1,t2:0.1):0.1,(t3:0.1,t4:0.1):0.1):0.05,"
                            "((t5:0.1,t6:0.1):0.1,(t7:0.1,t8:0.1):0.1):0.05);", "lr")
        self.assertSuccess(res)

    def test_user_tree_is_kept_when_topology_is_fixed(self):
        from newick import Newick, rf_distance
        text = "((t1:0.1,t5:0.1):0.1,(t3:0.1,t4:0.1):0.1,((t2:0.1,t6:0.1):0.1,(t7:0.1,t8:0.1):0.1):0.1);"
        res = self.run_tree(text, "lr")
        self.assertSuccess(res)
        self.assertEqual(rf_distance(res.tree("sim_jc_8.phy"), Newick.parse(text)), 0)

    def test_star_tree_with_topology_search(self):
        res = self.run_tree("(t1:0.1,t2:0.1,t3:0.1,t4:0.1,t5:0.1,t6:0.1,t7:0.1,t8:0.1);", "tlr")
        self.assertSuccess(res)

    def test_tree_without_branch_lengths(self):
        res = self.run_tree("((t1,t2),(t3,t4),((t5,t6),(t7,t8)));", "lr")
        self.assertSuccess(res)

    def test_tip_names_that_are_prefixes(self):
        """t1 must not be found inside t10 when the tree is read."""
        names, seqs = oracle.read_phylip(data("sim_jc_8.phy"))
        names = ["t10" if n == "t2" else n for n in names]
        with open(os.path.join(self.tmp, "t10.phy"), "w") as fh:
            fh.write("%d %d\n" % (len(names), len(seqs[0])))
            fh.writelines("%-10s %s\n" % (n, s) for n, s in zip(names, seqs))
        for text in ("((t10,t1),(t3,t4),((t5,t6),(t7,t8)));",
                     "((t1:0.1,t10:0.2):0.1,(t3:0.1,t4:0.1):0.1,((t5:0.1,t6:0.1):0.1,"
                     "(t7:0.1,t8:0.1):0.1):0.1);"):
            with self.subTest(tree=text):
                res = self.run_phyml("-i", "t10.phy", "-u", self.write_tree(text), "-m", "JC69",
                                     "-c", "1", "-o", "n", "-b", "0", timeout=60)
                self.assertSuccess(res)
                tree = res.tree("t10.phy")
                self.assertEqual(sorted(tree.leaf_names()), sorted(names))
                self.assertTrue(tree.has_split(["t1", "t10"]), tree.to_string())

    def test_internal_polytomy(self):
        res = self.run_tree("((t1:0.1,t2:0.1,t3:0.1):0.1,t4:0.1,"
                            "((t5:0.1,t6:0.1):0.1,(t7:0.1,t8:0.1):0.1):0.1);", "lr")
        self.assertSuccess(res)

    def test_star_tree_without_optimisation(self):
        res = self.run_tree("(t1:0.1,t2:0.1,t3:0.1,t4:0.1,t5:0.1,t6:0.1,t7:0.1,t8:0.1);", "n")
        self.assertSuccess(res)

    @known_bug("unbalanced parentheses crash the Newick reader instead of an error")
    def test_unbalanced_parentheses(self):
        self.stage("bad_tree.nwk")
        res = self.run_phyml("-i", "sim_jc_8.phy", "-u", "bad_tree.nwk", "-m", "JC69",
                             timeout=60)
        self.assertCleanFailure(res)


if __name__ == "__main__":
    unittest.main()
