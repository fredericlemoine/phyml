"""Memory-safety tests, run against the AddressSanitizer/UBSan build.

build/phyml-asan is compiled with -fsanitize-recover, and these tests run it
with halt_on_error=0, so PhyML keeps going after an error and every report
of a run is collected. Each test then checks that no report points into
PhyML's code.
"""

import os
import re
import unittest

from phymltest import PHYML_ASAN, PhymlTestCase, data

# detect_stack_use_after_return (on by default with GCC >= 13) makes PhyML
# about 40 times slower under ASan; leaks are out of scope here.
ASAN_ENV = dict(os.environ,
                ASAN_OPTIONS="detect_leaks=0:halt_on_error=0:abort_on_error=0:"
                             "detect_stack_use_after_return=0",
                UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=0")


def reports(output):
    """Sanitizer reports as (kind, first PhyML frame) pairs."""
    found = []
    blocks = re.split(r"(?m)^(?===\d+==ERROR: AddressSanitizer|.*runtime error)", output)
    for block in blocks:
        m = re.match(r"==\d+==ERROR: AddressSanitizer: ([\w-]+)", block)
        if m:
            kind = m.group(1)
        else:
            m = re.match(r".*?([\w./-]+\.c:\d+:\d+): runtime error: (.*)", block)
            if not m:
                continue
            found.append(("ubsan: " + m.group(2).strip(), m.group(1)))
            continue
        frame = None
        for fm in re.finditer(r"#\d+ 0x[0-9a-f]+ in (\w+) ([\w./-]*src/)?(\w+\.c:\d+)", block):
            frame = "%s %s" % (fm.group(1), fm.group(3))
            break
        found.append((kind, frame or "?"))
    return found


@unittest.skipIf(PHYML_ASAN is None, "sanitizer build (phyml-asan) not available")
class MemorySafety(PhymlTestCase):

    binary = PHYML_ASAN

    def run_asan(self, *args, timeout=300, stdin=None):
        return self.run_phyml(*args, timeout=timeout, env=ASAN_ENV, stdin=stdin)

    def other_reports(self, res):
        return reports(res.output)

    def assertNoMemoryErrors(self, res):
        self.assertFalse(res.timed_out, res.describe())
        bad = self.other_reports(res)
        self.assertEqual(bad, [], "sanitizer reports: %r\n%s" % (bad, res.describe()))

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        return name

    # --- regression guard: everyday analyses ---------------------------
    def test_everyday_analyses_are_clean(self):
        self.stage("sim_jc_8.phy", "dna_small.phy", "aa_small.phy", "oracle_dna.phy",
                   "oracle6.nwk", "dna_small.nex", "oracle_aa.phy", "oracle5.nwk")
        runs = [
            ["-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "0"],
            ["-i", "dna_small.phy"],
            ["-i", "aa_small.phy", "-d", "aa"],
            ["-i", "dna_small.nex", "-b", "-4"],
            ["-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "3"],
            ["-i", "sim_jc_8.phy", "-m", "JC69", "-c", "1", "-b", "3", "--tbe"],
            ["-i", "oracle_dna.phy", "-q", "-u", "oracle6.nwk", "-m", "GTR", "-o", "lr",
             "--ancestral", "--print_site_lnl"],
            ["-i", "dna_small.phy", "-m", "HKY85", "-v", "e", "--freerates"],
            # Amino-acid ancestral states: Integer_To_Bit shifted a signed int
            ["-i", "oracle_aa.phy", "-d", "aa", "-u", "oracle5.nwk", "-o", "n", "--ancestral"],
            ["-i", "oracle_aa.phy", "-d", "aa", "-m", "LG4M", "-u", "oracle5.nwk", "-o", "lr",
             "--ancestral"],
            ["-i", "oracle_aa.phy", "-d", "aa", "-m", "LG4X", "-b", "2"],
        ]
        for args in runs:
            with self.subTest(args=" ".join(args)):
                res = self.run_asan(*args)
                self.assertNoMemoryErrors(res)
                self.assertEqual(res.returncode, 0, res.describe())

    def test_fasta_input_is_clean(self):
        """FASTA reader: wrapped lines, descriptions, CRLF, growing buffers."""
        seqs = {"t%d" % i: ("ACGT" * 700)[i:i + 2500] for i in range(5)}
        text = "".join(">%s description\r\n%s\r\n%s\r\n" % (n, s[:1300], s[1300:])
                       for n, s in seqs.items())
        aln = self.write("a.fasta", text)
        res = self.run_asan("-i", aln, "-m", "JC69", "-c", "1", "-b", "0")
        self.assertNoMemoryErrors(res)
        self.assertEqual(res.returncode, 0, res.describe())

    # --- known defects ---------------------------------------------------
    def test_taxon_name_read_is_bounded(self):
        """The names were read with fscanf("%1000s") into 1000-byte buffers."""
        self.stage("sim_jc_8.phy")
        for opts in ((), ("-q",)):
            with self.subTest(opts=opts):
                self.stage("sim_jc_8.phy")
                res = self.run_asan("-i", "sim_jc_8.phy", *opts, "-m", "JC69", "-c", "1",
                                    "-b", "0")
                self.assertNoMemoryErrors(res)

    def test_thousand_character_taxon_name(self):
        name = "x" * 1000
        aln = self.write("long_name.phy", "4 10\n%s ACGTACGTAC\n" % name +
                         "".join("t%d ACGTACGTAC\n" % i for i in range(3)))
        res = self.run_asan("-i", aln, "-q", "-m", "JC69")
        self.assertNoMemoryErrors(res)
        self.assertIn("longer than 999 characters", res.output)

    def test_sequence_longer_than_header(self):
        names, seqs = [], []
        with open(data("oracle_dna.phy")) as fh:
            for line in fh.read().splitlines()[1:]:
                n, s = line.split()
                names.append(n)
                seqs.append(s + "ACGTACGTAC")
        aln = self.write("long_line.phy", "6 70\n" +
                         "".join("%s %s\n" % (n, s) for n, s in zip(names, seqs)))
        self.assertNoMemoryErrors(self.run_asan("-i", aln, "-q", "-m", "JC69"))

    def test_long_input_path(self):
        """Output names add a suffix and the run id to the input name; with
        200-byte buffers, a 180-character input path overflowed them."""
        rel = os.path.join("d" * 120, "a" * 55 + ".phy")
        self.assertEqual(len(rel), 180)
        os.makedirs(os.path.join(self.tmp, "d" * 120))
        with open(data("sim_jc_8.phy")) as src, open(os.path.join(self.tmp, rel), "w") as dst:
            dst.write(src.read())
        res = self.run_asan("-i", rel, "-m", "JC69", "-c", "1", "-b", "2", "--run_id", "r" * 99)
        self.assertNoMemoryErrors(res)
        self.assertEqual(res.returncode, 0, res.describe())
        self.assertTrue(os.path.exists(os.path.join(
            self.tmp, rel + "_phyml_boot_trees_" + "r" * 99 + ".txt")))

    def test_too_long_names_are_reported_cleanly(self):
        """The 'too long' message used to be built by strcat into a 200-byte
        buffer. Names up to 3840 characters are now accepted."""
        for opts, msg in ((["-i", "x" * 4000], "too long"),
                          (["-i", "missing.phy"], "does not exist"),
                          (["-i", "sim_jc_8.phy", "--run_id", "r" * 100], "run ID is longer")):
            with self.subTest(msg=msg):
                self.stage("sim_jc_8.phy")
                res = self.run_asan(*opts)
                self.assertIn(msg, res.output)
                self.assertNotEqual(res.returncode, 0)
                self.assertNoMemoryErrors(res)

    def test_long_line_in_the_menu(self):
        """The PHYLIP-like menu read up to 2 MB into 50- to 1000-byte buffers."""
        # No argument at all (memcheck=True: no --no_memory_check), so that
        # PhyML opens the menu
        res = self.run_phyml(stdin="y" * 5000 + "\n", timeout=60, env=ASAN_ENV,
                             memcheck=True)
        self.assertNoMemoryErrors(res)
        self.assertIn("longer than", res.output)

    def test_user_tree_without_branch_lengths(self):
        self.stage("sim_jc_8.phy")
        tree = self.write("t.nwk", "((t1,t2),(t3,t4),((t5,t6),(t7,t8)));\n")
        self.assertNoMemoryErrors(self.run_asan("-i", "sim_jc_8.phy", "-u", tree, "-m", "JC69",
                                                "-o", "lr", "-b", "0", timeout=120))

    def test_user_tree_with_unclosed_comment(self):
        self.stage("sim_jc_8.phy")
        tree = self.write("t.nwk", "((t1:1,t2:1)[x:1,(t3:1,t4:1):1,((t5:1,t6:1):1,(t7:1,t8:1):1):1);\n")
        res = self.run_asan("-i", "sim_jc_8.phy", "-u", tree, "-m", "JC69", "-b", "0", timeout=120)
        self.assertNoMemoryErrors(res)
        self.assertNotEqual(res.returncode, 0, res.describe())
        self.assertIn("not closed", res.output)

    def test_user_tree_with_comments(self):
        self.stage("sim_jc_8.phy")
        tree = self.write("t.nwk", "((t1:1[x],t2:1)[note]:1,(t3:1,t4:1)[a,b]:1,"
                                   "((t5:1,t6:1):[&k=v,w]1,(t7:1,t8:1)[]:1):1);\n")
        res = self.run_asan("-i", "sim_jc_8.phy", "-u", tree, "-m", "JC69", "-o", "n",
                            "-b", "0", timeout=120)
        self.assertNoMemoryErrors(res)
        self.assertEqual(res.returncode, 0, res.describe())


class ReportParser(unittest.TestCase):
    """The parser above must recognise the sanitizer output formats."""

    SAMPLE = """
==1==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x1 at pc 0x2
WRITE of size 1 at 0x1 thread T0
    #0 0x0000 in Read_One_Line_Seq /x/phyml/src/io.c:1721
    #1 0x0000 in Read_Seq_Sequential /x/phyml/src/io.c:1504
==1==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x1 at pc 0x2
    #0 0x0000 in scanf_common(void*, int, bool, char const*, char*)+0x7c8 (libclang_rt.asan)
    #1 0x0000 in fscanf+0x74 (libclang_rt.asan)
    #2 0x0000 in Read_Seq_Interleaved /x/phyml/src/io.c:1561
/x/phyml/src/cl.c:1882:21: runtime error: member access within null pointer
"""

    def test_parser(self):
        got = reports(self.SAMPLE)
        self.assertEqual(got[0], ("heap-buffer-overflow", "Read_One_Line_Seq io.c:1721"))
        self.assertEqual(got[1], ("heap-buffer-overflow", "Read_Seq_Interleaved io.c:1561"))
        self.assertEqual(got[2][1], "/x/phyml/src/cl.c:1882:21")
        self.assertEqual(len(got), 3)


if __name__ == "__main__":
    unittest.main()
