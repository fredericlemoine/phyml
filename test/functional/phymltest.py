"""Shared helpers for the PhyML functional tests.

The binary under test is taken from the PHYML environment variable, falling
back to test/build/phyml (built by `make -C test`). Optional variants built by
the test Makefile are picked up the same way:

  PHYML_GENERIC  binary built with -DDISABLE_NATIVE (scalar likelihood code)
  PHYML_ASAN     binary built with -fsanitize=address,undefined
  PHYML_NDEBUG   binary built with -DNDEBUG (assertions compiled out)
  PHYML_OMP      binary built with OpenMP (--threads)
"""

import os
import re
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(HERE)
ROOT = os.path.dirname(TEST_DIR)
DATA = os.path.join(TEST_DIR, "data")
BUILD = os.path.join(TEST_DIR, "build")

DEFAULT_TIMEOUT = float(os.environ.get("PHYML_TEST_TIMEOUT", "120"))


def binary(var="PHYML", default="phyml"):
    path = os.path.abspath(os.environ.get(var) or os.path.join(BUILD, default))
    return path if os.path.isfile(path) and os.access(path, os.X_OK) else None


PHYML = binary()
PHYML_GENERIC = binary("PHYML_GENERIC", "phyml-generic")
PHYML_ASAN = binary("PHYML_ASAN", "phyml-asan")
PHYML_NDEBUG = binary("PHYML_NDEBUG", "phyml-ndebug")
PHYML_OMP = binary("PHYML_OMP", "phyml-omp")


def data(name):
    return os.path.join(DATA, name)


class Result:
    def __init__(self, cwd, args, proc, timed_out=False):
        self.cwd = cwd
        self.args = args
        self.returncode = None if timed_out else proc.returncode
        self.stdout = proc.stdout if proc else ""
        self.stderr = proc.stderr if proc else ""
        self.timed_out = timed_out

    @property
    def output(self):
        return (self.stdout or "") + (self.stderr or "")

    @property
    def crashed(self):
        """True when the process died from a signal (segfault, abort...)."""
        return self.returncode is not None and (self.returncode < 0 or self.returncode > 128)

    def path(self, name):
        return os.path.join(self.cwd, name)

    def exists(self, name):
        return os.path.isfile(self.path(name))

    def read(self, name):
        with open(self.path(name)) as fh:
            return fh.read()

    def files(self):
        return sorted(os.listdir(self.cwd))

    def describe(self):
        tail = self.output[-1500:]
        return "phyml %s\nreturn code: %s%s\n--- output tail ---\n%s" % (
            " ".join(self.args), self.returncode,
            " (timed out)" if self.timed_out else "", tail)

    # ---- parsed outputs -------------------------------------------------
    def final_lnl(self):
        """Full-precision log-likelihood printed on stdout."""
        m = re.findall(r"Log likelihood of the current tree: (-?[\d.]+(?:e[-+]?\d+)?)",
                       self.stdout)
        if not m:
            raise AssertionError("no log-likelihood in output\n" + self.describe())
        return float(m[-1].rstrip("."))

    def stats(self, align):
        return Stats(self.read(stats_name(align, self)))

    def tree(self, align):
        from newick import Newick
        return Newick.parse(self.read(tree_name(align, self)).strip().split("\n")[0])

    def trees(self, name):
        from newick import Newick
        return [Newick.parse(l) for l in self.read(name).split("\n") if l.strip()]


def _suffix(result):
    if "--run_id" in result.args:
        return "_" + result.args[result.args.index("--run_id") + 1]
    return ""


def stats_name(align, result):
    return "%s_phyml_stats%s.txt" % (os.path.basename(align), _suffix(result))


def tree_name(align, result):
    return "%s_phyml_tree%s.txt" % (os.path.basename(align), _suffix(result))


class Stats:
    """Parser for *_phyml_stats.txt."""

    def __init__(self, text):
        self.text = text

    def _num(self, label):
        m = re.search(re.escape(label) + r"[^\n]*?(-?\d+(?:\.\d+)?(?:e[-+]?\d+)?)", self.text)
        return float(m.group(1)) if m else None

    @property
    def lnl(self):
        return self._num("Log-likelihood:")

    @property
    def parsimony(self):
        return self._num("Parsimony:")

    @property
    def tree_size(self):
        return self._num("Tree size:")

    @property
    def n_taxa(self):
        v = self._num("Number of taxa:")
        return int(v) if v is not None else None

    @property
    def alpha(self):
        return self._num("Gamma shape parameter:")

    @property
    def pinv(self):
        return self._num("Proportion of invariant:")

    @property
    def ts_tv(self):
        return self._num("Transition/transversion ratio:")

    @property
    def seed(self):
        v = self._num("Random seed:")
        return int(v) if v is not None else None

    @property
    def model(self):
        m = re.search(r"Model of (?:nucleotides|amino acids) substitution:\s*(\S+)", self.text)
        return m.group(1) if m else None

    def gamma_rates(self):
        return [float(x) for x in re.findall(r"Relative rate in class \d+:\s*([\d.]+)", self.text)]

    def nt_freqs(self):
        return [float(x) for x in re.findall(r"f\([ACGT]\)=\s*([\d.]+)", self.text)]

    def all_lnl(self):
        return [float(x) for x in re.findall(r"Log-likelihood:\s*(-?[\d.]+)", self.text)]


def stage_files(dest, *names):
    """Copy data files into `dest`; return their base names."""
    for n in names:
        shutil.copy(data(n), os.path.join(dest, os.path.basename(n)))
    return [os.path.basename(n) for n in names]


def run(args, cwd, exe=None, timeout=DEFAULT_TIMEOUT, memcheck=False, env=None, stdin=None):
    """Run PhyML with stdin closed, or reading the text `stdin`; never raises
    on failure or timeout."""
    exe = exe or PHYML
    args = [str(a) for a in args]
    if not memcheck and "--no_memory_check" not in args:
        args.append("--no_memory_check")
    try:
        proc = subprocess.run([exe] + args, cwd=cwd,
                              stdin=subprocess.DEVNULL if stdin is None else None,
                              input=stdin,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, errors="replace",
                              timeout=timeout, env=env)
        return Result(cwd, args, proc)
    except subprocess.TimeoutExpired as exc:
        def text(x):
            return x.decode(errors="replace") if isinstance(x, bytes) else (x or "")

        class P:  # minimal stand-in for CompletedProcess
            returncode = None
            stdout = text(exc.stdout)
            stderr = text(exc.stderr)
        return Result(cwd, args, P, timed_out=True)


class PhymlTestCase(unittest.TestCase):
    """Base class: each test runs PhyML in a fresh temporary directory.

    Classes that share expensive runs between tests do them in setUpClass,
    in the directory returned by cls.class_workdir().
    """

    binary = None  # overridden to test another build

    @classmethod
    def setUpClass(cls):
        if (cls.binary or PHYML) is None:
            raise unittest.SkipTest("PhyML binary not found: build it with "
                                    "`make -C test` or set PHYML=/path/to/phyml")

    @classmethod
    def class_workdir(cls):
        d = tempfile.mkdtemp(prefix="phyml-test-")
        cls.addClassCleanup(shutil.rmtree, d, True)
        return d

    @classmethod
    def class_run(cls, workdir, *args, **kw):
        kw.setdefault("exe", cls.binary)
        return run(args, workdir, **kw)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="phyml-test-")
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def stage(self, *names):
        """Copy data files into the working directory; return their names."""
        return stage_files(self.tmp, *names)

    def run_phyml(self, *args, binary=None, timeout=DEFAULT_TIMEOUT, cwd=None,
                  memcheck=False, env=None, stdin=None):
        return run(args, cwd or self.tmp, binary or self.binary, timeout, memcheck, env,
                   stdin)

    def assertSuccess(self, res):
        self.assertFalse(res.timed_out, res.describe())
        self.assertEqual(res.returncode, 0, res.describe())

    def assertCleanFailure(self, res):
        """Non-zero exit without crashing, hanging or tripping an assert."""
        self.assertFalse(res.timed_out, "PhyML hung\n" + res.describe())
        self.assertFalse(res.crashed, "PhyML crashed\n" + res.describe())
        self.assertNotEqual(res.returncode, 0, res.describe())
        self.assertNotIn("Assertion failed", res.output, res.describe())

    def assertLnlClose(self, a, b, tol=1e-6, msg=None):
        self.assertLessEqual(abs(a - b), tol * max(1.0, abs(b)),
                             msg or "%r != %r (tol %g)" % (a, b, tol))


def known_bug(reason):
    """Mark a test that checks correct behaviour PhyML does not have yet.

    It is reported as an expected failure while the bug exists, and as an
    unexpected success once it is fixed, which reminds the maintainer to drop
    the marker. `reason` names the defect.

    Set PHYML_SHOW_KNOWN_BUGS=1 to run these tests as ordinary tests and see
    how each one fails.
    """
    def deco(fn):
        fn.__known_bug__ = reason
        if os.environ.get("PHYML_SHOW_KNOWN_BUGS"):
            return fn
        return unittest.expectedFailure(fn)
    return deco
