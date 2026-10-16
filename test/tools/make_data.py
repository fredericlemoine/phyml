#!/usr/bin/env python3
"""Regenerate the data files used by the PhyML test suite.

The generated files are committed under test/data/, so this script only
needs to be run when a data set has to change. It is deterministic: it
derives small alignments from examples/ and simulates others with a fixed
seed (Python's Mersenne Twister gives the same stream on every platform).

Usage: python3 test/tools/make_data.py
"""

import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "test", "data")
EXAMPLES = os.path.join(ROOT, "examples")

sys.path.insert(0, os.path.join(ROOT, "test", "functional"))
from oracle import (DNA, AA, Newick, gtr_q, custom_aa_q, expm,  # noqa: E402
                    discrete_gamma_rates)


def read_phylip_interleaved(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    ntax, nsite = map(int, lines[0].split())
    names, seqs = [], []
    i = 1
    for _ in range(ntax):
        names.append(lines[i][:10].strip())
        seqs.append(lines[i][10:].replace(" ", ""))
        i += 1
    while any(len(s) < nsite for s in seqs):
        while lines[i].strip() == "":
            i += 1
        for k in range(ntax):
            seqs[k] += lines[i].replace(" ", "")
            i += 1
    assert all(len(s) == nsite for s in seqs), "bad parse of " + path
    return names, seqs


def write_phylip(path, names, seqs, interleaved=True, block=60):
    nsite = len(seqs[0])
    with open(path, "w") as fh:
        fh.write("%d %d\n" % (len(names), nsite))
        if not interleaved:
            for n, s in zip(names, seqs):
                fh.write("%-10s %s\n" % (n, s))
            return
        for start in range(0, nsite, block):
            for k, (n, s) in enumerate(zip(names, seqs)):
                prefix = "%-10s " % n if start == 0 else " " * 11
                chunk = s[start:start + block]
                fh.write(prefix + " ".join(chunk[j:j + 10]
                                           for j in range(0, len(chunk), 10)) + "\n")
            fh.write("\n")


def write_nexus(path, names, seqs, datatype="dna"):
    with open(path, "w") as fh:
        fh.write("#NEXUS\n\nBEGIN DATA;\n")
        fh.write("  DIMENSIONS NTAX=%d NCHAR=%d;\n" % (len(names), len(seqs[0])))
        fh.write("  FORMAT DATATYPE=%s MISSING=? GAP=-;\n" % datatype.upper())
        fh.write("  MATRIX\n")
        for n, s in zip(names, seqs):
            fh.write("    %-10s %s\n" % (n, s))
        fh.write("  ;\nEND;\n")


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def simulate(tree, nsite, alphabet, q, freqs, rng, alpha=None, ncat=4):
    """Simulate sequences down an unrooted tree (rooted at its first node)."""
    rates = [1.0] if alpha is None else discrete_gamma_rates(alpha, ncat)
    pcache = {}

    def pmat(length, r):
        key = (length, r)
        if key not in pcache:
            pcache[key] = expm(q, length * r)
        return pcache[key]

    def draw(probs):
        u = rng.random()
        acc = 0.0
        for i, p in enumerate(probs):
            acc += p
            if u < acc:
                return i
        return len(probs) - 1

    out = {name: [] for name in tree.leaf_names()}
    for _ in range(nsite):
        r = rates[rng.randrange(len(rates))]
        root_state = draw(freqs)

        def down(node, state):
            for child in node.children:
                p = pmat(child.length, r)
                s = draw(p[state])
                if child.children:
                    down(child, s)
                else:
                    out[child.name].append(alphabet[s])
        down(tree.root, root_state)
    names = tree.leaf_names()
    return names, ["".join(out[n]) for n in names]


def main():
    rng = random.Random(20260528)

    # --- Subsets of the bundled examples --------------------------------
    names, seqs = read_phylip_interleaved(os.path.join(EXAMPLES, "nucleic"))
    keep = [0, 1, 3, 4, 5, 6, 7, 8]
    sn = [names[i] for i in keep]
    ss = [seqs[i][:300] for i in keep]
    write_phylip(os.path.join(DATA, "dna_small.phy"), sn, ss, interleaved=True)
    write_phylip(os.path.join(DATA, "dna_small_seq.phy"), sn, ss, interleaved=False)
    write_nexus(os.path.join(DATA, "dna_small.nex"), sn, ss)

    names, seqs = read_phylip_interleaved(os.path.join(EXAMPLES, "proteic"))
    an = names[:8]
    aseq = [s[:200] for s in seqs[:8]]
    write_phylip(os.path.join(DATA, "aa_small.phy"), an, aseq, interleaved=True)

    # Three data sets in one file for -n 3.
    with open(os.path.join(DATA, "multi3.phy"), "w") as fh:
        for k in range(3):
            sub = [s[k * 100:(k + 1) * 100] for s in ss]
            fh.write("%d %d\n" % (len(sn), 100))
            for n, s in zip(sn, sub):
                fh.write("%-10s %s\n" % (n, s))
            fh.write("\n")

    # --- Simulated alignments on known trees ----------------------------
    t8 = Newick.parse("((t1:0.10,t2:0.12):0.08,(t3:0.09,t4:0.11):0.07,"
                      "((t5:0.10,t6:0.08):0.09,(t7:0.12,t8:0.10):0.06):0.08);")
    write(os.path.join(DATA, "sim_jc_8.true.nwk"), t8.to_string() + "\n")
    q = gtr_q([1] * 6, [0.25] * 4)
    n8, s8 = simulate(t8, 1500, DNA, q, [0.25] * 4, rng)
    write_phylip(os.path.join(DATA, "sim_jc_8.phy"), n8, s8)

    t12 = Newick.parse(
        "(((a1:0.06,a2:0.08):0.05,(a3:0.07,a4:0.05):0.06):0.05,"
        "((a5:0.08,a6:0.06):0.07,(a7:0.05,a8:0.09):0.05):0.06,"
        "((a9:0.07,a10:0.06):0.06,(a11:0.08,a12:0.07):0.05):0.05);")
    write(os.path.join(DATA, "sim_hky_12.true.nwk"), t12.to_string() + "\n")
    freqs = [0.3, 0.2, 0.2, 0.3]
    kappa = 4.0
    q = gtr_q([1, kappa, 1, 1, kappa, 1], freqs)
    n12, s12 = simulate(t12, 2000, DNA, q, freqs, rng, alpha=0.8)
    write_phylip(os.path.join(DATA, "sim_hky_12.phy"), n12, s12)

    # --- Small data for the likelihood oracle ---------------------------
    t6 = Newick.parse("((o1:0.05,o2:0.10):0.04,o3:0.20,((o4:0.08,o5:0.15):0.03,"
                      "o6:0.30):0.07);")
    write(os.path.join(DATA, "oracle6.nwk"), t6.to_string() + "\n")
    freqs = [0.1, 0.2, 0.3, 0.4]
    q = gtr_q([1.5, 3.0, 0.7, 1.2, 4.0, 1.0], freqs)
    n6, s6 = simulate(t6, 80, DNA, q, freqs, rng, alpha=0.7)
    write_phylip(os.path.join(DATA, "oracle_dna.phy"), n6, s6, interleaved=False)

    # Same data with ambiguity codes, gaps and missing data sprinkled in.
    amb = [list(s) for s in s6]
    codes = "NRYKMSWBDHV-?"
    for k in range(40):
        i = rng.randrange(len(amb))
        j = rng.randrange(len(amb[0]))
        amb[i][j] = codes[k % len(codes)]
    write_phylip(os.path.join(DATA, "oracle_dna_ambig.phy"), n6,
                 ["".join(a) for a in amb], interleaved=False)

    # Amino acids under a custom (PAML format) exchangeability matrix.
    exch = [[0.0] * 20 for _ in range(20)]
    for i in range(20):
        for j in range(i):
            exch[i][j] = exch[j][i] = round(rng.uniform(0.2, 3.0), 4)
    aaf = [rng.uniform(0.5, 1.5) for _ in range(20)]
    tot = sum(aaf)
    aaf = [round(f / tot, 6) for f in aaf]
    aaf[-1] = round(1.0 - sum(aaf[:-1]), 6)
    lines = []
    for i in range(1, 20):
        lines.append(" ".join("%.4f" % exch[i][j] for j in range(i)))
    lines.append("")
    lines.append(" ".join("%.6f" % f for f in aaf))
    write(os.path.join(DATA, "custom_aa.dat"), "\n".join(lines) + "\n")
    t5 = Newick.parse("((p1:0.10,p2:0.20):0.05,p3:0.15,(p4:0.25,p5:0.12):0.08);")
    write(os.path.join(DATA, "oracle5.nwk"), t5.to_string() + "\n")
    q = custom_aa_q(exch, aaf)
    n5, s5 = simulate(t5, 40, AA, q, aaf, rng)
    write_phylip(os.path.join(DATA, "oracle_aa.phy"), n5, s5, interleaved=False)

    # --- Trees for dna_small --------------------------------------------
    write(os.path.join(DATA, "dna_small.user.nwk"),
          "((tax1:0.02,tax2:0.03):0.01,tax4:0.05,((tax5:0.04,tax6:0.03):0.02,"
          "((tax7:0.05,tax8:0.04):0.01,tax9:0.06):0.02):0.01);\n")
    write(os.path.join(DATA, "sim_jc_8.constraint.nwk"), "((t1,t3),t2,t4,t5,t6,t7,t8);\n")

    # --- Inputs with special content ------------------------------------
    dn = list(n8[:6]) + ["dup1", "dup2"]
    ds = [s[:200] for s in s8[:6]] + [s8[0][:200], s8[0][:200]]
    write_phylip(os.path.join(DATA, "duplicates.phy"), dn, ds, interleaved=False)

    # --- Malformed inputs ------------------------------------------------
    base = [s[:50] for s in s8[:5]]
    bn = n8[:5]

    def phy(names_, seqs_, header=None):
        head = header if header is not None else "%d %d" % (len(names_), len(seqs_[0]))
        return head + "\n" + "".join("%-10s %s\n" % (n, s) for n, s in zip(names_, seqs_))

    write(os.path.join(DATA, "bad_short_seq.phy"),
          phy(bn, base[:4] + [base[4][:40]], header="5 50"))
    write(os.path.join(DATA, "bad_truncated.phy"), phy(bn[:3], base[:3], header="5 50"))
    write(os.path.join(DATA, "bad_char.phy"),
          phy(bn, base[:4] + [base[4][:20] + "@" + base[4][21:]]))
    write(os.path.join(DATA, "bad_dup_names.phy"), phy([bn[0]] * 2 + bn[2:], base))
    write(os.path.join(DATA, "bad_header.phy"), phy(bn, base, header="five fifty"))
    write(os.path.join(DATA, "bad_too_few_taxa.phy"), phy(bn[:2], base[:2]))
    write(os.path.join(DATA, "empty.phy"), "")
    write(os.path.join(DATA, "bad_tree.nwk"), "((t1:0.1,t2:0.1),(t3:0.1,t4:0.1);\n")
    write(os.path.join(DATA, "bad_tree_taxa.nwk"),
          "((t1:0.1,t2:0.1):0.1,t3:0.1,(t4:0.1,zz:0.1):0.1);\n")

    # --- Larger alignments for the multithreading tests -----------------
    # Threads only engage above 128 site patterns per thread. Own generator,
    # so that adding these did not change the data sets above.
    write_threads_data(random.Random(20261015))
    write_structure_data(random.Random(20261016))
    print("data written to", DATA)


def random_tree(names, rng, lo=0.02, hi=0.15):
    """Random unrooted tree (trifurcating at the top) as a Newick string."""
    def build(sub):
        if len(sub) == 1:
            return sub[0]
        k = rng.randint(1, len(sub) - 1)
        return "(%s:%.3f,%s:%.3f)" % (build(sub[:k]), rng.uniform(lo, hi),
                                      build(sub[k:]), rng.uniform(lo, hi))
    names = list(names)
    rng.shuffle(names)
    a, b = len(names) // 3, 2 * len(names) // 3
    return "(%s:%.3f,%s:%.3f,%s:%.3f);" % (build(names[:a]), rng.uniform(lo, hi),
                                           build(names[a:b]), rng.uniform(lo, hi),
                                           build(names[b:]), rng.uniform(lo, hi))


def write_threads_data(rng):
    t = Newick.parse(random_tree(["s%d" % i for i in range(1, 25)], rng))
    write(os.path.join(DATA, "threads_dna.true.nwk"), t.to_string() + "\n")
    freqs = [0.3, 0.2, 0.25, 0.25]
    q = gtr_q([1.2, 3.5, 0.8, 1.1, 4.2, 1.0], freqs)
    n, seqs = simulate(t, 3000, DNA, q, freqs, rng, alpha=0.6)
    # A few ambiguous characters and gaps, to exercise every kernel branch.
    seqs = [list(x) for x in seqs]
    for k in range(300):
        seqs[rng.randrange(len(seqs))][rng.randrange(3000)] = "NRY-?"[k % 5]
    write_phylip(os.path.join(DATA, "threads_dna.phy"), n, ["".join(x) for x in seqs])

    t = Newick.parse(random_tree(["p%d" % i for i in range(1, 13)], rng, 0.05, 0.3))
    aaf = [1.0 / 20] * 20
    q = custom_aa_q([[0.0 if i == j else 1.0 for j in range(20)] for i in range(20)], aaf)
    n, seqs = simulate(t, 800, AA, q, aaf, rng, alpha=1.0)
    write_phylip(os.path.join(DATA, "threads_aa.phy"), n, seqs)

    # Binary characters: '-d generic' uses the generic (non-SIMD) kernels.
    t = Newick.parse(random_tree(["b%d" % i for i in range(1, 17)], rng, 0.1, 0.4))
    q = gtr_q_n(2)
    n, seqs = simulate(t, 3000, "01", q, [0.5, 0.5], rng)
    write_phylip(os.path.join(DATA, "threads_binary.phy"), n, seqs)


STRUCT_CLASSES = ["BUR_EXT", "BUR_HEL", "BUR_OTH", "EXP_EXT", "EXP_HEL", "EXP_OTH"]
STRUCT_RATES = [0.761816796788931, 0.744425646802117, 0.532457759429489,
                1.5639387472863, 2.06403411829438, 1.43336795177594]


def read_structure_matrices():
    """EX_EHO matrices (Le & Gascuel 2010) from the authors' AllModels.txt."""
    mats, name, buf = {}, None, []
    with open(os.path.join(DATA, "aa_matrices", "structure_AllModels.txt")) as fh:
        for line in fh:
            s = line.strip()
            if s and s.replace("_", "").isalpha() and s.isupper():
                if name:
                    mats[name] = buf
                name, buf = s, []
            elif name and s and not s.startswith("-"):
                buf += [float(x) for x in s.split()]
    mats[name] = buf
    return mats


def write_structure_data(rng):
    """Proteins evolving under the EX_EHO classes of Le & Gascuel (2010), with
    noisy annotations: 20% of the sites get a random class, 5% none."""
    t = Newick.parse(random_tree(["a%d" % i for i in range(1, 17)], rng, 0.05, 0.3))
    write(os.path.join(DATA, "struct_aa.true.nwk"), t.to_string() + "\n")
    mats = read_structure_matrices()
    props = [0.2, 0.15, 0.2, 0.15, 0.15, 0.15]
    nsite = 600
    columns = []
    for k, cname in enumerate(STRUCT_CLASSES):
        v = mats[cname]
        exch = [[0.0] * 20 for _ in range(20)]
        idx = 0
        for i in range(1, 20):
            for j in range(i):
                exch[i][j] = exch[j][i] = v[idx]
                idx += 1
        freqs = [f / sum(v[190:210]) for f in v[190:210]]
        q = [[x * STRUCT_RATES[k] for x in row] for row in custom_aa_q(exch, freqs)]
        n, seqs = simulate(t, int(round(props[k] * nsite)), AA, q, freqs, rng, alpha=0.8)
        for c in range(len(seqs[0])):
            columns.append((k, [sq[c] for sq in seqs]))
    rng.shuffle(columns)
    seqs = ["".join(col[1][i] for col in columns) for i in range(len(n))]
    write_phylip(os.path.join(DATA, "struct_aa.phy"), n, seqs)
    expo, sstr = [], []
    for k, _ in columns:
        u = rng.random()
        if u < 0.05:
            expo.append("?")
            sstr.append("?")
            continue
        if u < 0.25:
            k = rng.randrange(6)
        expo.append("BE"[k // 3])
        sstr.append("EHO"[k % 3])
    write(os.path.join(DATA, "struct_aa.annot"),
          "# EX_EHO classes of struct_aa.phy (20%% noise, 5%% unknown)\n"
          "exposure: %s\nstructure: %s\n" % ("".join(expo), "".join(sstr)))
    write(os.path.join(DATA, "struct_aa.classes"),
          "".join(str(k) for k, _ in columns) + "\n")


def gtr_q_n(ns):
    """Equal-rates (Mk) generator for ns states, normalised to rate 1."""
    off = 1.0 / (ns - 1)
    return [[-1.0 if i == j else off for j in range(ns)] for i in range(ns)]


if __name__ == "__main__":
    main()
