"""Reference likelihood of LG4M and LG4X (Le, Dang & Gascuel 2012), built on
oracle.py, from the authors' matrix file (data/aa_matrices/LG4X_4M.txt).

    LG4M (eq. 5)  L_s = 1/4 sum_k L(x_s | Q_k, gamma_k)   (Gamma means, alpha)
    LG4X (eq. 6)  L_s = sum_k w_k L(x_s | Q_k, rho_k),    sum w = sum w rho = 1

Q_k uses the exchangeabilities and frequencies of matrix k (M1-M4 or X1-X4,
the first one the slowest), normalised to one substitution per unit time;
the published frequencies are divided by their sum.
"""

import math
import os
import re

import oracle
from phymltest import data

MATRIX_FILE = data(os.path.join("aa_matrices", "LG4X_4M.txt"))


def read_matrices():
    with open(MATRIX_FILE) as fh:
        lines = fh.read().splitlines()
    freqs, exch, cur = {}, {}, None
    for line in lines:
        parts = line.split()
        if len(parts) == 21 and re.fullmatch(r"[XM][1-4]", parts[0]):
            freqs[parts[0]] = [float(x) for x in parts[1:]]
        elif re.fullmatch(r"[XM][1-4]", line.strip()):
            cur = line.strip()
            exch[cur] = []
        elif cur and parts and len(exch[cur]) < 190:
            exch[cur] += [float(x) for x in parts]
    out = {}
    for name, v in exch.items():
        assert len(v) == 190, name
        e = [[0.0] * 20 for _ in range(20)]
        k = 0
        for i in range(1, 20):
            for j in range(i):
                e[i][j] = e[j][i] = v[k]
                k += 1
        total = sum(freqs[name])
        out[name] = (e, [f / total for f in freqs[name]])
    return out


def log_likelihood(tree, names, seqs, model, rates, weights):
    """model: 'LG4M' or 'LG4X'; rates and weights of the four categories."""
    mats = read_matrices()
    prefix = "M" if model == "LG4M" else "X"
    per_cat = []
    for k in range(4):
        exch, freqs = mats["%s%d" % (prefix, k + 1)]
        q = oracle.custom_aa_q(exch, freqs)
        per_cat.append(oracle.site_log_likelihoods(tree, names, seqs, q, freqs, oracle.AA,
                                                   rates=(rates[k],)))
    return sum(math.log(sum(w * math.exp(c[s]) for w, c in zip(weights, per_cat)))
               for s in range(len(seqs[0])))
