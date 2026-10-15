"""Independent reference implementation of the phylogenetic likelihood.

Written from the textbook definitions (Felsenstein pruning, time-reversible
Markov models, Yang's discrete Gamma) without looking at PhyML's numerical
code, so that agreement with PhyML is meaningful. Standard library only.
"""

import math

from newick import Newick  # noqa: F401  (re-exported for tools/make_data.py)

DNA = "ACGT"
AA = "ARNDCQEGHILKMFPSTWYV"

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG", "W": "AT",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT", "X": "ACGT", "-": "ACGT", "?": "ACGT", ".": "ACGT",
}


def tip_vector(char, alphabet):
    char = char.upper()
    if alphabet == DNA:
        allowed = IUPAC[char]
        return [1.0 if s in allowed else 0.0 for s in DNA]
    if char in AA:
        return [1.0 if s == char else 0.0 for s in AA]
    if char in "X-?":
        return [1.0] * 20
    if char == "B":
        return [1.0 if s in "ND" else 0.0 for s in AA]
    if char == "Z":
        return [1.0 if s in "QE" else 0.0 for s in AA]
    raise ValueError("unknown amino-acid character " + char)


# ----------------------------------------------------------------- models
def normalise(q, freqs):
    n = len(freqs)
    for i in range(n):
        q[i][i] = 0.0
        q[i][i] = -sum(q[i])
    mu = -sum(freqs[i] * q[i][i] for i in range(n))
    return [[x / mu for x in row] for row in q]


def gtr_q(rates, freqs):
    """rates in the order AC, AG, AT, CG, CT, GT."""
    ac, ag, at, cg, ct, gt = rates
    r = [[0, ac, ag, at], [ac, 0, cg, ct], [ag, cg, 0, gt], [at, ct, gt, 0]]
    q = [[r[i][j] * freqs[j] for j in range(4)] for i in range(4)]
    return normalise(q, freqs)


def custom_aa_q(exch, freqs):
    q = [[exch[i][j] * freqs[j] for j in range(20)] for i in range(20)]
    return normalise(q, freqs)


def read_paml_matrix(path):
    """Read a PAML-format amino-acid exchangeability file."""
    with open(path) as fh:
        nums = [float(x) for x in fh.read().split()]
    exch = [[0.0] * 20 for _ in range(20)]
    k = 0
    for i in range(1, 20):
        for j in range(i):
            exch[i][j] = exch[j][i] = nums[k]
            k += 1
    freqs = nums[k:k + 20]
    return exch, freqs


def matmul(a, b):
    n = len(a)
    bt = list(zip(*b))
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def expm(q, t):
    """exp(Q t) by scaling and squaring of a Taylor series."""
    n = len(q)
    a = [[x * t for x in row] for row in q]
    norm = max(sum(abs(x) for x in row) for row in a)
    s = max(0, int(math.ceil(math.log2(norm))) + 4) if norm > 0 else 0
    a = [[x / (2 ** s) for x in row] for row in a]
    result = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    term = [row[:] for row in result]
    for k in range(1, 20):
        term = matmul(term, a)
        term = [[x / k for x in row] for row in term]
        result = [[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(result, term)]
    for _ in range(s):
        result = matmul(result, result)
    return result


# ------------------------------------------------------- gamma functions
def _lower_reg_gamma(a, x):
    """Regularised lower incomplete gamma P(a, x)."""
    if x <= 0:
        return 0.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        term = total = 1.0 / a
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-16:
                break
        return total * math.exp(-x + a * math.log(x) - gln)
    # continued fraction for Q(a, x) (modified Lentz)
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return 1.0 - math.exp(-x + a * math.log(x) - gln) * h


def gamma_quantile(p, shape, rate):
    lo, hi = 0.0, 1.0
    while _lower_reg_gamma(shape, hi * rate) < p:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _lower_reg_gamma(shape, mid * rate) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def discrete_gamma_rates(alpha, ncat, median=False):
    """Yang (1994) equal-probability discrete Gamma rates with mean 1."""
    if ncat == 1:
        return [1.0]
    if median:
        r = [gamma_quantile((2 * i + 1) / (2.0 * ncat), alpha, alpha)
             for i in range(ncat)]
        m = sum(r) / ncat
        return [x / m for x in r]
    cuts = [gamma_quantile((i + 1.0) / ncat, alpha, alpha) for i in range(ncat - 1)]
    cdf = [0.0] + [_lower_reg_gamma(alpha + 1.0, c * alpha) for c in cuts] + [1.0]
    return [ncat * (cdf[i + 1] - cdf[i]) for i in range(ncat)]


# ------------------------------------------------------------ likelihood
def site_log_likelihoods(tree, names, seqs, q, freqs, alphabet,
                         rates=(1.0,), pinv=0.0, invar_rescale=True):
    """Per-site log-likelihoods of an alignment on a fixed tree.

    rates:   relative rates of equiprobable categories (mean 1).
    pinv:    proportion of invariable sites.
    invar_rescale: if True, rates of the variable classes are divided by
             (1 - pinv) so that branch lengths keep their meaning of expected
             substitutions per site over all sites.
    """
    seq_of = dict(zip(names, seqs))
    nsite = len(seqs[0])
    ns = len(alphabet)
    scale = 1.0 / (1.0 - pinv) if (pinv > 0 and invar_rescale) else 1.0
    pm = {}
    for node in tree.nodes():
        if node is tree.root:
            continue
        for r in rates:
            pm[(id(node), r)] = expm(q, node.length * r * scale)

    tips = {n: [tip_vector(c, alphabet) for c in s] for n, s in seq_of.items()}

    def partial(node, site, r):
        if node.is_leaf():
            return tips[node.name][site]
        out = [1.0] * ns
        for child in node.children:
            v = partial(child, site, r)
            p = pm[(id(child), r)]
            for i in range(ns):
                out[i] *= sum(p[i][j] * v[j] for j in range(ns))
        return out

    result = []
    for site in range(nsite):
        lk_var = 0.0
        for r in rates:
            v = partial(tree.root, site, r)
            lk_var += sum(freqs[i] * v[i] for i in range(ns)) / len(rates)
        lk = (1.0 - pinv) * lk_var
        if pinv > 0:
            inter = None
            for n in names:
                vec = tips[n][site]
                inter = vec if inter is None else [a * b for a, b in zip(inter, vec)]
            lk += pinv * sum(freqs[i] * inter[i] for i in range(ns))
        result.append(math.log(lk))
    return result


def ancestral_marginals(tree, names, seqs, q, freqs, alphabet, clade):
    """Marginal posterior state probabilities, per site, at the internal node
    whose descendants (in `tree` as written) are exactly `clade`.
    Single rate class, no invariable sites."""
    seq_of = dict(zip(names, seqs))
    ns = len(alphabet)
    # undirected adjacency: node -> [(neighbour, branch length)]
    adj = {}
    for node in tree.nodes():
        adj.setdefault(id(node), [])
        for c in node.children:
            adj[id(node)].append((c, c.length))
            adj.setdefault(id(c), []).append((node, c.length))
    below = {}

    def leaves_below(n):
        if n.is_leaf():
            below[id(n)] = frozenset([n.name])
        else:
            below[id(n)] = frozenset().union(*(leaves_below(c) for c in n.children))
        return below[id(n)]
    leaves_below(tree.root)
    target = [n for n in tree.internal_nodes() if below[id(n)] == frozenset(clade)]
    if len(target) != 1:
        raise ValueError("clade does not match exactly one internal node")
    target = target[0]
    pm = {}

    def pmat(length):
        if length not in pm:
            pm[length] = expm(q, length)
        return pm[length]

    def partial(node, came_from, site):
        if node.is_leaf() and came_from is not None:
            return tip_vector(seq_of[node.name][site], alphabet)
        out = [1.0] * ns
        for nb, length in adj[id(node)]:
            if nb is came_from:
                continue
            v = partial(nb, node, site)
            p = pmat(length)
            for i in range(ns):
                out[i] *= sum(p[i][j] * v[j] for j in range(ns))
        return out

    result = []
    for site in range(len(seqs[0])):
        v = partial(target, None, site)
        w = [freqs[i] * v[i] for i in range(ns)]
        tot = sum(w)
        result.append([x / tot for x in w])
    return result


def log_likelihood(*args, **kwargs):
    return sum(site_log_likelihoods(*args, **kwargs))


def read_phylip(path):
    """Read a sequential or interleaved PHYLIP file into (names, seqs)."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    ntax, nsite = map(int, lines[0].split())
    names, seqs = [], []
    for line in lines[1:1 + ntax]:
        parts = line.split()
        names.append(parts[0])
        seqs.append("".join(parts[1:]))
    k = 1 + ntax
    while any(len(s) < nsite for s in seqs):
        for i in range(ntax):
            seqs[i] += lines[k].replace(" ", "")
            k += 1
    return names, seqs


def empirical_freqs(seqs, alphabet):
    counts = [0] * len(alphabet)
    for s in seqs:
        for c in s:
            if c in alphabet:
                counts[alphabet.index(c)] += 1
    tot = float(sum(counts))
    return [c / tot for c in counts]
