"""Reference likelihood of the structure models of Le & Gascuel (2010),
built on oracle.py, independently of PhyML's implementation.

Each partition (EX, EHO, EX_EHO) has one normalised matrix per class, with a
global rate; LG is added for CONF/LG. The site likelihood is

    L_s = sum_k w_sk * mean_g L(Q_k, rate_k * r_g | site s)

with r_g the discrete Gamma rates and site weights w_sk:
    MIX       P_k
    PART      1 for the annotated class (P_k if unknown)
    CONF/LG   chi for the annotated class, 1 - chi for LG (LG if unknown)
    CONF/MIX  chi [k = class] + (1 - chi) P_k (P_k if unknown)
(paper equations 4, 2, 3 and 5).
"""

import math
import os
import re

import oracle
from phymltest import data

MATRIX_FILE = data(os.path.join("aa_matrices", "structure_AllModels.txt"))

CLASSES = {
    "EX": ["BURIED", "EXPOSED"],
    "EHO": ["EXTENDED", "HELIX", "OTHER"],
    "EX_EHO": ["BUR_EXT", "BUR_HEL", "BUR_OTH", "EXP_EXT", "EXP_HEL", "EXP_OTH"],
}
RATES = {  # EX, EHO: the authors' binary, mean 1 under PROPS; EX_EHO: paper/IQ-TREE
    "BURIED": 0.556729, "EXPOSED": 1.359756,
    "EXTENDED": 0.857583, "HELIX": 1.163010, "OTHER": 0.932396,
    "BUR_EXT": 0.761816796788931, "BUR_HEL": 0.744425646802117,
    "BUR_OTH": 0.532457759429489, "EXP_EXT": 1.5639387472863,
    "EXP_HEL": 2.06403411829438, "EXP_OTH": 1.43336795177594,
}
PROPS = {  # class proportions of the authors' binary (phyml-mixSS)
    "EX": [0.448, 0.552],
    "EHO": [0.207793, 0.360558, 0.431649],
    "EX_EHO": [0.132640, 0.170755, 0.144605, 0.075153, 0.189803, 0.287044],
}
EXPOSURE = {"B": 0, "E": 1}
_CACHE = {}  # site likelihoods of each component, shared by the models
STRUCTURE = {"E": 0, "H": 1, "O": 2}


def read_matrices():
    mats, name, buf = {}, None, []
    with open(MATRIX_FILE) as fh:
        for line in fh:
            s = line.strip()
            if re.match(r"^[A-Z_]+$", s):
                if name:
                    mats[name] = buf
                name, buf = s, []
            elif name and s and not s.startswith("-"):
                buf += [float(x) for x in s.split()]
    mats[name] = buf
    out = {}
    for name, v in mats.items():
        exch = [[0.0] * 20 for _ in range(20)]
        k = 0
        for i in range(1, 20):
            for j in range(i):
                exch[i][j] = exch[j][i] = v[k]
                k += 1
        freqs = v[190:210]
        total = sum(freqs)
        out[name] = (exch, [f / total for f in freqs])
    return out


def lg_matrix():
    """LG as built into PhyML (src/init.c, Init_Qmat_LG)."""
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "src", "init.c")
    with open(path) as fh:
        src = fh.read()
    st = src.index("int Init_Qmat_LG(")
    body = src[st:src.index("\n}\n", st)]
    exch = [[0.0] * 20 for _ in range(20)]
    for a, b, v in re.findall(r"daa\[(\d+) \* 20 \+ (\d+)\]\s*=\s*([-0-9.Ee+]+);", body):
        exch[int(a)][int(b)] = exch[int(b)][int(a)] = float(v)
    pi = [0.0] * 20
    for a, v in re.findall(r"pi\[(\d+)\]\s*=\s*([-0-9.Ee+]+);", body):
        pi[int(a)] = float(v)
    total = sum(pi)
    return exch, [f / total for f in pi]


def site_classes(partition, exposure=None, structure=None):
    """Class index of each site, or -1, from annotation strings."""
    n = len(exposure or structure)
    out = []
    for i in range(n):
        e = EXPOSURE.get(exposure[i], -1) if exposure else -1
        h = STRUCTURE.get(structure[i], -1) if structure else -1
        if partition == "EX":
            out.append(e)
        elif partition == "EHO":
            out.append(h)
        else:
            out.append(-1 if e < 0 or h < 0 else e * 3 + h)
    return out


def weights(model, n_class, cls, chi, props):
    """Weight of each component (classes, then LG for CONF/LG) for a site."""
    if model == "mix":
        return list(props)
    if model == "part":
        return list(props) if cls < 0 else [1.0 if k == cls else 0.0 for k in range(n_class)]
    if model == "conf_lg":
        if cls < 0:
            return [0.0] * n_class + [1.0]
        return [chi if k == cls else 0.0 for k in range(n_class)] + [1.0 - chi]
    if model == "conf_mix":
        if cls < 0:
            return list(props)
        return [(chi if k == cls else 0.0) + (1.0 - chi) * props[k] for k in range(n_class)]
    raise ValueError(model)


def log_likelihood(tree, names, seqs, partition, model, classes, alpha, ncat,
                   chi=0.8, props=None, lg=None):
    """lg: (exch, freqs) of LG, needed for CONF/LG. Defaults: PhyML's starting
    values (chi = 0.8, published proportions)."""
    mats = read_matrices()
    names_c = CLASSES[partition]
    n_class = len(names_c)
    props = props or PROPS[partition]
    gamma = oracle.discrete_gamma_rates(alpha, ncat)
    comps = [(mats[c], RATES[c]) for c in names_c]
    if model == "conf_lg":
        comps.append((lg, 1.0))
    key0 = (tree.to_string(), tuple(names), tuple(seqs), alpha, ncat)
    per_comp = []
    for name, ((exch, freqs), rate) in zip(names_c + ["LG"], comps):
        key = key0 + (name,)
        if key not in _CACHE:
            q = oracle.custom_aa_q(exch, freqs)
            per_cat = [oracle.site_log_likelihoods(tree, names, seqs, q, freqs, oracle.AA,
                                                   rates=(rate * r,))
                       for r in gamma]
            _CACHE[key] = [sum(math.exp(c[s]) for c in per_cat) / ncat
                           for s in range(len(seqs[0]))]
        per_comp.append(_CACHE[key])
    total = 0.0
    for s in range(len(seqs[0])):
        w = weights(model, n_class, classes[s], chi, props)
        total += math.log(sum(w[k] * per_comp[k][s] for k in range(len(per_comp))))
    return total
