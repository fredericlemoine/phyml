"""Minimal Newick parser used by the PhyML tests (standard library only)."""

import re


class Node:
    __slots__ = ("name", "length", "label", "children", "parent")

    def __init__(self, name=None, length=0.0, label=None):
        self.name = name
        self.length = length
        self.label = label          # internal-node label (e.g. a support value)
        self.children = []
        self.parent = None

    def is_leaf(self):
        return not self.children


class Newick:
    """A tree read from a Newick string. The root may have 2 or 3 children."""

    def __init__(self, root):
        self.root = root

    # ------------------------------------------------------------------ parse
    @classmethod
    def parse(cls, text):
        text = re.sub(r"\[[^\]]*\]", "", text.strip())  # drop [comments]
        if not text.endswith(";"):
            raise ValueError("Newick string must end with ';'")
        pos = 0

        def parse_label():
            nonlocal pos
            start = pos
            if pos < len(text) and text[pos] == "'":
                end = text.index("'", pos + 1)
                pos = end + 1
                return text[start + 1:end]
            while pos < len(text) and text[pos] not in ",():;":
                pos += 1
            return text[start:pos].strip()

        def parse_length():
            nonlocal pos
            if pos < len(text) and text[pos] == ":":
                pos += 1
                start = pos
                while pos < len(text) and text[pos] not in ",();":
                    pos += 1
                return float(text[start:pos])
            return None

        def parse_node():
            nonlocal pos
            node = Node()
            if text[pos] == "(":
                pos += 1
                while True:
                    child = parse_node()
                    child.parent = node
                    node.children.append(child)
                    if text[pos] == ",":
                        pos += 1
                        continue
                    if text[pos] == ")":
                        pos += 1
                        break
                    raise ValueError("unexpected %r at %d" % (text[pos], pos))
                lab = parse_label()
                node.label = lab or None
            else:
                node.name = parse_label()
            length = parse_length()
            node.length = 0.0 if length is None else length
            return node

        root = parse_node()
        if text[pos] != ";":
            raise ValueError("trailing characters in Newick string")
        return cls(root)

    @classmethod
    def read(cls, path):
        with open(path) as fh:
            return cls.parse(fh.read().strip().split("\n")[0])

    # ---------------------------------------------------------------- output
    def to_string(self):
        def rec(node):
            if node.is_leaf():
                s = node.name
            else:
                s = "(" + ",".join(rec(c) for c in node.children) + ")"
                if node.label:
                    s += node.label
            if node is not self.root:
                s += ":%.6g" % node.length
            return s
        return rec(self.root) + ";"

    # ------------------------------------------------------------- queries
    def nodes(self):
        out, stack = [], [self.root]
        while stack:
            n = stack.pop()
            out.append(n)
            stack.extend(reversed(n.children))
        return out

    def leaves(self):
        return [n for n in self.nodes() if n.is_leaf()]

    def leaf_names(self):
        return [n.name for n in self.leaves()]

    def internal_nodes(self):
        return [n for n in self.nodes() if not n.is_leaf()]

    def total_length(self):
        return sum(n.length for n in self.nodes() if n is not self.root)

    def splits(self):
        """Non-trivial bipartitions, each as the frozenset side that excludes
        the alphabetically first taxon. A rooted tree is treated as unrooted."""
        taxa = sorted(self.leaf_names())
        ref = taxa[0]
        full = frozenset(taxa)
        out = {}

        def rec(node):
            if node.is_leaf():
                return frozenset([node.name])
            s = frozenset().union(*(rec(c) for c in node.children))
            if node is not self.root:
                side = s if ref not in s else full - s
                if 1 < len(side) < len(full) - 1:
                    out[side] = node
            return s
        rec(self.root)
        return out

    def has_split(self, clade):
        """True if the tree has the bipartition separating `clade` from the rest."""
        clade = frozenset(clade)
        full = frozenset(self.leaf_names())
        ref = sorted(full)[0]
        side = clade if ref not in clade else full - clade
        return side in self.splits()

    def supports(self):
        """Map split -> numeric label of the corresponding internal node."""
        res = {}
        for split, node in self.splits().items():
            if node.label is not None:
                try:
                    res[split] = float(node.label)
                except ValueError:
                    pass
        return res


def rf_distance(t1, t2):
    """Unnormalised Robinson-Foulds distance between two trees."""
    s1, s2 = set(t1.splits()), set(t2.splits())
    return len(s1 ^ s2)


def transfer_support(ref, replicates):
    """Transfer Bootstrap Expectation (Lemoine et al. 2018) of each split of
    `ref`, computed by brute force from the replicate trees."""
    taxa = frozenset(ref.leaf_names())
    n = len(taxa)
    out = {}
    rep_splits = []
    for t in replicates:
        sp = set(t.splits())
        sp.update(frozenset([x]) for x in taxa)          # trivial splits too
        rep_splits.append(sp)
    for split in ref.splits():
        p = min(len(split), n - len(split))
        total = 0
        for sp in rep_splits:
            best = p - 1
            for other in sp:
                d = len(split ^ other)
                best = min(best, d, n - d)
            total += best
        out[split] = 1.0 - (total / len(replicates)) / (p - 1)
    return out
