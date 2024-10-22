"""encode view-from-nucleus rules

X is-present in celltype/tissue/compartment M at time T
"""
from functools import total_ordering

from . import vfg

_tissues = []
def _add_tissue(t):
    global _tissues
    _tissues.append(t)

def get_tissues():
    return list(sorted(_tissues))

def get_tissue_names():
    return list(sorted([ t.name for t in _tissues ]))

def reset():
    global _tissues
    _tissues = []


@total_ordering
class Tissue:
    def __init__(self, *, name=None):
        assert name, "Tissue must have a name"
        self.name = name
        self.neighbors = set()

        _add_tissue(self)

    def __repr__(self):
        return f"Tissue('{self.name}')"

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name

    def __lt__(self, other):
        return self.name < other.name

    def __hash__(self):
        return hash(self.name)

    def add_gene(self, *, gene=None, start=None, duration=None):
        assert gene
        assert gene in vfg._genes
        assert start is not None
        gene.is_present(start=start, duration=duration, where=self)

    def add_neighbor(self, *, neighbor=None):
        assert neighbor
        self.neighbors.add(neighbor)
