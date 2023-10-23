"""encode view-from-nucleus rules

X is-present in celltype/tissue/compartment M at time T
"""

_tissues = []
def _add_tissue(t):
    global _tissues
    _tissues.append(t)

def get_tissues():
    return list(_tissues)

def reset():
    global _tissues
    _tissues = []


class Tissue:
    def __init__(self, *, name=None):
        assert name
        self.name = name
        self.present = []
        self.neighbors = set()

        _add_tissue(self)

    def add_gene(self, *, gene=None, start=None, duration=None):
        assert gene
        assert start is not None
        self.present.append((gene, start, duration))

    def add_neighbor(self, *, neighbor=None):
        assert neighbor
        self.neighbors.add(neighbor)

    def all_active(self, *, at=None):
        assert at is not None
        for gene, start, duration in self.present:
            if at >= start:
                if duration is None:
                    yield gene
                elif at <= start + duration:
                    yield gene
