"""encode view-from-genome rules.

X binds-and-upregulates Y
X binds-and-represses Y
X directly-or-indirectly-upregulates Y
X directly-or-indirectly-represses Y

X binds-and-upregulates Y if A else binds-and-represses
"""

_rules = []
def _add_rule(ix):
    _rules.append(ix)

def get_rules():
    return list(_rules)

def reset():
    global _rules
    _rules = []


class Interactions:
    pass


class Interaction_Activates(Interactions):
    def __init__(self, *, source=None, dest=None):
        assert isinstance(source, Gene)
        assert isinstance(dest, Gene)
        self.src = source
        self.dest = dest

    def advance(self, *, present=[]):
        #print('XXX', self.src.name, self.dest.name, [p.name for p in present])
        if self.src in present:
            yield self.dest


class Interaction_AndNot(Interactions):
    def __init__(self, *, source=None, repressor=None, dest=None):
        assert isinstance(source, Gene)
        assert isinstance(repressor, Gene)
        assert isinstance(dest, Gene)
        self.src = source
        self.repressor = repressor
        self.dest = dest

    def advance(self, *, present=[]):
        #print('XXX', self.src.name, self.dest.name, [p.name for p in present])
        if self.src in present and not self.repressor in present:
            yield self.dest


class Gene:
    def __init__(self, *, name=None):
        assert name
        self.name = name

    def activated_by(self, *, source=None):
        assert isinstance(source, Gene)

        ix = Interaction_Activates(source=source, dest=self)
        _add_rule(ix)
    

    def and_not(self, *, activator=None, repressor=None):
        assert isinstance(repressor, Gene)
        assert isinstance(activator, Gene)

        ix = Interaction_AndNot(source=activator, repressor=repressor,
                                dest=self)
        _add_rule(ix)

    def is_present(self, *, where=None, start=None, duration=None):
        assert where
        assert start
        where.add(gene=self, start=start, duration=duration)
