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
        if self.src in present:
            yield self.dest


class Interaction_AndNot(Interactions):
    def __init__(self, *, source=None, repressor=None, dest=None):
        self.src = source
        self.repressor = repressor
        self.dest = dest

    def advance(self, *, present=[]):
        if self.src in present and not self.repressor in present:
            yield self.dest


class Interaction_And(Interactions):
    def __init__(self, *, sources=None, dest=None):
        self.sources = sources
        self.dest = dest

    def advance(self, *, present=[]):
        if all([ g in present for g in self.sources ]):
            yield self.dest


class Interaction_ToggleRepressed(Interactions):
    def __init__(self, *, tf=None, cofactor=None, dest=None):
        self.tf = tf
        self.cofactor = cofactor
        self.dest = dest

    def advance(self, *, present=[]):
        # tf has to be present...
        if self.tf in present:
            # if cofactor is present, repress.
            if self.cofactor not in present:
                yield self.dest


class Gene:
    def __init__(self, *, name=None):
        assert name
        self.name = name

    def activated_by(self, *, source=None):
        ix = Interaction_Activates(source=source, dest=self)
        _add_rule(ix)

    def and_not(self, *, activator=None, repressor=None):
        ix = Interaction_AndNot(source=activator, repressor=repressor,
                                dest=self)
        _add_rule(ix)

    def activated_by_and(self, *, sources):
        ix = Interaction_And(sources=sources, dest=self)
        _add_rule(ix)

    def toggle_repressed(self, *, tf=None, cofactor=None):
        ix = Interaction_ToggleRepressed(tf=tf, cofactor=cofactor,
                                         dest=self)
        _add_rule(ix)

    def is_present(self, *, where=None, start=None, duration=None):
        assert where
        assert start
        where.add(gene=self, start=start, duration=duration)
