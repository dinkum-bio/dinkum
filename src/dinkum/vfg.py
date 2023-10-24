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

    def advance(self, *, state=None, tissue=None):
        assert state
        assert tissue

        this_active = state[tissue]
        if self.src in this_active:
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_Or(Interactions):
    def __init__(self, *, sources=None, dest=None):
        for g in sources:
            assert isinstance(g, Gene)
        assert isinstance(dest, Gene)
        self.sources = sources
        self.dest = dest

    def advance(self, *, state=None, tissue=None):
        assert state
        assert tissue

        this_active = state[tissue]
        is_active = 0
        for g in self.sources:
            if g in this_active:
                is_active = 1

        yield self.dest, is_active


class Interaction_AndNot(Interactions):
    def __init__(self, *, source=None, repressor=None, dest=None):
        self.src = source
        self.repressor = repressor
        self.dest = dest

    def advance(self, *, state=None, tissue=None):
        assert state
        assert tissue

        this_active = state[tissue]
        if self.src in this_active and not self.repressor in this_active:
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_And(Interactions):
    def __init__(self, *, sources=None, dest=None):
        self.sources = sources
        self.dest = dest

    def advance(self, *, state=None, tissue=None):
        assert state
        assert tissue

        this_active = state[tissue]
        if all([ g in this_active for g in self.sources ]):
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_ToggleRepressed(Interactions):
    def __init__(self, *, tf=None, cofactor=None, dest=None):
        self.tf = tf
        self.cofactor = cofactor
        self.dest = dest

    def advance(self, *, state=None, tissue=None):
        assert state
        assert tissue

        this_active = state[tissue]
        activity = 0
        # tf has to be present...
        if self.tf in this_active:
            # if cofactor is present, repress.
            if self.cofactor not in this_active:
                activity = 1

        yield self.dest, activity


class Interaction_Ligand(Interactions):
    def __init__(self, *, activator=None, ligand=None, receptor=None):
        self.activator = activator
        self.ligand = ligand
        self.receptor = receptor

    def advance(self, *, state=None, tissue=None):
        assert state
        assert tissue

        activity = 0
        if self.activator in state[tissue]:
            #print(f'XXX {self.receptor.name} is present in {tissue.name} b/c {self.activator.name} is active')
            for neighbor in tissue.neighbors:
                neighbor_active = state[neighbor]
                if self.ligand in neighbor_active:
                    #print(f'XXX {self.receptor.name} is ACTIVE in {tissue.name} b/c {self.ligand.name} is in {neighbor.name}')
                    activity = 1

        yield self.receptor, activity


class Gene:
    def __init__(self, *, name=None):
        assert name
        self.name = name

    def active(self):           # present = active
        return 1

    def activated_by(self, *, source=None):
        ix = Interaction_Activates(source=source, dest=self)
        _add_rule(ix)

    def activated_or(self, *, sources=None):
        ix = Interaction_Or(sources=sources, dest=self)
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
        where.add_gene(gene=self, start=start, duration=duration)


class Receptor(Gene):
    def __init__(self, *, name=None):
        assert name
        self.name = name

    def ligand(self, *, activator=None, ligand=None):
        ix = Interaction_Ligand(activator=activator, ligand=ligand, receptor=self)
        _add_rule(ix)
