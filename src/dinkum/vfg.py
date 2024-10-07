"""encode view-from-genome rules.

X binds-and-upregulates Y
X binds-and-represses Y
X directly-or-indirectly-upregulates Y
X directly-or-indirectly-represses Y

X binds-and-upregulates Y if A else binds-and-represses
"""
from functools import total_ordering

_rules = []
_gene_names = []

def _add_rule(ix):
    _rules.append(ix)

def get_rules():
    return list(_rules)

def get_gene_names():
    return list(sorted(_gene_names))

def reset():
    global _rules
    global _gene_names
    _rules = []
    _gene_names = []


class Interactions:
    pass


class Interaction_Activates(Interactions):
    def __init__(self, *, source=None, dest=None, delay=1):
        assert isinstance(source, Gene)
        assert isinstance(dest, Gene)
        self.src = source
        self.dest = dest
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        The gene is active if its source was activate 'delay' ticks ago.
        """
        assert states
        assert tissue
        assert timepoint is not None

        if states.is_active(timepoint, self.delay, self.src, tissue):
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_Or(Interactions):
    def __init__(self, *, sources=None, dest=None, delay=1):
        for g in sources:
            assert isinstance(g, Gene)
        assert isinstance(dest, Gene)
        self.sources = sources
        self.dest = dest
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        The gene is active if any of its sources were activate 'delay'
        ticks ago.
        """
        assert states
        assert tissue

        source_active = [ states.is_active(timepoint, self.delay, g, tissue)
                          for g in self.sources ]

        if any(source_active):
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_AndNot(Interactions):
    def __init__(self, *, source=None, repressor=None, dest=None, delay=1):
        self.src = source
        self.repressor = repressor
        self.dest = dest
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        The gene is active if its activator was active 'delay' ticks ago,
        and its repressor was _not_ active then.
        """
        assert states
        assert tissue

        src_is_active = states.is_active(timepoint, self.delay,
                                         self.src, tissue)
        repressor_is_active = states.is_active(timepoint, self.delay,
                                              self.repressor, tissue)

        if src_is_active and not repressor_is_active:
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_And(Interactions):
    def __init__(self, *, sources=None, dest=None, delay=1):
        self.sources = sources
        self.dest = dest
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        The gene is active if all of its sources were active 'delay' ticks
        ago.
        """
        assert states
        assert tissue

        source_active = [ states.is_active(timepoint, self.delay, g, tissue)
                          for g in self.sources ]

        if all(source_active):
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_ToggleRepressed(Interactions):
    def __init__(self, *, tf=None, cofactor=None, dest=None, delay=1):
        self.tf = tf
        self.cofactor = cofactor
        self.dest = dest
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        The gene is active if the tf was active and the cofactor was active
        'delay' ticks ago.
        """
        assert states
        assert tissue

        tf_active = states.is_active(timepoint, self.delay,
                                     self.tf, tissue)
        cofactor_active = states.is_active(timepoint, self.delay,
                                           self.cofactor, tissue)


        if tf_active and not cofactor_active:
            yield self.dest, 1
        else:
            yield self.dest, 0


class Interaction_Ligand(Interactions):
    def __init__(self, *, activator=None, ligand=None, receptor=None, delay=1):
        self.activator = activator
        self.ligand = ligand
        self.receptor = receptor
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        A Ligand's next state is determined as follows:
        * its activator is ON
        * its ligand is currently ON in at least neighboring tissue
        """
        assert states
        assert tissue

        activator_is_active = states.is_active(timepoint, self.delay,
                                               self.activator, tissue)
        ligand_in_neighbors = []
        for neighbor in tissue.neighbors:
            neighbor_active = states.is_active(timepoint, self.delay,
                                               self.ligand, neighbor)
            ligand_in_neighbors.append(neighbor_active)

        activity = 0
        if activator_is_active and any(ligand_in_neighbors):
            activity = 1

        yield self.receptor, activity


class Gene:
    def __init__(self, *, name=None):
        global _gene_names

        assert name
        self.name = name

        _gene_names.append(name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name

    def __lt__(self, other):
        return self.name < other.name

    def __hash__(self):
        return hash(self.name)

    def active(self):           # present = active
        return 1

    def activated_by(self, *, source=None, delay=1):
        ix = Interaction_Activates(source=source, dest=self, delay=delay)
        _add_rule(ix)

    def activated_or(self, *, sources=None, delay=1):
        ix = Interaction_Or(sources=sources, dest=self, delay=delay)
        _add_rule(ix)

    def and_not(self, *, activator=None, repressor=None, delay=1):
        ix = Interaction_AndNot(source=activator, repressor=repressor,
                                dest=self, delay=delay)
        _add_rule(ix)

    def activated_by_and(self, *, sources, delay=1):
        ix = Interaction_And(sources=sources, dest=self, delay=delay)
        _add_rule(ix)

    def toggle_repressed(self, *, tf=None, cofactor=None, delay=1):
        ix = Interaction_ToggleRepressed(tf=tf, cofactor=cofactor,
                                         dest=self, delay=delay)
        _add_rule(ix)

    def is_present(self, *, where=None, start=None, duration=None):
        assert where
        assert start
        where.add_gene(gene=self, start=start, duration=duration)


class Receptor(Gene):
    def __init__(self, *, name=None):
        super().__init__(name=name)
        assert name
        self.name = name

    def ligand(self, *, activator=None, ligand=None):
        ix = Interaction_Ligand(activator=activator, ligand=ligand, receptor=self)
        _add_rule(ix)
