"""encode view-from-genome rules.

X binds-and-upregulates Y
X binds-and-represses Y
X directly-or-indirectly-upregulates Y
X directly-or-indirectly-represses Y

X binds-and-upregulates Y if A else binds-and-represses
"""
from functools import total_ordering
import inspect

_rules = []
_genes = []

def _add_rule(ix):
    _rules.append(ix)

def get_rules():
    return list(_rules)

def get_gene_names():
    return [ g.name for g in sorted(_genes) ]

def reset():
    global _rules
    global _genes
    _rules = []
    _genes = []


class Interactions:
    pass


class Interaction_Activates(Interactions):
    def __init__(self, *, source=None, dest=None, delay=1):
        assert isinstance(source, Gene), f"'{source}' must be a Gene (but is not)"
        assert isinstance(dest, Gene), f"'{dest}' must be a Gene (but is not)"
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
            assert isinstance(g, Gene), f"source '{g}' must be a Gene"
        assert isinstance(dest, Gene), f"dest '{dest}' must be a Gene"
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


class Interaction_Arbitrary(Interactions):
    def __init__(self, *, dest=None, state_fn=None, delay=1):
        assert dest
        assert state_fn

        self.dest = dest
        self.state_fn = state_fn
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        assert states
        assert tissue

        dep_gene_names = inspect.getfullargspec(self.state_fn).args
        dep_genes = []
        for name in dep_gene_names:
            found = False
            for g in _genes:
                if g.name == name:
                    dep_genes.append(g)
                    found = True
                    break
            if not found:
                raise Exception(f"no such gene: '{name}'")

        dep_state = [ states.is_active(timepoint, self.delay, g, tissue)
                      for g in dep_genes ]
        is_active = self.state_fn(*dep_state)

        if is_active:
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
        * its ligand is currently ON in at least one neighboring tissue
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
        global _genes

        assert name, "Gene must have a name"
        self.name = name

        _genes.append(self)
        self._set_ligand = None

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name

    def __lt__(self, other):
        return self.name < other.name

    def __hash__(self):
        return hash(self.name)

    def present(self):
        return 1

    def ligand_present(self):
        return (self._set_ligand is None) or 1 # @CTB

    def active(self):           # present = active
        return self.present() and self.ligand_present()

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

    def custom_activation(self, *, state_fn=None, delay=1):
        ix = Interaction_Arbitrary(dest=self, state_fn=state_fn, delay=delay)
        _add_rule(ix)


class Receptor(Gene):
    def __init__(self, *, name=None, ligand=None):
        super().__init__(name=name)
        assert name
        self._set_ligand = ligand

    def ligand(self, *, activator=None, ligand=None):
        if ligand is None:
            ligand = self._set_ligand
            if ligand is None:
                raise Exception("need to specify a ligand for this receptor, either at creation or here")

        ix = Interaction_Ligand(activator=activator, ligand=ligand, receptor=self)
        _add_rule(ix)

    def activated_by(self, *, activator=None, ligand=None):
        if ligand is None:
            ligand = self._set_ligand
            if ligand is None:
                raise Exception("need to specify a ligand for this receptor, either at creation or here")

        ix = Interaction_Ligand(activator=activator, ligand=ligand, receptor=self)
        _add_rule(ix)
