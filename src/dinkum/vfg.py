"""encode view-from-genome rules.

X binds-and-upregulates Y
X binds-and-represses Y
X directly-or-indirectly-upregulates Y
X directly-or-indirectly-represses Y

X binds-and-upregulates Y if A else binds-and-represses
"""
from functools import total_ordering
import inspect
import collections

GeneStateInfo = collections.namedtuple('GeneStateInfo', ['level', 'active'])
DEFAULT_OFF = GeneStateInfo(level=0, active=False)

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


def _retrieve_ligands(timepoint, states, tissue, delay):
    "Retrieve all ligands in neighboring tissues for the given timepoint/delay"
    #assert isinstance(tissue, Tissue)

    ligands = set()
    for gene in _genes:
        if gene._is_ligand:
            for neighbor in tissue.neighbors:
                if states.is_active(timepoint, delay, gene, neighbor):
                    ligands.add(gene)

    return ligands


class Interactions:
    multiple_allowed = False

    def check_ligand(self, timepoint, states, tissue, delay):
        if getattr(self.dest, '_set_ligand', None):
            #print(self.dest, "is a receptor w/ a ligand")
            ligands_in_neighbors = _retrieve_ligands(timepoint, states,
                                                     tissue, delay)
            if self.dest._set_ligand in ligands_in_neighbors:
                return True
            return False
        else:
            #print(self.dest, "is not a receptor")
            return True         # by default, not ligand => is active


class Interaction_IsPresent(Interactions):
    multiple_allowed = True

    def __init__(self, *, dest=None, start=None, duration=None, tissue=None,
                 level=None, decay=None):
        assert isinstance(dest, Gene), f"'{dest}' must be a Gene (but is not)"
        assert start is not None, "must provide start time"
        assert level is not None, "must provide level"
        assert decay is not None, "must provide decay"
        assert decay >= 1
        assert decay < 1e6
        assert tissue
        self.dest = dest
        self.tissue = tissue
        self.start = start
        self.duration = duration
        self.level = level
        self.decay = decay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        # ignore states
        if tissue == self.tissue:
            if timepoint >= self.start:
                if self.duration is None or \
                   timepoint < self.start + self.duration: # active!
                    if self.check_ligand(timepoint, states, tissue, delay=1):
                        yield self.dest, GeneStateInfo(level=self.level,
                                                       active=True)
                    else:
                        yield self.dest, GeneStateInfo(level=self.level,
                                                       active=False)
        # we have no opinion on activity outside our tissue!

        #print('XXX', self.level, timepoint)
        self.level = round(self.level / self.decay + 0.5)


class Interaction_Activates(Interactions):
    def __init__(self, *, source=None, dest=None, delay=1):
        assert isinstance(source, Gene), f"'{source}' must be a Gene (but is not)"
        assert isinstance(dest, Gene), f"'{dest}' must be a Gene (but is not)"
        self.src = source
        self.dest = dest
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        """
        The gene is active if its source was active 'delay' ticks ago.
        """
        if not states:
            return

        assert tissue
        assert timepoint is not None

        if states.is_active(timepoint, self.delay, self.src, tissue):
            is_active = self.check_ligand(timepoint, states, tissue, self.delay)
            yield self.dest, GeneStateInfo(level=100, active=is_active)
        else:
            yield self.dest, GeneStateInfo(level=0, active=False)


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
        if not states:
            return

        assert tissue

        source_active = [ states.is_active(timepoint, self.delay, g, tissue)
                          for g in self.sources ]

        if any(source_active):
            is_active = self.check_ligand(timepoint, states, tissue, self.delay)
            yield self.dest, GeneStateInfo(level=100, active=is_active)
        else:
            yield self.dest, GeneStateInfo(level=0, active=False)


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
        if not states:
            return

        assert tissue

        src_is_active = states.is_active(timepoint, self.delay,
                                         self.src, tissue)
        repressor_is_active = states.is_active(timepoint, self.delay,
                                              self.repressor, tissue)

        if src_is_active and not repressor_is_active:
            is_active = self.check_ligand(timepoint, states, tissue, self.delay)
            yield self.dest, GeneStateInfo(level=100, active=is_active)
        else:
            yield self.dest, GeneStateInfo(level=0, active=False)


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
        if not states:
            return

        assert tissue

        source_active = [ states.is_active(timepoint, self.delay, g, tissue)
                          for g in self.sources ]

        if all(source_active):
            is_active = self.check_ligand(timepoint, states, tissue, self.delay)
            yield self.dest, GeneStateInfo(level=100, active=is_active)
        else:
            yield self.dest, GeneStateInfo(level=0, active=False)


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
        if not states:
            return

        assert tissue

        tf_active = states.is_active(timepoint, self.delay,
                                     self.tf, tissue)
        cofactor_active = states.is_active(timepoint, self.delay,
                                           self.cofactor, tissue)


        # @CTB refigure for receptor/ligand, yah?
        if tf_active and not cofactor_active:
            is_active = self.check_ligand(timepoint, states, tissue, self.delay)
            yield self.dest, GeneStateInfo(level=100, active=is_active)
        else:
            yield self.dest, GeneStateInfo(level=0, active=False)


class Interaction_Arbitrary(Interactions):
    def __init__(self, *, dest=None, state_fn=None, delay=1):
        assert dest
        assert state_fn

        self.dest = dest
        self.state_fn = state_fn
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        if not states:
            return

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

        if is_active and self.check_ligand(timepoint,
                                           states,
                                           tissue,
                                           self.delay):
            yield self.dest, GeneStateInfo(level=100, active=True)
        else:
            yield self.dest, GeneStateInfo(level=0, active=False)


class Interaction_ArbitraryComplex(Interactions):
    """
    An interaction that supports arbitrary logic, + levels.
    """
    def __init__(self, *, dest=None, state_fn=None, delay=1):
        assert dest
        assert state_fn

        self.dest = dest
        self.state_fn = state_fn
        self.delay = delay

    def advance(self, *, timepoint=None, states=None, tissue=None):
        # 'states' is class States...
        # @CTB refactor for presence/activity
        if not states:
            return

        assert tissue

        # get the names of the genes on the function => and their activity
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

        # pass in their full GeneStateInfo
        delay = self.delay
        dep_state = [ states.get_gene_state_info(timepoint, delay, g, tissue)
                      for g in dep_genes ]

        level, is_active = self.state_fn(*dep_state)

        if is_active:
            is_active = self.check_ligand(timepoint,
                                          states,
                                          tissue,
                                          self.delay)

        yield self.dest, GeneStateInfo(level, is_active)


class Gene:
    def __init__(self, *, name=None):
        global _genes

        assert name, "Gene must have a name"
        self.name = name

        _genes.append(self)
        self._set_ligand = None
        self._is_ligand = None

    def __repr__(self):
        return f"Gene('{self.name}')"

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
        if self._is_ligand:
            return self.present() and self.ligand_present()
        else:
            return self.present()

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

    def is_present(self, *, where=None, start=None, duration=None, level=100,
                   decay=1):
        assert where
        assert start
        ix = Interaction_IsPresent(dest=self, start=start, duration=duration,
                                   tissue=where, level=level, decay=decay)
        _add_rule(ix)

    def custom_activation(self, *, state_fn=None, delay=1):
        ix = Interaction_Arbitrary(dest=self, state_fn=state_fn, delay=delay)
        _add_rule(ix)

    def custom_activation2(self, *, state_fn=None, delay=1):
        ix = Interaction_ArbitraryComplex(dest=self, state_fn=state_fn, delay=delay)
        _add_rule(ix)


class Ligand(Gene):
    def __init__(self, *, name=None):
        super().__init__(name=name)
        self._is_ligand = True

    def __repr__(self):
        return f"Ligand('{self.name}')"


class Receptor(Gene):
    def __init__(self, *, name=None, ligand=None):
        super().__init__(name=name)
        assert name
        self._set_ligand = ligand
        if ligand:
            ligand._is_ligand = True

    def __repr__(self):
        return f"Receptor('{self.name}')"

    def ligand(self, *, activator=None, ligand=None):
        # @CTB legacy
        if ligand is None:
            ligand = self._set_ligand
            if ligand is None:
                raise Exception("need to specify a ligand for this receptor, either at creation or here")
        else:
            ligand._is_ligand = True

        ix = Interaction_Activates(source=activator, dest=self, delay=1)
        _add_rule(ix)

    def activated_by(self, *, activator=None, source=None, delay=1):
        if activator is None:   # @CTB deprecated
            activator = source
        if activator is None:
            raise Exception("must supply an activator!")

        ix = Interaction_Activates(source=source, dest=self, delay=delay)
        _add_rule(ix)
