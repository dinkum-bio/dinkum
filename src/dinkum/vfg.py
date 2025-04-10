"""encode view-from-genome rules.
"""
from functools import total_ordering
import inspect
import collections

from .exceptions import *
from .vfn import check_is_valid_tissue

class GeneStateInfo:
    def __init__(self, level=0, active=False):
        self.level = level
        self.active = active

    def __bool__(self):
        return bool(self.level > 0 and self.active)

    def __iter__(self):
        return iter([self.level, self.active])

    def __repr__(self):
        return f"<level={self.level},active={self.active}>"

DEFAULT_OFF = GeneStateInfo(level=0, active=False)

_rules = []
_genes = []

def _add_rule(ix):
    _rules.append(ix)

    # check!
    seen = set()
    for r in _rules:
        if r.dest in seen and not r.multiple_allowed:
            raise DinkumException(f"multiple rules containing {r.dest} are not allowed!")
        if not r.multiple_allowed:
            seen.add(r.dest)


def get_rules():
    return list(_rules)

def get_gene_names():
    return [ g.name for g in sorted(_genes) ]

def get_gene(name):
    for g in sorted(_genes):
        if g.name == name:
            return g
    raise Exception(f"unknown genome name: '{name}'")

def reset():
    global _rules
    global _genes
    _rules = []
    _genes = []


def check_is_valid_gene(g):
    if not g in _genes:
        raise DinkumInvalidGene(f"{g.name} is invalid")

def check_is_tf(g):
    check_is_valid_gene(g)
    if not g.is_tf:
        raise DinkumNotATranscriptionFactor(f"{g.name} is not a transcription factor")


def _retrieve_ligands(timepoint, states, tissue, delay):
    "Retrieve all ligands in neighboring tissues for the given timepoint/delay"
    #assert isinstance(tissue, Tissue)

    ligands = set()
    for gene in _genes:
        if gene._is_ligand:
            for neighbor in tissue.neighbors:
                # for neighbor in (tissue, *tissue.neighbors):@CTB
                if states.is_active(timepoint, delay, gene, neighbor):
                    ligands.add(gene)

    return ligands


class CustomActivation(object):
    def __init__(self, *, input_genes=None):
        # only support explicit kwargs on __call__,
        # because otherwise we are a bit
        # fragile with respect to order of arguments.
        argspec = inspect.getfullargspec(self.__call__)
        argspec.args.remove('self')
        if argspec.args or argspec.varargs or argspec.varkw or \
           argspec.kwonlydefaults:
            raise DinkumInvalidActivationFunction("on __call__, must supply _only_ kwargs with no defaults")

        if input_genes is None: # retrieve from __call__
            input_genes = argspec.kwonlyargs

        self.input_genes = list(input_genes)


class Interactions:
    multiple_allowed = False

    def btp_autonomus_links(self):
        raise Unimplemented

    def btp_signal_links(self):
        raise Unimplemented

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

    def btp_autonomous_links(self):
        return []

    def btp_signal_links(self):
        return []

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

        self.level = round(self.level / self.decay + 0.5)


class Interaction_Activates(Interactions):
    def __init__(self, *, source=None, dest=None, delay=1):
        check_is_valid_gene(dest)
        check_is_tf(source)

        self.src = source
        self.dest = dest
        self.delay = delay

    def btp_autonomous_links(self):
        yield self.dest, self.src, "positive"

    def btp_signal_links(self):
        return []

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


class Interaction_Or(Interactions):
    def __init__(self, *, sources=None, dest=None, delay=1):
        check_is_valid_gene(dest)
        for g in sources:
            check_is_tf(g)

        for tf in sources:
            if not tf.is_tf:
                raise DinkumNotATranscriptionFactor(tf.name)
        self.sources = sources
        self.dest = dest
        self.delay = delay

    def btp_autonomous_links(self):
        for src in self.sources:
            return [self.dest, self.src, "positive"]

    def btp_signal_links(self):
        return []

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


class Interaction_AndNot(Interactions):
    def __init__(self, *, source=None, repressor=None, dest=None, delay=1):
        check_is_valid_gene(dest)
        check_is_tf(source)
        check_is_tf(repressor)

        self.src = source
        self.repressor = repressor
        self.dest = dest
        self.delay = delay
        for tf in (self.src, self.repressor):
            if not tf.is_tf:
                raise DinkumNotATranscriptionFactor(tf.name)

    def btp_autonomous_links(self):
        yield [self.dest, self.src, "positive"]
        yield [self.dest, self.repressor, "negative"]

    def btp_signal_links(self):
        return []

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


class Interaction_And(Interactions):
    def __init__(self, *, sources=None, dest=None, delay=1):
        check_is_valid_gene(dest)
        for g in sources:
            check_is_tf(g)

        self.sources = sources
        self.dest = dest
        self.delay = delay
        for tf in self.sources:
            if not tf.is_tf:
                raise DinkumNotATranscriptionFactor(tf.name)

    def btp_autonomous_links(self):
        for src in self.sources:
            yield [self.dest, src, "positive"]

    def btp_signal_links(self):
        return []

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


class Interaction_ToggleRepressed(Interactions):
    def __init__(self, *, tf=None, cofactor=None, dest=None, delay=1):
        check_is_valid_gene(dest)
        check_is_tf(tf)
        check_is_valid_gene(cofactor)

        self.tf = tf
        if not self.tf.is_tf:
            raise DinkumNotATranscriptionFactor(self.tf.name)
        self.cofactor = cofactor
        self.dest = dest
        self.delay = delay

    def btp_signal_links(self):
        return []

    def btp_autonomous_links(self):
        yield [self.dest, self.tf, "positive"]
        yield [self.dest, self.cofactor, "positive"] # @CTB toggle

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


class Interaction_Arbitrary(Interactions):
    """
    An interaction that supports arbitrary logic, + levels.
    """
    def __init__(self, *, dest=None, state_fn=None, delay=1):
        assert dest
        assert state_fn
        _check_gene_names = self._get_gene_names(state_fn)

        self.dest = dest
        self.state_fn = state_fn
        self.delay = delay


    def btp_autonomous_links(self):
        return []

    def btp_signal_links(self):
        return []

    def _get_gene_names(self, state_fn):
        # get the names of the genes on the function
        if isinstance(state_fn, CustomActivation):
            dep_gene_names = state_fn.input_genes

            # allow genes, or gene names
            dep_gene_names = [ g.name if isinstance(g, Gene) else g
                               for g in dep_gene_names ]
        else:
            # only support explicit kwargs, because otherwise we are a bit
            # fragile with respect to order of arguments.
            argspec = inspect.getfullargspec(state_fn)
            if argspec.args or argspec.varargs or argspec.varkw or \
               argspec.kwonlydefaults:
                raise DinkumInvalidActivationFunction("must supply kwargs only")
            dep_gene_names = argspec.kwonlyargs

        return dep_gene_names

    def _get_genes_for_activation_fn(self, state_fn):
        # get the activity of upstream genes
        dep_genes = []
        dep_gene_names = self._get_gene_names(state_fn)
        for name in dep_gene_names:
            found = False
            for g in _genes:
                if g.name == name:
                    check_is_tf(g)
                    dep_genes.append((name, g))
                    found = True
                    break
            if not found:
                raise DinkumInvalidGene(name)

        return dep_genes

    def advance(self, *, timepoint=None, states=None, tissue=None):
        # 'states' is class States...
        if not states:
            return

        assert tissue
        dep_genes = self._get_genes_for_activation_fn(self.state_fn)

        # pass in their full GeneStateInfo
        delay = self.delay
        dep_state = { name: states.get_gene_state_info(timepoint, delay, g, tissue)
                      for (name, g) in dep_genes }

        result = self.state_fn(**dep_state)
        if result is not None:
            if not isinstance(result, GeneStateInfo):
                try:
                    if len(tuple(result)) == 2:
                        result = GeneStateInfo(int(result[0]), bool(result[1]))
                except:
                    pass
            if not isinstance(result, GeneStateInfo):
                raise DinkumInvalidActivationResult(f"result '{result}' of custom activation function is not a GeneStateInfo tuple (and cannot be converted)")

        if result is not None:
            level, is_active = result

            if is_active:
                # @CTB ...check ligand should be made more general...
                is_active = self.check_ligand(timepoint,
                                              states,
                                              tissue,
                                              self.delay)

            yield self.dest, GeneStateInfo(level, is_active)


class Gene:
    is_receptor = False
    is_tf = True

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
        check_is_valid_gene(self)
        check_is_valid_gene(source)
        ix = Interaction_Activates(source=source, dest=self, delay=delay)
        _add_rule(ix)

    def activated_by_or(self, *, sources=None, delay=1):
        ix = Interaction_Or(sources=sources, dest=self, delay=delay)
        _add_rule(ix)
    activated_or = activated_by_or

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
        check_is_valid_gene(self)
        check_is_valid_tissue(where)
        ix = Interaction_IsPresent(dest=self, start=start, duration=duration,
                                   tissue=where, level=level, decay=decay)
        _add_rule(ix)

    def custom_activation(self, *, state_fn=None, delay=1):
        ix = Interaction_Arbitrary(dest=self, state_fn=state_fn, delay=delay)
        _add_rule(ix)


class Ligand(Gene):
    is_tf = False

    def __init__(self, *, name=None):
        super().__init__(name=name)
        self._is_ligand = True

    def __repr__(self):
        return f"Ligand('{self.name}')"


class Receptor(Gene):
    is_receptor = True

    def __init__(self, *, name=None, ligand=None):
        super().__init__(name=name)
        assert name
        self._set_ligand = ligand
        if ligand:
            ligand._is_ligand = True

    def __repr__(self):
        return f"Receptor('{self.name}')"

    def ligand(self, *, activator=None, ligand=None):
        # assert 0, "legacy method!"
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
            # assert 0
            activator = source
        if activator is None:
            raise Exception("must supply an activator!")

        ix = Interaction_Activates(source=source, dest=self, delay=delay)
        _add_rule(ix)
