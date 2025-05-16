"""encode view-from-genome rules."""

from functools import total_ordering
import inspect
import collections

from .exceptions import *
from .vfn import check_is_valid_tissue
from .vfg_functions import *


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
            raise DinkumMultipleRules(
                f"multiple rules containing {r.dest} are not allowed!"
            )
        if not r.multiple_allowed:
            seen.add(r.dest)


def get_rules():
    return list(_rules)


def get_gene_names():
    return [g.name for g in sorted(_genes)]


def get_gene(name):
    for g in sorted(_genes):
        if g.name == name:
            return g
    raise DinkumInvalidGene(f"unknown gene name: '{name}'")


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
    # assert isinstance(tissue, Tissue)

    ligands = set()
    for gene in _genes:
        if gene._is_ligand:
            for neighbor in tissue.neighbors:
                if states.is_active(timepoint, delay, gene, neighbor):
                    is_juxtacrine = getattr(gene, "is_juxtacrine", False)

                    # if it's juxtacrine, it only activates other, not self.
                    if is_juxtacrine:
                        if neighbor != tissue:
                            ligands.add(gene)
                    else:
                        ligands.add(gene)

    return ligands


def check_ligand(*, dest, timepoint, states, tissue, delay):
    """If this is a receptor, check that ligand is in neighboring tissue,
    and return True if it is, False else.

    Otherwise, return True.
    """
    if dest.is_receptor:
        ligands_in_neighbors = _retrieve_ligands(timepoint, states, tissue, delay)
        if dest._set_ligand in ligands_in_neighbors:
            return True
        return False
    else:
        return True  # by default, not ligand => is active


class CustomActivation:
    def __init__(self, *, input_genes=None):
        # only support explicit kwargs on __call__,
        # because otherwise we are a bit
        # fragile with respect to order of arguments.
        argspec = inspect.getfullargspec(self.__call__)
        argspec.args.remove("self")
        if argspec.args or argspec.varargs or argspec.varkw or argspec.kwonlydefaults:
            raise DinkumInvalidActivationFunction(
                "on __call__, must supply _only_ kwargs with no defaults"
            )

        if input_genes is None:  # retrieve from __call__
            input_genes = argspec.kwonlyargs

        self.input_genes = list(input_genes)

    def __call__(self, *args, **kw):
        raise NotImplementedError


class Interactions:
    multiple_allowed = False

    def btp_autonomus_links(self):
        raise Unimplemented

    def btp_signal_links(self):
        raise Unimplemented

    def check_ligand(self, timepoint, states, tissue, delay):
        if getattr(self.dest, "_set_ligand", None):
            ligands_in_neighbors = _retrieve_ligands(timepoint, states, tissue, delay)
            if self.dest._set_ligand in ligands_in_neighbors:
                return True
            return False
        else:
            return True  # by default, not ligand => is active


class Interaction_IsPresent(Interactions):
    # @CTB recode as custom2?
    multiple_allowed = True

    def __init__(
        self,
        *,
        dest=None,
        start=None,
        duration=None,
        tissue=None,
        level=None,
        decay=None,
    ):
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
                if (
                    self.duration is None or timepoint < self.start + self.duration
                ):  # active!
                    if self.check_ligand(timepoint, states, tissue, delay=1):
                        yield self.dest, GeneStateInfo(level=self.level, active=True)
                    else:
                        yield self.dest, GeneStateInfo(level=self.level, active=False)
        # we have no opinion on activity outside our tissue!

        self.level = round(self.level / self.decay + 0.5)


class Interaction_Custom(Interactions):
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
            dep_gene_names = [
                g.name if isinstance(g, Gene) else g for g in dep_gene_names
            ]
        else:
            # only support explicit kwargs, because otherwise we are a bit
            # fragile with respect to order of arguments.
            argspec = inspect.getfullargspec(state_fn)
            if (
                argspec.args
                or argspec.varargs
                or argspec.varkw
                or argspec.kwonlydefaults
            ):
                raise DinkumInvalidActivationFunction("must supply kwargs only")
            dep_gene_names = argspec.kwonlyargs

        return dep_gene_names

    def _get_genes_for_activation_fn(self, state_fn):
        # get the activity of upstream genes
        dep_genes = []
        dep_gene_names = self._get_gene_names(state_fn)
        for name in dep_gene_names:
            found = False
            g = get_gene(name)
            check_is_tf(g)
            dep_genes.append((name, g))

        return dep_genes

    def advance(self, *, timepoint=None, states=None, tissue=None):
        # 'states' is class States...
        if not states:
            #            assert 0            # what is this if for??
            return

        assert tissue
        dep_genes = self._get_genes_for_activation_fn(self.state_fn)

        # pass in their full GeneStateInfo
        delay = self.delay

        dep_state = {}
        for name, gene in dep_genes:
            gene_state = states.get_gene_state_info(
                timepoint=timepoint, delay=delay, gene=gene, tissue=tissue
            )
            if gene_state is None:
                gene_state = DEFAULT_OFF

            dep_state[name] = gene_state

        result = self.state_fn(**dep_state)
        if result is not None:
            if not isinstance(result, GeneStateInfo):
                if len(tuple(result)) == 2:
                    result = GeneStateInfo(int(result[0]), bool(result[1]))

            if not isinstance(result, GeneStateInfo):
                raise DinkumInvalidActivationResult(
                    f"result '{result}' of custom activation function '{self.state_fn.__name__}' is not a GeneStateInfo tuple (and cannot be converted)"
                )

        if result is not None:
            level, is_active = result

            if is_active:
                is_active = self.check_ligand(timepoint, states, tissue, self.delay)

            yield self.dest, GeneStateInfo(level, is_active)


class Interaction_CustomObj(Interactions):
    """
    An interaction that supports even more powerful arbitrary logic & levels.
    """

    def __init__(self, *, dest=None, obj=None):
        assert dest is not None
        assert obj is not None
        self.dest = dest
        self.obj = obj

    def btp_autonomous_links(self):
        return []

    def btp_signal_links(self):
        return []

    def advance(self, *, timepoint=None, states=None, tissue=None):
        assert tissue

        result = self.obj.advance(timepoint, states, tissue)
        if result is not None:
            target, gsi = result
            if not isinstance(gsi, GeneStateInfo):
                raise DinkumInvalidActivationResult(
                    f"result '{result}' of custom2 activation function '{self.obj} is not a GeneStateInfo tuple"
                )
            yield target, gsi


class Gene:
    is_receptor = False
    is_ligand = False
    is_tf = True

    def __init__(self, *, name=None):
        global _genes

        assert name, "Gene must have a name"
        self.name = name

        _genes.append(self)
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

    def active(self):  # present = active
        if self._is_ligand:
            return self.present() and self._set_ligand
        else:
            return self.present()

    def activated_by(self, *, source=None, delay=1):
        check_is_valid_gene(self)
        check_is_tf(source)
        self.custom_obj(Activator(rate=100, activator_name=source.name, delay=delay))

    def activated_by_or(self, *, sources=None, delay=1):
        for src in sources:
            check_is_tf(src)
        check_is_valid_gene(self)
        names = [src.name for src in sources]
        weights = [1] * len(names)  # OR
        self.custom_obj(
            LogisticMultiActivator(
                activator_names=names, weights=weights, delay=delay, rate=100
            )
        )

    activated_or = activated_by_or

    def and_not(self, *, activator=None, repressor=None, delay=1):
        check_is_valid_gene(self)
        check_is_tf(activator)
        check_is_tf(repressor)
        self.custom_obj(
            Repressor(
                activator_name=activator.name,
                repressor_name=repressor.name,
                delay=delay,
                repressor_rate=100,
                activator_rate=100,
            )
        )

    def activated_by_and(self, *, sources, delay=1):
        for src in sources:
            check_is_tf(src)
        check_is_valid_gene(self)
        names = [src.name for src in sources]
        weights = [1 / len(names)] * len(names)  # AND
        self.custom_obj(
            LogisticMultiActivator(
                activator_names=names,
                weights=weights,
                delay=delay,
                rate=100,
                midpoint=99,
            )
        )

    def is_present(self, *, where=None, start=None, duration=None, level=100, decay=1):
        assert where
        assert start
        check_is_valid_gene(self)
        check_is_valid_tissue(where)
        ix = Interaction_IsPresent(
            dest=self,
            start=start,
            duration=duration,
            tissue=where,
            level=level,
            decay=decay,
        )
        _add_rule(ix)

    def custom_fn(self, *, state_fn=None, delay=1):
        ix = Interaction_Custom(dest=self, state_fn=state_fn, delay=delay)
        _add_rule(ix)

    def custom_obj(self, obj):
        obj.set_gene(self)
        ix = Interaction_CustomObj(dest=self, obj=obj)
        _add_rule(ix)


class Ligand(Gene):
    is_tf = False
    is_ligand = True

    def __init__(self, *, name=None, is_juxtacrine=False):
        super().__init__(name=name)
        self.is_juxtacrine = is_juxtacrine
        self._is_ligand = True

    def __repr__(self):
        return f"Ligand('{self.name}')"


class Receptor(Gene):
    is_receptor = True

    def __init__(self, *, name=None, ligand=None):
        super().__init__(name=name)
        assert name
        self._set_ligand = ligand
        if ligand and not isinstance(ligand, Ligand):
            raise DinkumNotALigand(f"gene {ligand.name} is not a Ligand")

    def __repr__(self):
        return f"Receptor('{self.name}')"
