"""
Top-level dinkum module.

Contains core execution information.
"""

import sys
from importlib.metadata import version
import pandas as pd

__version__ = version("dinkum-bio")

import itertools
import collections

from . import vfg
from .vfg import GeneStateInfo, DEFAULT_OFF, get_gene, Gene
from . import vfn
from .vfn import get_tissue, Tissue
from . import observations
from . import utils
from .exceptions import *


def reset(*, verbose=True):
    vfg.reset()
    vfn.reset()
    observations.reset()
    if verbose:
        print(f"initializing: dinkum v{__version__}")


def run_and_display_df(
    *,
    start=1,
    stop=10,
    gene_names=None,
    tissue_names=None,
    verbose=False,
    save_image=None,
    trace_fn=None,
):
    """
    Run and display the circuit model; for use in Jupyter notebooks.

    Key parameter:
    - start (default: 1)
    - stop (default 10)

    Other parameters:
    - 'tissue_names' - a list of tissue names to display (default: all).
    - 'gene_names' - a list of gene names to display (default: all).
    - 'verbose' - display more text output.
    - 'save_image' - save image to this file.
    - 'canvas_type' - 'ipycanvas' or 'pillow' (default: 'pillow')
    """
    from dinkum.display import MultiTissuePanel

    if not gene_names:
        gene_names = vfg.get_gene_names()
    if not tissue_names:
        tissue_names = vfn.get_tissue_names()

    try:
        tc = _run(start=start, stop=stop, verbose=verbose, trace_fn=trace_fn)
    except DinkumException as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        print("Halting execution.", file=sys.stderr)
        return None

    states = tc.get_states()

    level_df, active_df = states.to_dataframe(gene_names)

    mp = MultiTissuePanel(
        states=states,
        tissue_names=tissue_names,
        gene_names=gene_names,
        save_image=save_image,
    )
    return mp.draw(states), level_df, active_df


def run_and_display(*args, **kwargs):
    display_obj, _level_df, _active_df = run_and_display_df(*args, **kwargs)
    return display_obj


def run_df(*args, **kwargs):
    display_obj, _level_df, _active_df = run_and_display_df(*args, **kwargs)
    return display_obj


class OnlyGeneStates:
    """Gene states.

    A class to contain/manage GeneStateInfo objects for multiple genes,
    at a particular time and in a particular tissue.
    """

    def __init__(self):
        self.genes_by_name = {}

    def __repr__(self):
        return repr(self.genes_by_name)

    def __getitem__(self, gene_name):
        return self.get_gene_state(gene_name)

    def set_gene_state(self, *, gene=None, state_info=None):
        assert gene is not None
        assert state_info is not None
        self.genes_by_name[gene.name] = state_info

    def is_present(self, gene_name):
        state_info = self.genes_by_name.get(gene_name)
        if state_info is None or state_info.level == 0:
            return False

        return True

    def is_active(self, gene_name):
        state_info = self.genes_by_name.get(gene_name, DEFAULT_OFF)
        return state_info.level > 0 and state_info.active

    def get_level(self, gene_name):
        state_info = self.genes_by_name.get(gene_name, DEFAULT_OFF)
        return state_info.level

    def get_gene_state(self, gene_name):
        assert not isinstance(gene_name, vfg.Gene)
        state_info = self.genes_by_name.get(gene_name, DEFAULT_OFF)
        return state_info

    def __contains__(self, gene):
        return self.is_active(gene.name)

    def report_activity(self):
        rl = []
        for k, v in sorted(self.genes_by_name.items()):
            level = v.level
            active = 1 if v.active else 0
            rl.append(f"{k}={level} ({active})")

        return rl


class TissueAndGeneStateAtTime:
    """
    Hold the gene activity state for multiple tissues at a particular tp.

    Holds multiple tissues, each with their own OnlyGeneStates object.

    dict interface supports getting and setting gene activity (value) by
    tissue (key).

    `tissues` attribute provides enumeration of available tissues

    `is_active(gene, tissue)` returns True/False around activity.
    """

    def __init__(self, *, tissues=None, time=None):
        assert tissues is not None
        assert time is not None

        self._tissues = list(tissues)
        self._tissues_by_name = {}  # @CTB do we need to set from tissues?
        self.time = time

    def __setitem__(self, tissue, genes):
        "set Tissue object by name."
        assert tissue in self._tissues
        assert isinstance(genes, OnlyGeneStates)
        self._tissues_by_name[tissue.name] = genes

    def __getitem__(self, tissue):
        "get Tissue object."
        return self._tissues_by_name[tissue.name]

    def get(self, tissue):
        return self._tissues_by_name.get(tissue.name)

    def get_by_tissue_name(self, tissue_name):
        assert not isinstance(tissue_name, vfn.Tissue)
        return self._tissues_by_name[tissue_name]

    @property
    def tissues(self):
        return self._tissues

    def is_active(self, gene, tissue):
        ts = self[tissue]
        if ts and gene in ts:
            return True
        return False

    def get_gene_state_info(self, gene, tissue):
        ts = self[tissue]
        return ts[gene.name]


class TissueGeneStates(collections.UserDict):
    """
    Contains (potentially incomplete) set of tissue/gene states for many
    timepoints.

    The top-level container.

    Indexed by timepoint (integer), returns TissueAndGeneStateAtTime objects.
    """

    def __init__(self):
        self.data = {}

    def is_active(self, current_tp, delay, gene, tissue):
        # @CTB deprecate
        from .vfg import Gene

        assert int(current_tp)
        delay = int(delay)
        assert isinstance(gene, Gene)

        check_tp = current_tp - delay
        state = self.get(check_tp)
        if state and state.is_active(gene, tissue):
            return True
        return False

    def get_gene_state_info(self, *, timepoint, delay=0, gene, tissue):
        from .vfg import Gene

        assert int(timepoint)
        delay = int(delay)
        assert isinstance(gene, Gene)
        assert isinstance(tissue, Tissue)

        check_tp = timepoint - delay
        time_state = self.get(check_tp)
        if time_state:
            return time_state.get_gene_state_info(gene, tissue)
        return None

    def set_gene_state(
        self, *, timepoint=None, tissue_name=None, gene_name=None, state_info=None
    ):
        assert timepoint is not None
        assert tissue_name is not None
        assert gene_name is not None
        assert state_info is not None
        timepoint = int(timepoint)

        tissue = vfn.get_tissue(tissue_name)
        gene = vfg.get_gene(gene_name)

        time_state = self.get(timepoint)
        if time_state is None:
            time_state = TissueAndGeneStateAtTime(tissues=[tissue], time=timepoint)
            self.data[timepoint] = time_state

        gene_state = time_state.get(tissue)
        if gene_state is None:
            gene_state = OnlyGeneStates()
            time_state[tissue] = gene_state

        gene_state.set_gene_state(gene=gene, state_info=state_info)

    def to_dataframe(self, gene_names=None):
        """
        Convert to a pandas DataFrame.
        """
        level_rows = []
        active_rows = []

        # extract gene names?
        all_gene_names = set()
        if gene_names is None:
            for timepoint, state in self.items():
                for tissue in state.tissues:
                    tissue_name = tissue.name
                    activity = state.get_by_tissue_name(tissue_name)
                    all_gene_names.update(activity.genes_by_name)

            gene_names = list(sorted(all_gene_names))

        for timepoint, tissue_and_gene_sat in self.items():
            timepoint_str = f"t={timepoint}"

            # for each tissue, get level of each gene
            for tissue in tissue_and_gene_sat.tissues:
                level_d = dict(
                    tissue=tissue.name, timepoint=timepoint, timepoint_str=timepoint_str
                )
                active_d = dict(
                    tissue=tissue.name, timepoint=timepoint, timepoint_str=timepoint_str
                )
                for gene_name in gene_names:
                    gene = get_gene(gene_name)
                    gsi = self.get_gene_state_info(
                        timepoint=timepoint, delay=0, gene=gene, tissue=tissue
                    )
                    level_d[gene_name] = gsi.level
                    active_d[gene_name] = gsi.active
                level_rows.append(level_d)
                active_rows.append(active_d)

        level_df = pd.DataFrame.from_dict(level_rows).set_index("timepoint")
        active_df = pd.DataFrame.from_dict(active_rows).set_index("timepoint")

        return level_df, active_df


class Timecourse:
    """
    Run and record a time course for a system b/t two time points,
    start and stop.
    """

    def __init__(self, *, start=None, stop=None, trace_fn=None):
        assert start is not None
        assert stop is not None

        print(f"start={start} stop={stop}")

        self.start = start
        self.stop = stop
        self.states_d = TissueGeneStates()
        self.trace_fn = trace_fn

    def reset(self):
        self.states_d = TissueGeneStates()

    def keys(self):
        return self.states_d.keys()

    def __iter__(self):
        return iter(self.states_d.values())

    def __len__(self):
        return len(self.states_d)

    def run(self, *, verbose=False):
        start = self.start
        stop = self.stop

        tissues = vfn.get_tissues()
        if verbose:
            print(f"got {len(tissues)} tissues.")
            for t in tissues:
                print(f"\ttissue {t.name}")
            print("")

        # advance one tick at a time
        this_state = {}
        trace_fn = self.trace_fn
        for tp in range(start, stop + 1):
            next_state = TissueAndGeneStateAtTime(tissues=tissues, time=tp)

            for tissue in tissues:
                seen = set()
                next_active = OnlyGeneStates()
                for r in vfg.get_rules():
                    # advance state of all genes based on last state
                    for gene, state_info in r.advance(
                        timepoint=tp, states=self.states_d, tissue=tissue
                    ):

                        next_active.set_gene_state(gene=gene, state_info=state_info)
                        if trace_fn:
                            trace_fn(
                                tp=tp, tissue=tissue, gene=gene, state_info=state_info
                            )

                next_state[tissue] = next_active
                if verbose:
                    print(tp, tissue.name, next_active)

            # advance => next state
            self.states_d[tp] = next_state
            this_state = next_state

    def check(self):
        "Test all of the observations for all of the states."
        for state in iter(self):
            if not observations.test_observations(state):
                raise DinkumObservationFailed(state.time)

    def get_states(self):
        return self.states_d


def _run(*, start, stop, trace_fn=None, verbose=False):
    "Run a time course. No output by default."
    tc = Timecourse(start=start, stop=stop, trace_fn=trace_fn)
    tc.run(verbose=verbose)
    tc.check()
    return tc


def run(start, stop, *, verbose=False, trace_fn=None):
    """Run a time course in 'headless' mode - minimal output.

    Use for Python script/test execution.
    """
    tc = _run(start=start, stop=stop, verbose=verbose, trace_fn=trace_fn)

    for state in tc:
        print(f"time={state.time}")
        for ti in state.tissues:
            present = state[ti]
            print(f"\ttissue={ti.name}, {present.report_activity()}")
        if not observations.test_observations(state):
            raise DinkumObservationFailed(state.time)

    tc.check()

    return tc
