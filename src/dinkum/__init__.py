import sys
from . import vfg
from . import vfn
from . import observations

import itertools
import collections

class DinkumException(Exception):
    pass

class DinkumObservationFailed(DinkumException):
    pass


def reset():
    vfg.reset()
    vfn.reset()
    observations.reset()


def run_and_display(*, start=1, stop=10, gene_names=None, tissue_names=None,
            verbose=False):
    "@CTB document."
    from dinkum.display import MultiTissuePanel, tc_record_activity
    if not gene_names:
        gene_names = vfg.get_gene_names()
    if not tissue_names:
        tissue_names = vfn.get_tissue_names()

    try:
        states, tissues, is_active_fn = tc_record_activity(start=start,
                                                           stop=stop,
                                                           gene_names=gene_names,
                                                           verbose=verbose)
    except DinkumException as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        print("Halting execution.", file=sys.stderr)
        return

    #@CTB what is tissues here and how might it differ from tissue_names?
    # perhaps just set it...
    #print(tissues)

    mp = MultiTissuePanel(states=states, tissue_names=tissue_names,
                          genes_by_name=gene_names)
    return mp.draw(is_active_fn)


class GeneActivity:
    def __init__(self):
        self.genes_by_name = {}

    def set_activity(self, *, gene=None, active=None):
        assert gene
        assert active is not None
        self.genes_by_name[gene.name] = active

    def is_active(self, gene_name):
        return self.genes_by_name.get(gene_name, False)

    def __contains__(self, gene):
        return self.is_active(gene.name)

    def report_activity(self):
        rl = []
        for k, v in sorted(self.genes_by_name.items()):
            if v:
                v = 1
            else:
                v = 0
            rl.append(f"{k}={v}")

        return rl


class State:
    """
    Hold the gene activity state for multiple tissues.

    Holds multiple tissue, each with their own GeneActivity object.

    dict interface supports getting and setting gene activity (value) by
    tissue (key).

    `tissues` attribute provides enumeration of available tissues

    `is_active(gene, tissue)` returns True/False around activity.
    """
    def __init__(self, *, tissues=None, time=None):
        assert tissues is not None
        assert time is not None

        self._tissues = list(tissues)
        self._tissues_by_name = {}
        self.time = time

    def __setitem__(self, tissue, genes):
        assert tissue in self._tissues
        assert isinstance(genes, GeneActivity)
        self._tissues_by_name[tissue.name] = genes

    def __getitem__(self, tissue):
        return self._tissues_by_name[tissue.name]

    def get_by_tissue_name(self, tissue_name):
        return self._tissues_by_name[tissue_name]

    @property
    def tissues(self):
        return self._tissues

    def is_active(self, gene, tissue):
        ts = self[tissue]
        if ts and gene in ts:
            return True
        return False


class States(collections.UserDict):
    """
    Contains (potentially incomplete) set of tissue/gene states for many
    timepoints.
    """
    def __init__(self):
        self.data = {}

    def is_active(self, current_tp, delay, gene, tissue):
        from .vfg import Gene

        assert int(current_tp)
        assert int(delay)
        assert isinstance(gene, Gene)

        check_tp = current_tp - delay
        state = self.get(check_tp)
        if state and state.is_active(gene, tissue):
            return True
        return False


class Timecourse:
    """
    Run a time course for a system b/t two time points, start and stop.
    """
    def __init__(self, *, start=None, stop=None):
        assert start is not None
        assert stop is not None

        print(f'start={start} stop={stop}')

        self.start = start
        self.stop = stop
        self.states_d = States()

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
            print('')

        # initialize state at start for all tissues
        this_state = State(tissues=tissues, time=start)
        for t in tissues:
            this_active = GeneActivity()
            for g in t.all_active(at=start):
                this_active.set_activity(gene=g, active=True)

            this_state[t] = this_active

        self.states_d[start] = this_state

        # advance one tick at a time
        for tp in range(start + 1, stop + 1):
            next_state = State(tissues=tissues, time=tp)

            for t in tissues:
                this_active = this_state[t]

                next_active = GeneActivity()

                # bring in maternal/always active
                seen = set()
                for g in t.all_active(at=tp):
                    next_active.set_activity(gene=g, active=1)
                    seen.add(g.name)

                for r in vfg.get_rules():
                    # advance state of all genes based on last state
                    for g, activity in r.advance(timepoint=tp,
                                                 states=self.states_d,
                                                 tissue=t):
                        if g.name in seen:
                            raise DinkumException(f"multiple rules containing {g.name}")
                        next_active.set_activity(gene=g, active=activity)
                        seen.add(g.name)


                next_state[t] = next_active

            # advance => next state
            self.states_d[tp] = next_state
            this_state = next_state

    def check(self):
        "Test all of the observations for all of the states."
        for state in iter(self):
            if not observations.test_observations(state):
                raise DinkumObservationFailed(state.time)


def run(start, stop, *, verbose=False):
    # run time course
    tc = Timecourse(start=start, stop=stop)
    tc.run(verbose=verbose)

    for state in tc:
        print(f"time={state.time}")
        for ti in state.tissues:
            present = state[ti]
            print(f"\ttissue={ti.name}, {present.report_activity()}")
        if not observations.test_observations(state):
            raise DinkumObservationFailed(state.time)

    return tc
