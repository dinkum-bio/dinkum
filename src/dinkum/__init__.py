#
from . import vfg
from . import vfn
from . import observations

import itertools

class DinkumException(Exception):
    pass

class DinkumObservationFailed(DinkumException):
    pass


def reset():
    vfg.reset()
    vfn.reset()
    observations.reset()


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


class Timecourse:
    def iterate(self, *, start=None, stop=None, verbose=False):
        assert start is not None
        assert stop is not None

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

        # advance one tick at a time
        for i in range(start + 1, stop + 1):
            yield this_state
            next_state = State(tissues=tissues, time=i)

            for t in tissues:
                this_active = this_state[t]

                next_active = GeneActivity()

                # bring in maternal/always active
                seen = set()
                for g in t.all_active(at=i):
                    next_active.set_activity(gene=g, active=1)
                    seen.add(g.name)

                for r in vfg.get_rules():
                    # advance state of all genes based on last state
                    for g, activity in r.advance(state=this_state,
                                                 tissue=t):
                        if g.name in seen:
                            raise DinkumException(f"multiple rules containing {g.name}")
                        next_active.set_activity(gene=g, active=activity)
                        seen.add(g.name)


                next_state[t] = next_active

            # advance => next state
            this_state = next_state

        yield this_state


def run(start, stop):
    # run time course
    tc = Timecourse()
    for state in tc.iterate(start=start, stop=stop):
        print(f"time={state.time}")
        for ti in state.tissues:
            present = state[ti]
            print(f"\ttissue={ti.name}, {present.report_activity()}")
        if not observations.test_observations(state):
            raise DinkumObservationFailed(state.time)
