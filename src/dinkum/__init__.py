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


class State:
    def __init__(self, *, tissues=None, time=None):
        assert tissues is not None
        assert time is not None

        self._tissues = list(tissues)
        self._tissues_by_name = {}
        self.time = time

    def __setitem__(self, tissue, genes):
        assert tissue in self._tissues
        self._tissues_by_name[tissue.name] = set(genes)

    def __getitem__(self, tissue):
        return self._tissues_by_name[tissue.name]

    def get_by_tissue_name(self, tissue_name):
        return self._tissues_by_name[tissue_name]

    @property
    def tissues(self):
        return self._tissues


class Timecourse:
    def iterate(self, *, start=None, stop=None):
        assert start is not None
        assert stop is not None

        tissues = vfn.get_tissues()
        print(f"got {len(tissues)} tissues.")
        for t in tissues:
            print(f"\ttissue {t.name}")

        # initialize state at start
        this_state = State(tissues=tissues, time=start)
        for t in tissues:
            this_present = set()
            this_present.update(t.all_present(at=start))

            this_state[t] = this_present

        for i in range(start + 1, stop + 1):
            yield this_state
            next_state = State(tissues=tissues, time=i)

            for t in tissues:
                this_present = this_state[t]

                next_present = set()
                for r in vfg.get_rules():
                    for g in r.advance(present=this_present):
                        next_present.add(g)
                    next_present.update(t.all_present(at=i))

                next_state[t] = next_present

            this_state = next_state

def run(start, stop):
    # run time course
    tc = Timecourse()
    for state in tc.iterate(start=start, stop=stop):
        print(f"time={state.time}")
        for ti in state.tissues:
            present = state[ti]
            print(f"\ttissue={ti.name}, {[ g.name for g in present ]}")
        if not observations.test_observations(state):
            raise DinkumObservationFailed(state.time)
