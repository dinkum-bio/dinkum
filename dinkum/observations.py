"""encode observations

X is-transcribed in celltype/tissue/compartment M at time T
Z is-present in M at time T
X is-up in M at time T
Y is-down in M at time T
"""

_obs = []
def _add_obs(ob):
    assert isinstance(ob, Observation)
    global _obs
    _obs.append(ob)

def get_obs():
    return list(_obs)


class Observation:
    pass

class Obs_IsPresent(Observation):
    def __init__(self, *, gene=None, time=None, tissue=None):
        assert gene
        assert time is not None
        assert tissue is not None
        self.gene_name = gene
        self.time = time
        self.tissue_name = tissue

    def check(self, state):
        # not applicable
        if state.time != self.time:
            return None

        tissue_state = state.get_by_tissue_name(self.tissue_name)
        for g in tissue_state:
            if g.name == self.gene_name:
                return True
        return False

def check_is_present(*, gene=None, time=None, tissue=None):
    ob = Obs_IsPresent(gene=gene, time=time, tissue=tissue)
    _add_obs(ob)



class Obs_IsNotPresent(Observation):
    def __init__(self, *, gene=None, time=None, tissue=None):
        assert gene
        assert time is not None
        assert tissue is not None
        self.gene_name = gene
        self.time = time
        self.tissue_name = tissue

    def check(self, state):
        # not applicable
        if state.time != self.time:
            return None

        tissue_state = state.get_by_tissue_name(self.tissue_name)
        for g in tissue_state:
            if g.name == self.gene_name:
                return False
        return True

def check_is_not_present(*, gene=None, time=None, tissue=None):
    ob = Obs_IsNotPresent(gene=gene, time=time, tissue=tissue)
    _add_obs(ob)


def test_observations(state):
    for ob in get_obs():
        check = ob.check(state)
        if check is None:
            pass
        elif check:
            print('passed', ob)
        else:
            print('failed', ob)
