"""encode observations

X is-transcribed in celltype/tissue/compartment M at time T
Z is-present in M at time T
X is-up in M at time T
Y is-down in M at time T
"""

_obs = []
def _add_obs(ob):
    assert isinstance(ob, Observation), f"{ob} must be an Observation"
    global _obs
    _obs.append(ob)

def get_obs():
    return list(_obs)

def reset():
    global _obs
    _obs = []

class Observation:
    pass

class Obs_IsPresent(Observation):
    def __init__(self, *, gene=None, time=None, tissue=None):
        assert gene, "gene must be set"
        assert time is not None, "time must be set"
        assert tissue is not None, "tissue must be set"
        self.gene_name = gene
        self.time = time
        self.tissue_name = tissue

    def check(self, state):
        # not applicable
        if state.time != self.time:
            return None

        tissue_state = state.get_by_tissue_name(self.tissue_name)
        return tissue_state.is_active(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is ON in tissue {self.tissue_name} at time {self.time}"

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
        return not tissue_state.is_active(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is NOT ON in tissue {self.tissue_name} at time {self.time}"

def check_is_not_present(*, gene=None, time=None, tissue=None):
    ob = Obs_IsNotPresent(gene=gene, time=time, tissue=tissue)
    _add_obs(ob)



class Obs_IsNeverPresent(Observation):
    def __init__(self, *, gene=None, tissue=None):
        assert gene
        assert tissue is not None
        self.gene_name = gene
        self.tissue_name = tissue

    def check(self, state):
        tissue_state = state.get_by_tissue_name(self.tissue_name)
        return not tissue_state.is_active(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is NEVER ON in tissue {self.tissue_name}"

def check_is_never_present(*, gene=None, tissue=None):
    ob = Obs_IsNeverPresent(gene=gene, tissue=tissue)
    _add_obs(ob)


class Obs_IsAlwaysPresent(Observation):
    def __init__(self, *, gene=None, tissue=None):
        assert gene
        assert tissue is not None
        self.gene_name = gene
        self.tissue_name = tissue

    def check(self, state):
        tissue_state = state.get_by_tissue_name(self.tissue_name)
        return tissue_state.is_active(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is ALWAYS ON in tissue {self.tissue_name}"


def check_is_always_present(*, gene=None, tissue=None):
    ob = Obs_IsAlwaysPresent(gene=gene, tissue=tissue)
    _add_obs(ob)


def test_observations(state):
    succeed = True
    for ob in get_obs():
        check = ob.check(state)
        if check is None:
            pass
        elif check:
            print('passed:', ob.render())
        else:
            print('** FAILED:', ob.render())
            succeed = False

    return succeed
