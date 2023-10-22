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
    return list(obs)


class Observation:
    pass

class Obs_IsPresent(Observation):
    def __init__(self, *, gene=None, time=None, tissue=None):
        self.gene = gene
        self.time = time
        self.tissue = None

    def check(self, state):
        # not applicable
        if state.time != self.time:
            return None

        # applicable! check.
        tissue_state = state[self.tissue]
        if gene in tissue_state.genes:
            return True
        else:
            return False
