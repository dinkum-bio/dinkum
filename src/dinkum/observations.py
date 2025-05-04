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
        return tissue_state.is_present(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is PRESENT in tissue {self.tissue_name} at time {self.time}"


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
        return not tissue_state.is_present(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is NOT PRESENT in tissue {self.tissue_name} at time {self.time}"


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
        return not tissue_state.is_present(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is NEVER PRESENT in tissue {self.tissue_name}"


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
        return tissue_state.is_present(gene_name=self.gene_name)

    def render(self):
        return f"{self.gene_name} is ALWAYS PRESENT in tissue {self.tissue_name}"


def check_is_always_present(*, gene=None, tissue=None):
    ob = Obs_IsAlwaysPresent(gene=gene, tissue=tissue)
    _add_obs(ob)


class Obs_IsActive(Observation):
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
        return f"{self.gene_name} is ACTIVE in tissue {self.tissue_name} at time {self.time}"


def check_is_active(*, gene=None, time=None, tissue=None):
    ob = Obs_IsActive(gene=gene, time=time, tissue=tissue)
    _add_obs(ob)


class Obs_IsNotActive(Observation):
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
        return f"{self.gene_name} is NOT ACTIVE in tissue {self.tissue_name} at time {self.time}"


def check_is_not_active(*, gene=None, time=None, tissue=None):
    ob = Obs_IsNotActive(gene=gene, time=time, tissue=tissue)
    _add_obs(ob)


class Obs_LevelIsBetween(Observation):
    def __init__(
        self, *, gene=None, time=None, tissue=None, min_level=None, max_level=None
    ):
        assert gene
        assert time is not None
        assert tissue is not None
        assert min_level is not None
        assert max_level is not None
        self.gene_name = gene
        self.time = time
        self.tissue_name = tissue
        self.min_level = min_level
        self.max_level = max_level

    def check(self, state):
        # not applicable
        if state.time != self.time:
            return None

        tissue_state = state.get_by_tissue_name(self.tissue_name)
        level = tissue_state.get_level(self.gene_name)
        if level >= self.min_level and level <= self.max_level:
            return True
        return False

    def render(self):
        return f"{self.gene_name} has level NOT between {self.min_level} and {self.max_level} in tissue {self.tissue_name} at time {self.time}"


def check_level_is_between(
    *, gene=None, time=None, tissue=None, min_level=None, max_level=None
):
    ob = Obs_LevelIsBetween(
        gene=gene, time=time, tissue=tissue, min_level=min_level, max_level=max_level
    )
    _add_obs(ob)


def test_observations(state):
    succeed = True
    for ob in get_obs():
        check = ob.check(state)
        if check is None:
            pass
        elif check:
            print("passed:", ob.render())
        else:
            print("** FAILED: it is not true that:", ob.render())
            succeed = False

    return succeed
