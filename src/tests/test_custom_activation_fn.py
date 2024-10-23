import pytest

import dinkum
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations

from dinkum.exceptions import *
from dinkum.log import *


def test_custom_1():
    # basic does-it-work
    dinkum.reset()

    x = Gene(name='X')
    y = Gene(name='Y')
    m = Tissue(name='M')

    def activator_fn(*, X):
        return X

    x.is_present(where=m, start=1, duration=1)
    y.custom_activation(state_fn=activator_fn, delay=1)

    # set observations
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='X', time=2, tissue='M')
    observations.check_is_not_present(gene='Y', time=1, tissue='M')
    observations.check_is_present(gene='Y', time=2, tissue='M')

    observations.check_is_not_present(gene='X', time=3, tissue='M')
    observations.check_is_not_present(gene='Y', time=3, tissue='M')

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_2():
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name='X')
    y = Gene(name='Y')
    z = Gene(name='Z')
    m = Tissue(name='M')

    def activator_fn(*, Z, X):     # allow order independence
        return X

    x.is_present(where=m, start=1, duration=1)
    z.is_present(where=m, start=1, duration=1)
    y.custom_activation(state_fn=activator_fn, delay=1)

    # set observations
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='X', time=2, tissue='M')
    observations.check_is_not_present(gene='Y', time=1, tissue='M')
    observations.check_is_present(gene='Y', time=2, tissue='M')

    observations.check_is_not_present(gene='X', time=3, tissue='M')
    observations.check_is_not_present(gene='Y', time=3, tissue='M')

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_1_fail():
    # no gene Z
    dinkum.reset()

    x = Gene(name='X')
    y = Gene(name='Y')
    m = Tissue(name='M')

    def activator_fn(*, Z):
        return X

    x.is_present(where=m, start=1, duration=1)
    y.custom_activation(state_fn=activator_fn, delay=1)

    # run time course
    with pytest.raises(DinkumInvalidGene):
        tc = dinkum.run(1, 5)


def test_custom_bad_defn():
    # invalid activation function defns - must be kwargs
    dinkum.reset()

    x = Gene(name='X')
    y = Gene(name='Y')
    m = Tissue(name='M')

    # must be explicitly named, not positional
    def activator_fn(Z):
        return X
    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_activation(state_fn=activator_fn, delay=1)

    # no defaults
    def activator_fn(Z=None):
        return X
    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_activation(state_fn=activator_fn, delay=1)

    # no defaults 2
    def activator_fn(*, Z=None):
        return X
    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_activation(state_fn=activator_fn, delay=1)

    # no general kwargs
    def activator_fn(**kwargs):
        return X
    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_activation(state_fn=activator_fn, delay=1)

    # no general args
    def activator_fn(*args):
        return X
    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_activation(state_fn=activator_fn, delay=1)
