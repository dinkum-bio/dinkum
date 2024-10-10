import pytest

import dinkum
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations


def test_custom_1():
    dinkum.reset()

    x = Gene(name='X')
    y = Gene(name='Y')
    m = Tissue(name='M')

    def activator_fn(X):
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


def test_custom_1_fail():
    dinkum.reset()

    x = Gene(name='X')
    y = Gene(name='Y')
    m = Tissue(name='M')

    def activator_fn(Z):
        return X

    x.is_present(where=m, start=1, duration=1)
    y.custom_activation(state_fn=activator_fn, delay=1)

    # run time course
    with pytest.raises(Exception): # @CTB change to Dinkum exception
        tc = dinkum.run(1, 5)
