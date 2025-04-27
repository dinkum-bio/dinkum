import pytest

import dinkum
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations
from dinkum.exceptions import *


def test_1():
    dinkum.reset()

    x = Gene(name='X')
    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5

    states = tc.get_states()
    df = states.to_dataframe()

    # @CTB test something here!
