import pytest

import dinkum
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations


def test_maternal():
    dinkum.reset()

    x = Gene(name='X')
    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    # set observations
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='X', time=2, tissue='M')

    # run time course
    dinkum.run(1, 5)


def test_maternal_fail():
    dinkum.reset()

    x = Gene(name='X')
    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    # set observations
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_present(gene='X', time=2, tissue='M')

    with pytest.raises(dinkum.DinkumObservationFailed):
        dinkum.run(1, 5)


def test_activation():
    dinkum.reset()

    # set it all up!
    x = Gene(name='X')
    y = Gene(name='Y')

    y.activated_by(source=x)
    x.activated_by(source=y)

    m = Tissue(name='M')
    x.is_present(where=m, start=1)

    # set observations
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='Y', time=1, tissue='M')
    observations.check_is_present(gene='X', time=2, tissue='M')
    observations.check_is_present(gene='Y', time=2, tissue='M')

    # run!
    dinkum.run(start=1, stop=5)


def test_simple_repression():
    # establish preconditions
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='Y', time=1, tissue='M')

    observations.check_is_present(gene='X', time=2, tissue='M')
    observations.check_is_present(gene='Y', time=2, tissue='M')
    observations.check_is_present(gene='Z', time=2, tissue='M')

    observations.check_is_not_present(gene='Z', time=3, tissue='M')

    # set it all up!
    x = Gene(name='X')
    y = Gene(name='Y')
    z = Gene(name='Z')

    y.activated_by(source=x)
    x.activated_by(source=y)

    z.and_not(activator=x, repressor=y)

    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    # run!
    dinkum.run(1, 5)


def test_simple_multiple_tissues():
    ## tissue M
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='Y', time=1, tissue='M')

    observations.check_is_present(gene='X', time=2, tissue='M')
    observations.check_is_present(gene='Y', time=2, tissue='M')
    observations.check_is_present(gene='Z', time=2, tissue='M')

    observations.check_is_not_present(gene='Z', time=3, tissue='M')

    ## tissue N
    observations.check_is_present(gene='Y', time=1, tissue='N')
    observations.check_is_not_present(gene='X', time=1, tissue='N')
    observations.check_is_present(gene='X', time=2, tissue='N')
    observations.check_is_present(gene='Y', time=2, tissue='N')

    # set it all up!

    ## VFG
    x = Gene(name='X')
    y = Gene(name='Y')
    z = Gene(name='Z')

    y.activated_by(source=x)
    x.activated_by(source=y)

    z.and_not(activator=x, repressor=y)

    ## VFN
    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    n = Tissue(name='N')
    y.is_present(where=n, start=1, duration=1)

    # run!
    dinkum.run(1, 5)
