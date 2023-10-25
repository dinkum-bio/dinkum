import pytest

import dinkum
from dinkum.vfg import Gene, Receptor
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations


def test_signaling():
    dinkum.reset()

    #observations.check_is_present(gene='R', tissue='M', time=2)
    observations.check_is_present(gene='X', tissue='N', time=2)
    observations.check_is_not_present(gene='Y', tissue='M', time=2)
    observations.check_is_present(gene='Y', tissue='M', time=4)

    m = Tissue(name='M')
    n = Tissue(name='N')

    m.add_neighbor(neighbor=n)
    assert n in m.neighbors     # should this be bidirectional? probably.

    x = Gene(name='X')
    a = Gene(name='A')
    r = Receptor(name='R')
    y = Gene(name='Y')

    r.ligand(activator=a, ligand=x)

    y.activated_by(source=r)

    # receptor is always present in M
    m.add_gene(gene=a, start=1)

    # x is present in N at time >= 2
    n.add_gene(gene=x, start=2)

    dinkum.run(1, 5)


def test_community_effect():
    # transient input in one cell => mutual lock on via positive feedback/
    # signalling
    # good question for students: how long does pulse need to be, and why??
    # CTB: _active_ vs _present_
    dinkum.reset()

    observations.check_is_present(gene='A', tissue='M', time=6)
    observations.check_is_not_present(gene='A', tissue='M', time=7)

    observations.check_is_not_present(gene='L', tissue='M', time=1) # not act
    observations.check_is_present(gene='L', tissue='M', time=2) # active
    observations.check_is_present(gene='R', tissue='N', time=3) # active

    observations.check_is_present(gene='Y', tissue='N', time=4)
    observations.check_is_present(gene='L', tissue='N', time=5)

    observations.check_is_present(gene='R', tissue='M', time=6) # active

    # both on (& staying on)
    observations.check_is_present(gene='Y', tissue='N', time=7)
    observations.check_is_present(gene='Y', tissue='M', time=7)

    # two tissues
    m = Tissue(name='M')
    n = Tissue(name='N')

    # neighbors
    m.add_neighbor(neighbor=n)
    n.add_neighbor(neighbor=m)
    assert n in m.neighbors     # should this be bidirectional? probably. CTB.
    assert m in n.neighbors

    # VFN:
    a = Gene(name='A')          # transient input in M to turn on ligand L
    m.add_gene(gene=a, start=1, duration=6)

    b = Gene(name='B')          # permanent input in M and N to turn on receptor R
    m.add_gene(gene=b, start=1)
    n.add_gene(gene=b, start=1)

    # VFG: 
    ligand = Gene(name='L')          # ligand
    r = Receptor(name='R')      # receptor
    y = Gene(name='Y')          # activated by R

    ligand.activated_or(sources=[a, y])

    r.ligand(activator=b, ligand=ligand) # expression driven by B,
                                         # activated by ligand

    y.activated_by(source=r)    # transcription of Y turned on by activated R.

    # so,
    # pulse of A in M turns on L in M.
    # receptor R is always expressed in N b/c B turns it on.
    # M is neighbor of N, so N sees L.
    # R is activated by L in N.
    # downstream of R, R activates L in N.
    # L then activates R in M.
    # R in M then activates L in M.

    dinkum.run(1, 12)
