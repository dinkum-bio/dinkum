import pytest

import dinkum
from dinkum.exceptions import *
from dinkum.vfg import Gene, Receptor, Ligand
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations


def test_neighbors():
    dinkum.reset()

    m = Tissue(name='M')
    n = Tissue(name='N')

    assert m in m.neighbors     # always self!
    m.add_neighbor(neighbor=n)
    assert n in m.neighbors
    assert m in n.neighbors     # check bidirectional!


def test_neighbors_one_way():
    dinkum.reset()

    m = Tissue(name='M')
    n = Tissue(name='N')

    assert m in m.neighbors     # always self!
    m.add_neighbor(neighbor=n, bidirectional=False)
    assert n in m.neighbors
    assert m not in n.neighbors     # check bidirectional!


def test_signaling_orig_api():
    # CTB deprecate?
    dinkum.reset()

    #observations.check_is_present(gene='R', tissue='M', time=2)
    observations.check_is_present(gene='X', tissue='N', time=2)
    observations.check_is_not_present(gene='Y', tissue='M', time=2)
    observations.check_is_present(gene='Y', tissue='M', time=4)

    m = Tissue(name='M')
    n = Tissue(name='N')

    m.add_neighbor(neighbor=n)

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


def test_signaling_new_api():
    # CTB deprecate?
    # use newer API for Receptor
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
    r = Receptor(name='R', ligand=x)
    y = Gene(name='Y')

    r.ligand(activator=a)

    y.activated_by(source=r)

    # receptor is always present in M
    m.add_gene(gene=a, start=1)

    # x is present in N at time >= 2
    n.add_gene(gene=x, start=2)

    dinkum.run(1, 5)


def test_signaling_new_api_3():
    # use newEST API for Receptor.
    dinkum.reset()

    #observations.check_is_present(gene='R', tissue='M', time=2)
    observations.check_is_present(gene='X', tissue='N', time=2)
    observations.check_is_not_present(gene='Y', tissue='M', time=2)
    observations.check_is_present(gene='Y', tissue='M', time=4)

    m = Tissue(name='M')
    n = Tissue(name='N')

    m.add_neighbor(neighbor=n)
    assert n in m.neighbors     # should this be bidirectional? probably.

    x = Ligand(name='X')
    a = Gene(name='A')
    r = Receptor(name='R', ligand=x)
    y = Gene(name='Y')

    r.ligand(activator=a)

    y.activated_by(source=r)

    # receptor is always present in M
    m.add_gene(gene=a, start=1)

    # x is present in N at time >= 2
    n.add_gene(gene=x, start=2)

    dinkum.run(1, 5)


def test_signaling_ligand_is_not_direct():
    # check that ligands are not allowed to directly regulate
    dinkum.reset()

    x = Ligand(name='X')
    a = Gene(name='A')

    with pytest.raises(DinkumNotATranscriptionFactor):
        a.activated_by(source=x)

    with pytest.raises(DinkumNotATranscriptionFactor):
        a.activated_by_or(sources=[x])

    with pytest.raises(DinkumNotATranscriptionFactor):
        a.activated_by_and(sources=[x])

    with pytest.raises(DinkumNotATranscriptionFactor):
        a.and_not(activator=x, repressor=a)

    with pytest.raises(DinkumNotATranscriptionFactor):
        a.and_not(activator=a, repressor=x)

    with pytest.raises(DinkumNotATranscriptionFactor):
        a.toggle_repressed(tf=x)


def test_signaling_ligand_is_not_direct_custom():
    # check that ligands can't directly activate in custom activation fn
    dinkum.reset()

    m = Tissue(name='M')
    x = Ligand(name='X')
    a = Gene(name='A')

    m.add_gene(gene=x, start=1)

    def activator_fn(*, X):
        return X

    a.custom_activation(state_fn=activator_fn, delay=1)
    with pytest.raises(DinkumNotATranscriptionFactor):
        dinkum.run(1, 12)


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

    ligand.activated_by_or(sources=[a, y])

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
