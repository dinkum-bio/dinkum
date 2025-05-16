import pytest

import dinkum
from dinkum.exceptions import *
from dinkum.vfg import Gene, Receptor, Ligand
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations


def test_neighbors():
    dinkum.reset()

    m = Tissue(name="M")
    n = Tissue(name="N")

    assert m in m.neighbors  # always self!
    m.add_neighbor(neighbor=n)
    assert n in m.neighbors
    assert m in n.neighbors  # check bidirectional!


def test_neighbors_one_way():
    dinkum.reset()

    m = Tissue(name="M")
    n = Tissue(name="N")

    assert m in m.neighbors  # always self!
    m.add_neighbor(neighbor=n, bidirectional=False)
    assert n in m.neighbors
    assert m not in n.neighbors  # check bidirectional!


def test_ligand_receptor():
    dinkum.reset()

    x = Gene(name="X")
    with pytest.raises(DinkumNotALigand):
        r = Receptor(name="R", ligand=x)


def test_signaling():
    dinkum.reset()

    # receptor present, but not active at time 2 in M:
    observations.check_is_present(gene="R", tissue="M", time=2)
    observations.check_is_not_active(gene="R", tissue="M", time=2)

    # ligand present at time 2 in N
    observations.check_is_present(gene="X", tissue="N", time=2)

    # receptor present _and_ active at time 3 in M:
    observations.check_is_present(gene="R", tissue="M", time=3)
    observations.check_is_active(gene="R", tissue="M", time=3)

    # Y is off at t=2, on at t=4
    observations.check_is_not_present(gene="Y", tissue="M", time=2)
    observations.check_is_present(gene="Y", tissue="M", time=4)

    m = Tissue(name="M")
    n = Tissue(name="N")

    m.add_neighbor(neighbor=n)

    x = Ligand(name="X")
    a = Gene(name="A")
    r = Receptor(name="R", ligand=x)
    y = Gene(name="Y")

    r.activated_by(source=a)
    y.activated_by(source=r)

    # receptor is always present in M
    m.add_gene(gene=a, start=1)

    # x is present in N at time >= 2
    n.add_gene(gene=x, start=2)

    def trace_me(**kw):
        print(kw)

    dinkum.run(1, 5, trace_fn=trace_me)


def test_signaling_2():
    dinkum.reset()

    observations.check_is_present(gene="R", tissue="M", time=2)
    observations.check_is_present(gene="X", tissue="N", time=2)
    observations.check_is_not_present(gene="Y", tissue="M", time=2)
    observations.check_is_present(gene="Y", tissue="M", time=4)

    m = Tissue(name="M")
    n = Tissue(name="N")

    m.add_neighbor(neighbor=n)
    assert n in m.neighbors  # should this be bidirectional? probably.

    x = Ligand(name="X")
    a = Gene(name="A")
    r = Receptor(name="R", ligand=x)
    y = Gene(name="Y")

    r.activated_by(source=a)
    y.activated_by(source=r)

    # receptor is always present in M
    m.add_gene(gene=a, start=1)

    # x is present in N at time >= 2
    n.add_gene(gene=x, start=2)

    def trace_me(**kw):
        print(kw)

    dinkum.run(1, 5, trace_fn=trace_me)


def test_signaling_3():
    # use newest API for Receptor.
    dinkum.reset()

    # observations.check_is_present(gene='R', tissue='M', time=2)
    observations.check_is_present(gene="X", tissue="N", time=2)
    observations.check_is_not_present(gene="Y", tissue="M", time=2)
    observations.check_is_present(gene="Y", tissue="M", time=4)

    m = Tissue(name="M")
    n = Tissue(name="N")

    m.add_neighbor(neighbor=n)
    assert n in m.neighbors  # should this be bidirectional? probably.

    x = Ligand(name="X")
    a = Gene(name="A")
    r = Receptor(name="R", ligand=x)
    y = Gene(name="Y")

    x.activated_by(source=a)
    y.activated_by(source=r)

    m.add_gene(gene=a, start=1)

    # receptor is always present in M
    m.add_gene(gene=r, start=1)

    # x is present in N at time >= 2
    n.add_gene(gene=x, start=2)

    def trace_me(**kw):
        print(kw)

    dinkum.run(1, 5, trace_fn=trace_me)


def test_signaling_new_api_3_default_no_juxtacrine():
    # use newest API for Receptor; check that by default signals
    # signal back to receptors in same cell.
    dinkum.reset()

    observations.check_is_not_present(gene="Y", tissue="N", time=1)
    observations.check_is_not_present(gene="Y", tissue="M", time=1)
    observations.check_is_present(gene="Y", tissue="N", time=3)
    observations.check_is_present(gene="Y", tissue="M", time=3)

    m = Tissue(name="M")
    n = Tissue(name="N")

    m.add_neighbor(neighbor=n)
    assert n in m.neighbors

    x = Ligand(name="X")
    r = Receptor(name="R", ligand=x)
    y = Gene(name="Y")

    y.activated_by(source=r)

    # receptor is always present in M and N
    m.add_gene(gene=r, start=1)
    n.add_gene(gene=r, start=1)

    # add ligand at t=1
    m.add_gene(gene=x, start=1)

    # now, y should be turned on by time 3 by activated receptor

    dinkum.run(1, 5)


def test_signaling_new_api_3_juxtacrine():
    # use newest API for Receptor; check that by default signals
    # signal back to receptors in same cell.
    dinkum.reset()

    observations.check_is_not_present(gene="Y", tissue="N", time=1)
    observations.check_is_not_present(gene="Y", tissue="M", time=1)
    observations.check_is_present(gene="Y", tissue="N", time=3)
    observations.check_is_not_present(gene="Y", tissue="M", time=3)

    m = Tissue(name="M")
    n = Tissue(name="N")

    m.add_neighbor(neighbor=n)
    assert n in m.neighbors

    x = Ligand(name="X", is_juxtacrine=True)
    r = Receptor(name="R", ligand=x)
    y = Gene(name="Y")

    y.activated_by(source=r)

    # receptor is always present in M and N
    m.add_gene(gene=r, start=1)
    n.add_gene(gene=r, start=1)

    # add ligand at t=1
    m.add_gene(gene=x, start=1)

    # now, y should be turned on by time 3 by activated receptor

    dinkum.run(1, 5)


def test_signaling_ligand_is_not_direct():
    # check that ligands are not allowed to directly regulate
    dinkum.reset()

    x = Ligand(name="X")
    a = Gene(name="A")

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


def test_signaling_ligand_is_not_direct_custom():
    # check that ligands can't directly activate in custom activation fn
    dinkum.reset()

    m = Tissue(name="M")
    x = Ligand(name="X")
    a = Gene(name="A")

    m.add_gene(gene=x, start=1)

    def activator_fn(*, X):
        return X

    a.custom_fn(state_fn=activator_fn, delay=1)
    with pytest.raises(DinkumNotATranscriptionFactor):
        dinkum.run(1, 12)


def test_community_effect():
    # transient input in one cell => mutual lock on via positive feedback/
    # signalling
    # good question for students: how long does pulse need to be, and why??
    dinkum.reset()

    observations.check_is_present(gene="A", tissue="M", time=6)
    observations.check_is_not_present(gene="A", tissue="M", time=7)

    observations.check_is_not_present(gene="L", tissue="M", time=1)  # not act
    observations.check_is_present(gene="L", tissue="M", time=2)  # active
    observations.check_is_active(gene="R", tissue="N", time=3)  # active

    observations.check_is_present(gene="Y", tissue="N", time=4)
    observations.check_is_present(gene="L", tissue="N", time=5)

    observations.check_is_active(gene="R", tissue="M", time=6)  # active

    # both on (& staying on)
    observations.check_is_present(gene="Y", tissue="N", time=7)
    observations.check_is_present(gene="Y", tissue="M", time=7)

    # two tissues
    m = Tissue(name="M")
    n = Tissue(name="N")

    # neighbors
    m.add_neighbor(neighbor=n)
    n.add_neighbor(neighbor=m)
    assert n in m.neighbors
    assert m in n.neighbors

    # VFN:
    a = Gene(name="A")  # transient input in M to turn on ligand L
    m.add_gene(gene=a, start=1, duration=6)

    b = Gene(name="B")  # permanent input in M and N to turn on receptor R
    m.add_gene(gene=b, start=1)
    n.add_gene(gene=b, start=1)

    # VFG:
    ligand = Ligand(name="L")  # ligand
    r = Receptor(name="R", ligand=ligand)  # receptor
    y = Gene(name="Y")  # activated by R

    ligand.activated_by_or(sources=[a, y])

    y.activated_by(source=r)  # transcription of Y turned on by activated R.

    r.activated_by(source=b)

    # so,
    # pulse of A in M turns on L in M.
    # receptor R is always expressed in N b/c B turns it on.
    # M is neighbor of N, so N sees L.
    # R is activated by L in N.
    # downstream of R, R activates L in N.
    # L then activates R in M.
    # R in M then activates L in M.

    def trace_me(**kw):
        print(kw)

    dinkum.run(1, 12, trace_fn=trace_me)
