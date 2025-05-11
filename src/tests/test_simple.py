import pytest

import dinkum
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations
from dinkum.exceptions import *


def test_maternal():
    dinkum.reset()

    x = Gene(name="X")
    m = Tissue(name="M")
    x.is_present(where=m, start=1, duration=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="X", time=2, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_maternal_fail():
    dinkum.reset()

    x = Gene(name="X")
    m = Tissue(name="M")
    x.is_present(where=m, start=1, duration=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_present(gene="X", time=2, tissue="M")

    with pytest.raises(dinkum.DinkumObservationFailed):
        dinkum.run(1, 5, verbose=True)


def test_activation():
    dinkum.reset()

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")

    y.activated_by(source=x)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")

    # run!
    dinkum.run(start=1, stop=5)


def test_activation_fail():
    dinkum.reset()

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")

    y.activated_by(source=x)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # run!
    tc = dinkum.run(start=1, stop=5)

    # set observations
    observations.check_is_not_present(gene="X", time=1, tissue="M")
    with pytest.raises(dinkum.DinkumObservationFailed):
        tc.check()
    observations.reset()

    observations.check_is_present(gene="Y", time=1, tissue="M")
    with pytest.raises(dinkum.DinkumObservationFailed):
        tc.check()
    observations.reset()

    observations.check_is_not_present(gene="X", time=2, tissue="M")
    with pytest.raises(dinkum.DinkumObservationFailed):
        tc.check()
    observations.reset()

    observations.check_is_not_present(gene="Y", time=2, tissue="M")
    with pytest.raises(dinkum.DinkumObservationFailed):
        tc.check()
    observations.reset()


def test_simple_repression():
    dinkum.reset()

    # establish preconditions
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")

    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")
    observations.check_is_present(gene="Z", time=2, tissue="M")

    observations.check_is_not_present(gene="Z", time=3, tissue="M")

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")

    y.activated_by(source=x)

    z.and_not(activator=x, repressor=y)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # run!
    dinkum.run(1, 5, verbose=True)


def test_simple_multiple_tissues():
    dinkum.reset()

    ## tissue M
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_never_present(gene="Y", tissue="M")

    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Z", time=3, tissue="M")

    ## tissue N
    observations.check_is_present(gene="Y", time=2, tissue="N")
    observations.check_is_present(gene="X", time=2, tissue="N")
    observations.check_is_present(gene="X", time=3, tissue="N")
    observations.check_is_present(gene="Y", time=3, tissue="N")
    observations.check_is_never_present(gene="Z", tissue="N")

    # set it all up!

    ## VFG
    a = Gene(name="A")
    b = Gene(name="B")

    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")

    x.activated_by(source=a)
    y.activated_by(source=b)
    z.and_not(activator=x, repressor=y)

    ## VFN
    m = Tissue(name="M")
    a.is_present(where=m, start=1)

    n = Tissue(name="N")
    a.is_present(where=n, start=1)
    b.is_present(where=n, start=1)

    # run!
    dinkum.run(1, 5)


def test_simple_positive_feedback():
    dinkum.reset()

    # establish preconditions
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=2, tissue="M")

    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")

    observations.check_is_present(gene="X", time=4, tissue="M")
    observations.check_is_present(gene="Y", time=4, tissue="M")

    # set it all up!
    a = Gene(name="A")
    x = Gene(name="X")
    y = Gene(name="Y")

    y.activated_by(source=x)
    x.activated_by_or(sources=[a, y])

    m = Tissue(name="M")
    a.is_present(where=m, start=1, duration=2)

    # run!
    dinkum.run(1, 5, verbose=True)


def test_simple_positive_feedback_old_name_activated_or():
    dinkum.reset()

    # establish preconditions
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=2, tissue="M")

    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")

    observations.check_is_present(gene="X", time=4, tissue="M")
    observations.check_is_present(gene="Y", time=4, tissue="M")

    # set it all up!
    a = Gene(name="A")
    x = Gene(name="X")
    y = Gene(name="Y")

    y.activated_by(source=x)
    # this is an old name for 'activated_by_or'. Check that it works.
    x.activated_or(sources=[a, y])

    m = Tissue(name="M")
    a.is_present(where=m, start=1, duration=2)

    # run!
    dinkum.run(1, 5)


def test_simple_coherent_feed_forward():
    dinkum.reset()

    # establish preconditions
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_not_present(gene="Z", time=1, tissue="M")

    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")
    observations.check_is_not_present(gene="Z", time=2, tissue="M")

    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")
    observations.check_is_present(gene="Z", time=3, tissue="M")

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")

    y.activated_by(source=x)
    z.activated_by_and(sources=[x, y])

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # run!
    dinkum.run(1, 5, verbose=True)


def test_simple_incoherent_feed_forward():
    dinkum.reset()

    # establish preconditions
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_not_present(gene="Z", time=1, tissue="M")

    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")
    observations.check_is_present(gene="Z", time=2, tissue="M")

    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")
    observations.check_is_not_present(gene="Z", time=3, tissue="M")

    observations.check_is_present(gene="X", time=4, tissue="M")
    observations.check_is_present(gene="Y", time=4, tissue="M")
    observations.check_is_not_present(gene="Z", time=4, tissue="M")

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")

    y.activated_by(source=x)
    z.and_not(activator=x, repressor=y)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # run!
    dinkum.run(1, 5)


def test_simple_incoherent_feed_forward_2_tissues():
    dinkum.reset()

    # establish preconditions
    observations.check_is_not_present(gene="Z", time=1, tissue="M")
    observations.check_is_present(gene="Z", time=2, tissue="M")
    observations.check_is_present(gene="Z", time=3, tissue="M")
    observations.check_is_present(gene="Z", time=4, tissue="M")

    observations.check_is_never_present(gene="Y", tissue="M")

    observations.check_is_not_present(gene="Z", time=1, tissue="N")
    observations.check_is_present(gene="Z", time=2, tissue="N")
    observations.check_is_not_present(gene="Z", time=3, tissue="N")
    observations.check_is_not_present(gene="Z", time=4, tissue="N")

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")
    s = Gene(name="S")  # switches spec states b/t tissues M and N

    y.activated_by_and(sources=[x, s])
    z.and_not(activator=x, repressor=y)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    m = Tissue(name="N")
    x.is_present(where=m, start=1)
    s.is_present(where=m, start=1)

    # run!
    dinkum.run(1, 5, verbose=True)


def test_simple_mutual_repression():
    dinkum.reset()

    # establish preconditions
    observations.check_is_never_present(gene="X", tissue="N")
    observations.check_is_never_present(gene="Y", tissue="M")

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")
    a = Gene(name="A")
    b = Gene(name="B")

    x.and_not(activator=a, repressor=y)
    y.and_not(activator=b, repressor=x)

    m = Tissue(name="M")
    a.is_present(where=m, start=1)

    n = Tissue(name="N")
    b.is_present(where=n, start=1)

    # run!
    dinkum.run(1, 5)


def test_delayed_activation():
    dinkum.reset()

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")

    y.activated_by(source=x, delay=2)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=2, tissue="M")
    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")

    # run!
    dinkum.run(start=1, stop=5)


def test_delayed_activation_trace():
    # test trace function
    dinkum.reset()

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")

    y.activated_by(source=x, delay=2)

    m = Tissue(name="M")
    x.is_present(where=m, start=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=2, tissue="M")
    observations.check_is_present(gene="X", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")

    x = []

    def trace_me(*, gene=None, state_info=None, tp=None, tissue=None):
        assert gene is not None
        assert state_info is not None
        assert tp is not None
        assert tissue is not None

        x.append((tp, gene, state_info, tissue))

    # run!
    dinkum.run(start=1, stop=5, trace_fn=trace_me)

    x.sort()
    print(x)
    assert len(x) == 10, len(x)

    tp, gene, state_info, tissue = x[0]
    assert gene.name == "X"
    assert str(state_info) == "<level=100,active=True>"
    assert tissue.name == "M"
    assert tp == 1


def test_multiple_rules_error():
    dinkum.reset()

    # set it all up!
    x = Gene(name="X")
    y = Gene(name="Y")

    x.activated_by(source=x, delay=2)
    with pytest.raises(DinkumMultipleRules):
        x.activated_by(source=y, delay=2)
