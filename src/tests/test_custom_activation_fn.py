import pytest

import dinkum
from dinkum.vfg import Gene, CustomActivation
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations

from dinkum.exceptions import *
from dinkum.log import *


def test_custom_fn_1():
    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    def activator_fn(*, X):
        if X.level == 0 or not X.active:
            assert bool(X) == False
        return X

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=activator_fn, delay=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")

    observations.check_is_not_present(gene="X", time=3, tissue="M")
    observations.check_is_not_present(gene="Y", time=3, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_fn_1_delay_2():
    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    def activator_fn(*, X):
        if X.level == 0 or not X.active:
            assert bool(X) == False
        return X

    x.is_present(where=m, start=1, duration=2)
    y.custom_fn(state_fn=activator_fn, delay=2)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="X", time=3, tissue="M")

    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_not_present(gene="Y", time=2, tissue="M")
    observations.check_is_present(gene="Y", time=3, tissue="M")
    observations.check_is_present(gene="Y", time=4, tissue="M")
    observations.check_is_not_present(gene="Y", time=5, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_fn_2():
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")
    m = Tissue(name="M")

    def activator_fn(*, Z, X):  # allow order independence
        return X

    x.is_present(where=m, start=1, duration=1)
    z.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=activator_fn, delay=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")

    observations.check_is_not_present(gene="X", time=3, tissue="M")
    observations.check_is_not_present(gene="Y", time=3, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_1_fail():
    # no gene Z
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    def activator_fn(*, Z):
        return X

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=activator_fn, delay=1)

    # run time course
    with pytest.raises(DinkumInvalidGene):
        tc = dinkum.run(1, 5)


def test_custom_bad_defn():
    # invalid activation function defns - must be kwargs
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    # must be explicitly named, not positional
    def activator_fn(Z):
        return X

    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_fn(state_fn=activator_fn, delay=1)

    # no defaults
    def activator_fn(Z=None):
        return X

    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_fn(state_fn=activator_fn, delay=1)

    # no defaults 2
    def activator_fn(*, Z=None):
        return X

    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_fn(state_fn=activator_fn, delay=1)

    # no general kwargs
    def activator_fn(**kwargs):
        return X

    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_fn(state_fn=activator_fn, delay=1)

    # no general args
    def activator_fn(*args):
        return X

    with pytest.raises(DinkumInvalidActivationFunction):
        y.custom_fn(state_fn=activator_fn, delay=1)


def test_custom_fail_bad_return():
    # basic does-it-work - test bad return value
    set_debug(True)

    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    def activator_fn(*, X):  # allow order independence
        return 100, False, "something else"

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=activator_fn, delay=1)

    # run time course; expect error
    with pytest.raises(DinkumInvalidActivationResult):
        dinkum.run(1, 5)


def test_custom_class_1():
    # does it work with a custom class? give a list of gene names
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    class ActivateMe(CustomActivation):
        def __call__(self, *, X):
            return X

    state_fn = ActivateMe(input_genes=["X"])

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=state_fn, delay=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")

    observations.check_is_not_present(gene="X", time=3, tissue="M")
    observations.check_is_not_present(gene="Y", time=3, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_class_1_bad_args():
    # test that custom classes fail with bad activator fn __call__
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    class ActivateMe(CustomActivation):
        def __call__(self, foo, *, X):
            return X

    with pytest.raises(DinkumInvalidActivationFunction):
        state_fn = ActivateMe(input_genes=["X"])


def test_custom_class_2():
    # does it work with a custom class? give a list of genes.
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    class ActivateMe(CustomActivation):
        def __call__(self, *, X):
            return X

    state_fn = ActivateMe(input_genes=[x])

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=state_fn, delay=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")

    observations.check_is_not_present(gene="X", time=3, tissue="M")
    observations.check_is_not_present(gene="Y", time=3, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_class_3():
    # does it work with a custom class? retrieve gene names from __call__
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    class ActivateMe(CustomActivation):
        def __call__(self, *, X):
            return X

    state_fn = ActivateMe()

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=state_fn, delay=1)

    # set observations
    observations.check_is_present(gene="X", time=1, tissue="M")
    observations.check_is_not_present(gene="X", time=2, tissue="M")
    observations.check_is_not_present(gene="Y", time=1, tissue="M")
    observations.check_is_present(gene="Y", time=2, tissue="M")

    observations.check_is_not_present(gene="X", time=3, tissue="M")
    observations.check_is_not_present(gene="Y", time=3, tissue="M")

    # run time course
    tc = dinkum.run(1, 5)
    assert len(tc) == 5


def test_custom_class_1_fail():
    # various failure modes
    set_debug(True)

    # basic does-it-work
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    m = Tissue(name="M")

    class ActivateMe(CustomActivation):
        def __call__(self, *, Z):
            return X

    state_fn = ActivateMe()

    x.is_present(where=m, start=1, duration=1)
    y.custom_fn(state_fn=state_fn, delay=1)

    # run time course
    with pytest.raises(DinkumInvalidGene):
        tc = dinkum.run(1, 5)
