import pytest
from lmfit import Parameters

import dinkum
from dinkum.exceptions import DinkumInvalidGene
from dinkum.vfg import Gene, GeneStateInfo
from dinkum.vfn import Tissue
from dinkum.vfg_functions import (
    Growth,
    Decay,
    LinearCombination,
    GeneTimecourse,
    run_lmfit,
    LogisticRepressor,
    LogisticActivator,
    calc_response_1d,
    calc_response_2d,
    LogisticMultiRepressor,
    LogisticRepressor2,
)

from dinkum import observations


def test_basic_gene_timecourse():
    dinkum.reset()

    z = Gene(name="Z")
    m = Tissue(name="M")

    z.custom_obj(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))

    observations.check_level_is_between(
        gene="Z", time=2, tissue="M", min_level=33, max_level=33
    )
    observations.check_level_is_between(
        gene="Z", time=3, tissue="M", min_level=66, max_level=66
    )
    observations.check_level_is_between(
        gene="Z", time=4, tissue="M", min_level=100, max_level=100
    )

    dinkum.run(1, 5, verbose=True)


def test_decay_defaults():
    dinkum.reset()

    m = Tissue(name="M")
    fn = Decay(rate=1.2, tissue=m)
    assert fn.initial_level == 100
    assert fn.start_time == 1


def test_basic_decay():
    dinkum.reset()

    x = Gene(name="X")
    m = Tissue(name="M")

    x.custom_obj(Decay(start_time=1, rate=1.2, initial_level=100, tissue=m))

    observations.check_level_is_between(
        gene="X", time=1, tissue="M", min_level=100, max_level=100
    )
    observations.check_level_is_between(
        gene="X", time=2, tissue="M", min_level=83, max_level=84
    )
    observations.check_level_is_between(
        gene="X", time=3, tissue="M", min_level=69, max_level=70
    )
    observations.check_level_is_between(
        gene="X", time=4, tissue="M", min_level=57, max_level=58
    )
    observations.check_level_is_between(
        gene="X", time=5, tissue="M", min_level=48, max_level=49
    )

    dinkum.run(1, 5, verbose=True)


def test_growth_defaults():
    dinkum.reset()

    m = Tissue(name="M")
    fn = Growth(rate=0.25, tissue=m)
    assert fn.initial_level == 0
    assert fn.start_time == 1


def test_basic_growth():
    dinkum.reset()

    y = Gene(name="Y")
    m = Tissue(name="M")

    y.custom_obj(Growth(start_time=1, rate=0.5, initial_level=0, tissue=m))

    observations.check_level_is_between(
        gene="Y", time=1, tissue="M", min_level=0, max_level=0
    )
    observations.check_level_is_between(
        gene="Y", time=2, tissue="M", min_level=50, max_level=50
    )
    observations.check_level_is_between(
        gene="Y", time=3, tissue="M", min_level=75, max_level=75
    )
    observations.check_level_is_between(
        gene="Y", time=4, tissue="M", min_level=87, max_level=87
    )
    observations.check_level_is_between(
        gene="Y", time=5, tissue="M", min_level=93, max_level=93
    )

    dinkum.run(1, 5, verbose=True)


def test_basic_linear_combination():
    dinkum.reset()

    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    z.custom_obj(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))
    linear_combination = LinearCombination(weights=[1], gene_names=["Z"])
    out.custom_obj(linear_combination)

    observations.check_level_is_between(
        gene="Z", time=2, tissue="M", min_level=33, max_level=33
    )
    observations.check_level_is_between(
        gene="Z", time=3, tissue="M", min_level=66, max_level=66
    )
    observations.check_level_is_between(
        gene="Z", time=4, tissue="M", min_level=100, max_level=100
    )

    # should exactly mirror 'Z', just one tick later
    observations.check_level_is_between(
        gene="out", time=3, tissue="M", min_level=33, max_level=33
    )
    observations.check_level_is_between(
        gene="out", time=4, tissue="M", min_level=66, max_level=66
    )
    observations.check_level_is_between(
        gene="out", time=5, tissue="M", min_level=100, max_level=100
    )

    dinkum.run(1, 5, verbose=True)


def test_basic_linear_combination_no_such_gene():
    dinkum.reset()

    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    z.custom_obj(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))
    linear_combination = LinearCombination(
        weights=[1, 0.33, 0.33], gene_names=["X", "Y", "Z"]
    )
    out.custom_obj(linear_combination)

    with pytest.raises(DinkumInvalidGene):
        dinkum.run(1, 5, verbose=True)


def test_logistic_activator_defaults():
    dinkum.reset()

    logit = LogisticActivator(activator_name="Z")
    assert logit.rate == 11
    assert logit.midpoint == 50


def test_logistic_activator():
    dinkum.reset()

    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    out.custom_obj(LogisticActivator(rate=100, midpoint=50, activator_name="Z"))

    xvals, yvals = calc_response_1d(
        timepoint=2, target_gene_name="out", variable_gene_name="Z"
    )
    assert yvals[0] == 0
    assert yvals[47] == 0
    assert yvals[48] == 1
    assert yvals[49] == 9
    assert yvals[50] == 50
    assert yvals[51] == 91
    assert yvals[52] == 99
    assert yvals[53] == 100
    assert yvals[100] == 100


def test_logistic_repressor_defaults():
    dinkum.reset()

    logit = LogisticRepressor(activator_name="Z", repressor_name="X")
    assert logit.rate == 11
    assert logit.midpoint == 50


def test_logistic_repressor():
    dinkum.reset()

    x = Gene(name="X")
    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    out.custom_obj(
        LogisticRepressor(rate=100, midpoint=50, activator_name="X", repressor_name="Z")
    )

    xvals, yvals = calc_response_1d(
        timepoint=2,
        target_gene_name="out",
        variable_gene_name="Z",
        fixed_gene_states={"X": GeneStateInfo(100, True)},
    )
    print(yvals)
    assert yvals[0] == 100
    assert yvals[47] == 100
    assert yvals[48] == 99
    assert yvals[49] == 91
    assert yvals[50] == 50
    assert yvals[51] == 9
    assert yvals[52] == 1
    assert yvals[53] == 0
    assert yvals[100] == 0


def test_logistic_repressor_2d():
    dinkum.reset()

    x = Gene(name="X")
    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    out.custom_obj(
        LogisticRepressor(rate=100, midpoint=50, activator_name="X", repressor_name="Z")
    )

    arr = calc_response_2d(
        timepoint=2, target_gene_name="out", x_gene_name="X", y_gene_name="Z"
    )

    # test values? @CTB


def test_logistic_repressor2_defaults():
    dinkum.reset()

    logit = LogisticRepressor2(activator_name="Z", repressor_name="X")
    assert logit.activator_rate == 11
    assert logit.activator_midpoint == 25
    assert logit.repressor_rate == 11
    assert logit.repressor_midpoint == 75


def test_logistic_repressor2():
    dinkum.reset()

    x = Gene(name="X")
    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    out.custom_obj(
        LogisticRepressor2(
            activator_rate=100,
            activator_name="X",
            repressor_rate=100,
            repressor_name="Z",
        )
    )

    xvals, yvals = calc_response_1d(
        timepoint=2,
        target_gene_name="out",
        variable_gene_name="Z",
        fixed_gene_states={"X": GeneStateInfo(100, True)},
    )
    print(yvals)
    assert yvals[50] == 100
    assert yvals[80] == 0


def test_logistic_repressor2_2d():
    dinkum.reset()

    x = Gene(name="X")
    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    out.custom_obj(
        LogisticRepressor2(
            activator_rate=100,
            activator_name="X",
            repressor_rate=100,
            repressor_name="Z",
        )
    )

    arr = calc_response_2d(
        timepoint=2, target_gene_name="out", x_gene_name="X", y_gene_name="Z"
    )
    assert arr[0, 0] == 0
    assert arr[25, 30] == 100
    assert arr[100, 100] == 0
    assert arr[100, 25] == 0

    # test values? @CTB


def test_logistic_repressor_multi():
    dinkum.reset()

    x = Gene(name="X")
    z = Gene(name="Z")
    out = Gene(name="out")
    m = Tissue(name="M")

    out.custom_obj(
        LogisticMultiRepressor(
            rate=100, midpoint=50, activator_name="X", repressor_names=["Z"]
        )
    )

    xvals, yvals = calc_response_1d(
        timepoint=2,
        target_gene_name="out",
        variable_gene_name="Z",
        fixed_gene_states={"X": GeneStateInfo(100, True)},
    )
    print(yvals)
    assert yvals[0] == 100
    assert yvals[47] == 100
    assert yvals[48] == 99
    assert yvals[49] == 91
    assert yvals[50] == 50
    assert yvals[51] == 9
    assert yvals[52] == 1
    assert yvals[53] == 0
    assert yvals[100] == 0


def test_fit():
    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")
    z = Gene(name="Z")
    o = Gene(name="out")
    m = Tissue(name="M")

    x.custom_obj(Decay(start_time=1, rate=1.2, initial_level=100, tissue=m))
    y.custom_obj(Growth(start_time=1, rate=0.5, initial_level=0, tissue=m))
    z.custom_obj(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))

    linear_combination = LinearCombination(gene_names=["X", "Y", "Z"])
    o.custom_obj(linear_combination)

    fit_values = [0, 100, 83, 69, 57]

    run_lmfit(1, 5, fit_values=fit_values, fit_genes=[o], debug=True)

    p = Parameters()
    linear_combination.get_params(p)
    assert round(p["out_wX"].value, 2) == 1.00
    assert round(p["out_wY"].value, 2) == 0.00
    assert round(p["out_wZ"].value, 2) == -0.01


def test_fit_2_logistic_activator():
    # can we fit parameters to given output for LogisticActivator?

    # first, calculate a response curve
    dinkum.reset()

    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))
    out.custom_obj(LogisticActivator(rate=92, midpoint=58, activator_name="X"))

    tc = dinkum.run(start=1, stop=20, verbose=True)
    states = tc.get_states()
    level_df, _ = states.to_dataframe()

    print(level_df)
    fit_to_vals = list(level_df["out"])
    print(fit_to_vals)

    # ok, now, reset and run fit on growth -
    dinkum.reset()

    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))

    logit = LogisticActivator(activator_name="X")
    out.custom_obj(logit)

    # fit!!
    res = run_lmfit(
        1, 20, fit_values=fit_to_vals, fit_genes=[out], debug=True, method="brute"
    )

    # extract parameters -
    print(res.params)
    logit.set_params(res.params)

    # have some tolerance...
    assert int(logit.rate) in range(80, 95)
    assert int(logit.midpoint) in range(55, 65)


def test_fit_2_growth():
    # can we fit parameters to given output for Growth?

    # first, calculate a response curve
    dinkum.reset()

    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))
    out.custom_obj(LogisticActivator(rate=11, midpoint=58, activator_name="X"))

    tc = dinkum.run(start=1, stop=20, verbose=True)
    states = tc.get_states()
    level_df, _ = states.to_dataframe()

    print(level_df)
    fit_to_vals = list(level_df["out"])
    print(fit_to_vals)

    # ok, now, reset and run fit on growth -
    dinkum.reset()

    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    # use incorrect params for Growth
    growth = Growth(start_time=1, rate=1, initial_level=50, tissue=m)
    x.custom_obj(growth)

    logit = LogisticActivator(rate=11, midpoint=58, activator_name="X")
    out.custom_obj(logit)

    # fit!!
    res = run_lmfit(
        1, 20, fit_values=fit_to_vals, fit_genes=[x], debug=True, method="brute"
    )

    # extract parameters -
    print(res.params)
    growth.set_params(res.params)
    assert growth.rate == 0.1
    assert growth.initial_level == 0


def test_fit_2_decay():
    # can we fit parameters to given output for Decay?

    # first, calculate a response curve
    dinkum.reset()

    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    x.custom_obj(Decay(start_time=1, rate=1.2, initial_level=100, tissue=m))
    out.custom_obj(LogisticActivator(rate=11, midpoint=58, activator_name="X"))

    tc = dinkum.run(start=1, stop=20, verbose=True)
    states = tc.get_states()
    level_df, _ = states.to_dataframe()

    print(level_df)
    fit_to_vals = list(level_df["out"])
    print(fit_to_vals)

    # ok, now, reset and run fit on growth -
    dinkum.reset()

    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    # use incorrect params for Growth
    growth = Decay(start_time=1, rate=1, initial_level=50, tissue=m)
    x.custom_obj(growth)

    logit = LogisticActivator(rate=11, midpoint=58, activator_name="X")
    out.custom_obj(logit)

    # fit!!
    res = run_lmfit(
        1, 20, fit_values=fit_to_vals, fit_genes=[x], debug=True, method="brute"
    )

    # extract parameters -
    print(res.params)
    growth.set_params(res.params)
    assert round(growth.rate, 1) == 1.2
    assert int(growth.initial_level) == 75


def test_fit_2_logistic_repressor():
    # can we fit parameters to given output for LogisticRepressor?

    # first, calculate a response curve
    dinkum.reset()

    ubiq = Gene(name="ubiq")
    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    ubiq.is_present(where=m, start=1)
    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))
    out.custom_obj(
        LogisticRepressor(
            rate=92, midpoint=58, activator_name="ubiq", repressor_name="X"
        )
    )

    tc = dinkum.run(start=1, stop=20, verbose=True)
    states = tc.get_states()
    level_df, _ = states.to_dataframe()

    print(level_df)
    fit_to_vals = list(level_df["out"])
    print(fit_to_vals)

    # ok, now, reset and run fit on growth -
    dinkum.reset()

    ubiq = Gene(name="ubiq")
    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    ubiq.is_present(where=m, start=1)
    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))

    logit = LogisticRepressor(activator_name="ubiq", repressor_name="X")
    out.custom_obj(logit)

    # fit!!
    res = run_lmfit(
        1, 20, fit_values=fit_to_vals, fit_genes=[out], debug=True, method="brute"
    )

    # extract parameters -
    print(res.params)
    logit.set_params(res.params)

    # have some tolerance...
    assert int(logit.rate) in range(80, 95)
    assert int(logit.midpoint) in range(55, 65)


def test_fit_2_logistic_multi_repressor():
    # can we fit parameters to given output for LogisticMultiRepressor?

    # first, calculate a response curve
    dinkum.reset()

    ubiq = Gene(name="ubiq")
    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    ubiq.is_present(where=m, start=1)
    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))
    out.custom_obj(
        LogisticMultiRepressor(
            rate=92, midpoint=58, activator_name="ubiq", repressor_names=["X"]
        )
    )

    tc = dinkum.run(start=1, stop=20, verbose=True)
    states = tc.get_states()
    level_df, _ = states.to_dataframe()

    print(level_df)
    fit_to_vals = list(level_df["out"])
    print(fit_to_vals)

    # ok, now, reset and run fit on growth -
    dinkum.reset()

    ubiq = Gene(name="ubiq")
    x = Gene(name="X")
    out = Gene(name="out")
    m = Tissue(name="M")

    ubiq.is_present(where=m, start=1)
    x.custom_obj(Growth(start_time=1, rate=0.1, initial_level=0, tissue=m))

    logit = LogisticMultiRepressor(activator_name="ubiq", repressor_names=["X"])
    out.custom_obj(logit)

    # fit!!
    res = run_lmfit(
        1, 20, fit_values=fit_to_vals, fit_genes=[out], debug=False, method="brute"
    )

    # extract parameters -
    print(res.params)
    logit.set_params(res.params)

    assert int(logit.rate) == 19  # could change?
    assert round(logit.weights[0], 1) == 0.9
