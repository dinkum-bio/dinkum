import pytest
from lmfit import Parameters

import dinkum
from dinkum.exceptions import DinkumInvalidGene
from dinkum.vfg import Gene, GeneStateInfo
from dinkum.vfn import Tissue
from dinkum.vfg2 import (Growth, Decay, LinearCombination, GeneTimecourse,
                         run_lmfit, LogisticRepressor, LogisticActivator,
                         calc_response_1d, calc_response_2d)

from dinkum import observations


def test_basic_gene_timecourse():
    dinkum.reset()

    z = Gene(name='Z')
    m = Tissue(name='M')

    z.custom2(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))

    observations.check_level_is_between(gene='Z',
                                        time=2,
                                        tissue='M',
                                        min_level=33,
                                        max_level=33)
    observations.check_level_is_between(gene='Z',
                                        time=3,
                                        tissue='M',
                                        min_level=66,
                                        max_level=66)
    observations.check_level_is_between(gene='Z',
                                        time=4,
                                        tissue='M',
                                        min_level=100,
                                        max_level=100)

    dinkum.run(1, 5, verbose=True)


def test_basic_decay():
    dinkum.reset()

    x = Gene(name='X')
    m = Tissue(name='M')

    x.custom2(Decay(start_time=1, rate=1.2, initial_level=100, tissue=m))

    observations.check_level_is_between(gene='X',
                                        time=1,
                                        tissue='M',
                                        min_level=100,
                                        max_level=100)
    observations.check_level_is_between(gene='X',
                                        time=2,
                                        tissue='M',
                                        min_level=83,
                                        max_level=84)
    observations.check_level_is_between(gene='X',
                                        time=3,
                                        tissue='M',
                                        min_level=69,
                                        max_level=70)
    observations.check_level_is_between(gene='X',
                                        time=4,
                                        tissue='M',
                                        min_level=57,
                                        max_level=58)
    observations.check_level_is_between(gene='X',
                                        time=5,
                                        tissue='M',
                                        min_level=48,
                                        max_level=49)

    dinkum.run(1, 5, verbose=True)


def test_basic_growth():
    dinkum.reset()

    y = Gene(name='Y')
    m = Tissue(name='M')

    y.custom2(Growth(start_time=1, rate=0.5, initial_level=0, tissue=m))

    observations.check_level_is_between(gene='Y',
                                        time=1,
                                        tissue='M',
                                        min_level=0,
                                        max_level=0)
    observations.check_level_is_between(gene='Y',
                                        time=2,
                                        tissue='M',
                                        min_level=50,
                                        max_level=50)
    observations.check_level_is_between(gene='Y',
                                        time=3,
                                        tissue='M',
                                        min_level=75,
                                        max_level=75)
    observations.check_level_is_between(gene='Y',
                                        time=4,
                                        tissue='M',
                                        min_level=87,
                                        max_level=87)
    observations.check_level_is_between(gene='Y',
                                        time=5,
                                        tissue='M',
                                        min_level=93,
                                        max_level=93)

    dinkum.run(1, 5, verbose=True)


def test_basic_linear_combination():
    dinkum.reset()

    z = Gene(name='Z')
    out = Gene(name='out')
    m = Tissue(name='M')

    z.custom2(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))
    linear_combination = LinearCombination(weights=[1],
                                           gene_names=['Z'])
    out.custom2(linear_combination)

    observations.check_level_is_between(gene='Z',
                                        time=2,
                                        tissue='M',
                                        min_level=33,
                                        max_level=33)
    observations.check_level_is_between(gene='Z',
                                        time=3,
                                        tissue='M',
                                        min_level=66,
                                        max_level=66)
    observations.check_level_is_between(gene='Z',
                                        time=4,
                                        tissue='M',
                                        min_level=100,
                                        max_level=100)

    # should exactly mirror 'Z', just one tick later
    observations.check_level_is_between(gene='out',
                                        time=3,
                                        tissue='M',
                                        min_level=33,
                                        max_level=33)
    observations.check_level_is_between(gene='out',
                                        time=4,
                                        tissue='M',
                                        min_level=66,
                                        max_level=66)
    observations.check_level_is_between(gene='out',
                                        time=5,
                                        tissue='M',
                                        min_level=100,
                                        max_level=100)

    dinkum.run(1, 5, verbose=True)


def test_basic_linear_combination_no_such_gene():
    dinkum.reset()

    z = Gene(name='Z')
    out = Gene(name='out')
    m = Tissue(name='M')

    z.custom2(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))
    linear_combination = LinearCombination(weights=[1, .33, .33],
                                           gene_names=['X', 'Y', 'Z'])
    out.custom2(linear_combination)


    with pytest.raises(DinkumInvalidGene):
        dinkum.run(1, 5, verbose=True)


def test_logistic_activator():
    dinkum.reset()

    z = Gene(name='Z')
    out = Gene(name='out')
    m = Tissue(name='M')

    out.custom2(LogisticActivator(rate=100, midpoint=50, activator_name='Z'))

    xvals, yvals = calc_response_1d(timepoint=2,
                                    target_gene_name='out',
                                    variable_gene_name='Z')
    assert yvals[0] == 0
    assert yvals[47] == 0
    assert yvals[48] == 1
    assert yvals[49] == 9
    assert yvals[50] == 50
    assert yvals[51] == 91
    assert yvals[52] == 99
    assert yvals[53] == 100
    assert yvals[100] == 100


def test_logistic_repressor():
    dinkum.reset()

    x = Gene(name='X')
    z = Gene(name='Z')
    out = Gene(name='out')
    m = Tissue(name='M')

    out.custom2(LogisticRepressor(rate=100, midpoint=50,
                                  activator_name='X',
                                  repressor_name='Z'))

    xvals, yvals = calc_response_1d(timepoint=2,
                                    target_gene_name='out',
                                    variable_gene_name='Z',
                                    fixed_gene_states={ 'X': GeneStateInfo(100, True) })
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

    x = Gene(name='X')
    y = Gene(name='Y')
    z = Gene(name='Z')
    o = Gene(name='out')
    m = Tissue(name='M')

    x.custom2(Decay(start_time=1, rate=1.2, initial_level=100, tissue=m))
    y.custom2(Growth(start_time=1, rate=0.5, initial_level=0, tissue=m))
    z.custom2(GeneTimecourse(start_time=2, tissue=m, values=[100, 200, 300]))

    linear_combination = LinearCombination(weights=[-0.5, .33, .33],
                                           gene_names=['X', 'Y', 'Z'])
    o.custom2(linear_combination)

    fit_values = [0, 100, 83, 69, 57]

    run_lmfit(1, 5, fit_values, [o])

    p = Parameters()
    linear_combination.get_params(p)
    assert round(p['out_wX'].value, 2) == 1.00
    assert round(p['out_wY'].value, 2) == 0.00
    assert round(p['out_wZ'].value, 2) == -0.01
