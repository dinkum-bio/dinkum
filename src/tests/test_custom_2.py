import dinkum
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum.vfg2 import Growth, Decay, LinearCombination, GeneTimecourse, \
    run_lmfit

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

    d = dict(linear_combination.get_params())
    assert round(d['out_wX'], 2) == 1.00
    assert round(d['out_wY'], 2) == 0.00
    assert round(d['out_wZ'], 2) == -0.01

