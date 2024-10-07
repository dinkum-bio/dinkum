import dinkum
from dinkum.vfg import Gene, Receptor
from dinkum.vfn import Tissue


def test_tissue_cmp():
    dinkum.reset()
    m = Tissue(name='M')
    n = Tissue(name='N')

    assert m < n


def test_gene_cmp():
    dinkum.reset()
    a = Gene(name='M')
    b = Gene(name='N')

    assert a < b

