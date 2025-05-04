import dinkum
from dinkum import vfg
from dinkum import vfn
from dinkum.vfg import Gene, Receptor
from dinkum.vfn import Tissue


def test_tissue_cmp():
    dinkum.reset()
    m = Tissue(name="M")
    n = Tissue(name="N")

    assert m < n


def test_gene_cmp():
    dinkum.reset()
    a = Gene(name="M")
    b = Gene(name="N")

    assert a < b


def test_gene_names_with_gene():
    dinkum.reset()
    a = Gene(name="M")
    b = Gene(name="N")
    assert vfg.get_gene_names() == ["M", "N"]


def test_gene_names_with_genes_and_receptor():
    dinkum.reset()
    c = Receptor(name="R")
    a = Gene(name="M")
    b = Gene(name="N")
    assert vfg.get_gene_names() == ["M", "N", "R"]
