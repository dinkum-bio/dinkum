import pytest

import dinkum
from dinkum.exceptions import *
from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations


def test_delete_gene():
    # check that 'reset' clears the gene.
    dinkum.reset()

    x = Gene(name="X")

    dinkum.reset()

    m = Tissue(name="M")
    y = Gene(name="Y")

    # these should both fail
    with pytest.raises(DinkumInvalidGene):
        x.is_present(where=m, start=1)

    with pytest.raises(DinkumInvalidGene):
        y.activated_by(source=x)


def test_delete_tissue():
    # check that 'reset' clears the tissue
    dinkum.reset()

    m = Tissue(name="M")

    dinkum.reset()

    x = Gene(name="X")
    y = Gene(name="Y")

    with pytest.raises(DinkumInvalidTissue):
        x.is_present(where=m, start=1)
