import dinkum
from dinkum.vfg import Gene, Receptor, Ligand
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import vfg, vfn


def define_model():
    dinkum.reset()

    # gataE, gcm, foxA

    # set it all up!
    pmar1 = Gene(name="pmar1")
    hesC = Gene(name="hesC")
    alx1 = Gene(name="alx1")
    delta = Ligand(name="delta")
    notch = Receptor(name="su(h)", ligand=delta)
    tbr = Gene(name="tbr")
    ets1 = Gene(name="ets1")
    gataE = Gene(name="gataE")
    gcm = Gene(name="gcm")
    foxA = Gene(name="foxA")

    embryo = Tissue(name="rest of embryo")
    micromere = Tissue(name="micromeres")
    embryo.add_neighbor(neighbor=micromere)

    # maternal genes
    early_ubiq = Gene(name="ub1")
    late_ubiq = Gene(name="ub2")

    ## set up maternal gene expression

    # early ubiq
    early_ubiq.is_present(where=micromere, start=1)
    early_ubiq.is_present(where=embryo, start=1)

    # late ubiq
    late_ubiq.is_present(where=micromere, start=4)
    late_ubiq.is_present(where=embryo, start=4)

    # pmar1: maternal in micromeres only
    pmar1.is_present(where=micromere, start=1)

    # hesC: present in both at the beginning
    hesC.is_present(where=micromere, start=1, duration=1)
    hesC.is_present(where=embryo, start=1, duration=1)

    # notch: present everywhere (?)
    notch.is_present(where=micromere, start=1)
    notch.is_present(where=embryo, start=1)

    ## set up all downstream genes

    # hesC: early, if not for pmar1
    hesC.and_not(activator=early_ubiq, repressor=pmar1)

    # alx &c: late, if not for hesC
    alx1.and_not(activator=late_ubiq, repressor=hesC)
    delta.and_not(activator=late_ubiq, repressor=hesC)
    tbr.and_not(activator=late_ubiq, repressor=hesC)
    ets1.and_not(activator=late_ubiq, repressor=hesC)

    # gataE, gcm, foxA: all activated by notch
    gataE.activated_by_and(sources=[late_ubiq, notch])
    gcm.activated_by_and(sources=[late_ubiq, notch])
    foxA.activated_by_and(sources=[late_ubiq, notch])


def test_output_double_neg():
    filename = "/tmp/out.btp.csv"
    define_model()
    print(f"writing to {filename}")
    dinkum.utils.output_biotapestry_csv("/tmp/out.btp.csv")
