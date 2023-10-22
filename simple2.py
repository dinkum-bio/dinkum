import sys
import dinkum
import argparse

from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse

def main():
    x = Gene(name='X')
    y = Gene(name='Y')
    z = Gene(name='Z')

    y.activated_by(source=x)
    x.activated_by(source=y)

    z.and_not(activator=x, repressor=y)

    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    t = Timecourse()
    t.iterate(start=1, stop=5)


if __name__ == '__main__':
    sys.exit(main())
