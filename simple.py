import sys
import dinkum
import argparse

from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse

def main():
    x = Gene(name='X')
    y = Gene(name='Y')

    y.activated_by(source=x)
    x.activated_by(source=y)

    m = Tissue(name='M')
    x.is_present(where=m, start=1)

    tc = Timecourse()
    for state in tc.iterate(start=1, stop=5):
        print(f"time={state.time}")
        for ti in state.tissues:
            present = state[ti]
            print(f"\ttissue={ti.name}, {[ g.name for g in present ]}")


if __name__ == '__main__':
    sys.exit(main())
