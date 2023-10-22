#! /usr/bin/env python
import sys
import dinkum
import argparse

from dinkum.vfg import Gene
from dinkum.vfn import Tissue
from dinkum import Timecourse
from dinkum import observations

def main():
    # establish preconditions
    observations.check_is_present(gene='X', time=1, tissue='M')
    observations.check_is_not_present(gene='Y', time=1, tissue='M')

    observations.check_is_present(gene='X', time=2, tissue='M')
    observations.check_is_present(gene='Y', time=2, tissue='M')
    observations.check_is_present(gene='Z', time=2, tissue='M')

    observations.check_is_not_present(gene='Z', time=3, tissue='M')

    # set it all up!
    x = Gene(name='X')
    y = Gene(name='Y')
    z = Gene(name='Z')

    y.activated_by(source=x)
    x.activated_by(source=y)

    z.and_not(activator=x, repressor=y)

    m = Tissue(name='M')
    x.is_present(where=m, start=1, duration=1)

    tc = Timecourse()
    succ = True
    for state in tc.iterate(start=1, stop=5):
        print(f"time={state.time}")
        for ti in state.tissues:
            present = state[ti]
            print(f"\ttissue={ti.name}, {[ g.name for g in present ]}")
        observations.test_observations(state)
        succ = succ and observations.test_observations(state)

    if not succ:
        print("** FAILED **")

    return 0 if succ else 1


if __name__ == '__main__':
    sys.exit(main())
