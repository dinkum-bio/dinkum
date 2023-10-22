#
from . import vfg
from . import vfn

import itertools

class Timecourse:
    def iterate(self, *, start=None, stop=None):
        assert start is not None
        assert stop is not None

        tissues = vfn.get_tissues()
        print(f"got {len(tissues)} tissues.")
        for t in tissues:
            print(f"\ttissue {t.name}")

        next_present = set()
        for i in range(start, stop + 1):
            print(f"time {i}")
            for t in tissues:
                this_present = set()
                this_present.update(t.all_present(at=i))
                this_present.update(next_present)

                yield i, t, this_present

                next_present = set()
                for r in vfg.get_rules():
                    for g in r.advance(present=this_present):
                        #print('ZZZ', g.name)
                        next_present.add(g)
