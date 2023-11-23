from .widgets import *

from .. import Timecourse

def tc_record_activity(*, start=1, stop=10, gene_names=None, verbose=False):
    tc = Timecourse(start=start, stop=stop)

    state_record = []     # (tp_name, state)

    time_points = {}      # time_point_name => index
    all_tissues = set()   # all tissues across all time points

    tc.run()

    for n, state in enumerate(iter(tc)):
        tp = f"t={state.time}"
        if verbose:
            print(tp)
        time_points[tp] = n

        for ti in state.tissues:
            all_tissues.add(ti.name)
            present = state[ti]
            if verbose:
                print(f"\ttissue={ti.name}, {present.report_activity()}")

        state_record.append((tp, state))

    def is_active(tissue_name, time_point, gene):
        time_idx = time_points[time_point]
        ga = state_record[time_idx]
        active = ga[1].get_by_tissue_name(tissue_name).is_active(gene)
        return bool(active)
    
    return state_record, list(all_tissues), is_active
