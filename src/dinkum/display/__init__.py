"Generate basic activity recording."

from .widgets import *

from .. import Timecourse

def tc_record_activity(*, start=1, stop=10, gene_names=None, verbose=False,
                       trace_fn=None):
    tc = Timecourse(start=start, stop=stop, trace_fn=trace_fn)

    state_record = []     # (tp_name, state)

    time_points = {}      # time_point_name => index
    all_tissues = set()   # all tissues across all time points

    tc.run()

    # iterate over timecourses, pulling out state information.
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

    # build a function (closure-ish) that returns the gene state.
    def get_gene_state(tissue_name, time_point, gene):
        time_idx = time_points[time_point]
        ga = state_record[time_idx]
        gs = ga[1].get_by_tissue_name(tissue_name).get_gene_state(gene)
        return gs

    return state_record, list(all_tissues), get_gene_state
