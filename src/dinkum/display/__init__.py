"Generate basic activity recording."

from .widgets import *

from .. import Timecourse, TissueGeneStates

def tc_record_activity(*, start=1, stop=10, gene_names=None, verbose=False,
                       trace_fn=None):
    # @CTB deprecate this approach!
    tc = Timecourse(start=start, stop=stop, trace_fn=trace_fn)

    states = TissueGeneStates()
    state_record = []     # (tp_name, state)

    time_points = {}      # time_point_name => index
    all_tissues = set()   # all tissues across all time points

    tc.run()
    tc.check()

    # iterate over timecourses, pulling out state information.
    for n, state in enumerate(iter(tc)):
        tp = f"t={state.time}"
        if verbose:
            print(tp)

        for ti in state.tissues:
            all_tissues.add(ti.name)
            present = state[ti]
            if verbose:
                print(f"\ttissue={ti.name}, {present.report_activity()}")

        states[state.time] = state

    return states
