"Generate basic activity recording."

from .widgets import *

from .. import Timecourse, TissueGeneStates, _run


def tc_record_activity(
    *, start=1, stop=10, gene_names=None, verbose=False, trace_fn=None
):
    """Execute time course and return states dictionary.

    Legacy notebook function.
    """
    tc = _run(start=start, stop=stop, verbose=verbose, trace_fn=trace_fn)

    states = TissueGeneStates()

    tc.run()
    tc.check()

    # iterate over timecourses, pulling out state information.
    for n, state in enumerate(iter(tc)):
        tp = f"t={state.time}"
        if verbose:
            print(tp)

        for ti in state.tissues:
            present = state[ti]
            if verbose:
                print(f"\ttissue={ti.name}, {present.report_activity()}")

        states[state.time] = state

    return states
