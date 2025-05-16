import numpy as np
from lmfit import minimize, Parameters
import math

import dinkum
from dinkum import vfg, vfn

# from dinkum import vfg, Timecourse, TissueGeneStates, get_tissue, get_gene
# from dinkum.vfg import GeneStateInfo, check_ligand
from itertools import chain

# @CTB prevent set_gene from being called multiple times


def logistic_output(*, rate, input_level, midpoint):
    "calc logistic function, centered at midpoint, with k = log(rate/10)"
    rate = math.log(rate / 10)
    expon = -rate * (input_level - midpoint)
    expon = min(expon, 50)
    denom = 1 + math.exp(expon)
    output = round(100 / denom)

    return output


class Decay:
    def __init__(self, *, start_time=1, rate, initial_level=100, tissue, delay=1):
        self.start_time = start_time
        self.rate = rate
        self.initial_level = initial_level
        self.tissue = tissue
        self.level = None
        self.delay = delay  # 'delay' is only for checking for ligands

    def get_params(self, params_obj):
        target_name = self.target.name
        decay_name = f"{target_name}_decay"
        initial_name = f"{target_name}_initial"

        params_obj.add(decay_name, value=self.rate, brute_step=0.1, min=0.01)
        params_obj.add(initial_name, value=self.initial_level, brute_step=5)

    def set_params(self, params_obj):
        target_name = self.target.name
        decay_name = f"{target_name}_decay"
        initial_name = f"{target_name}_initial"

        self.rate = params_obj[decay_name].value
        self.initial_level = params_obj[initial_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        # not right tissue, or not started yet? no opinion.
        if tissue != self.tissue or timepoint < self.start_time:
            return None

        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=self.delay,
        )
        if timepoint == self.start_time:
            # start!!
            self.level = self.initial_level
            return self.target, vfg.GeneStateInfo(self.level, active)
        else:
            self.level /= self.rate
            return self.target, vfg.GeneStateInfo(self.level, active)


class Growth:
    def __init__(self, *, start_time=1, rate, initial_level=0, tissue, delay=1):
        self.start_time = start_time
        self.rate = rate
        self.initial_level = initial_level
        self.tissue = tissue
        self.level = None
        self.delay = delay

    def get_params(self, params_obj):
        target_name = self.target.name
        growth_name = f"{target_name}_growth"
        initial_name = f"{target_name}_initial"

        params_obj.add(growth_name, value=self.rate, brute_step=0.1)
        params_obj.add(initial_name, value=self.initial_level, brute_step=5)

    def set_params(self, params_obj):
        target_name = self.target.name
        growth_name = f"{target_name}_growth"
        self.rate = params_obj[growth_name].value

        initial_name = f"{target_name}_initial"
        self.initial_level = params_obj[initial_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        # not right tissue, or not started yet? no opinion.
        if tissue != self.tissue or timepoint < self.start_time:
            return None

        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=self.delay,
        )
        if timepoint == self.start_time:
            # start!!
            self.level = self.initial_level
            return self.target, vfg.GeneStateInfo(self.level, active)
        else:
            self.level += int(100 - self.level) * self.rate
            level = min(self.level, 100.0)
            level = max(level, 0)
            return self.target, vfg.GeneStateInfo(int(level), active)


class GeneTimecourse:
    def __init__(self, *, start_time, tissue, values, normalize=True, delay=1):
        self.start_time = start_time
        self.tissue = tissue
        self.delay = delay  # only for checking for ligands

        values = list(values)
        if normalize:
            max_val = max(values) or 1
            values = [x / max_val * 100 for x in values]

        assert min(values) >= 0, min(values)
        assert max(values) <= 100, max(values)
        self.values = values

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        # not right tissue, or not started yet? no opinion.
        if tissue != self.tissue or timepoint < self.start_time:
            return None

        index = timepoint - self.start_time
        if 0 <= index < len(self.values):
            val = self.values[index]
            active = True
        else:
            val = 0
            active = False

        active = active and vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=self.delay,
        )
        return self.target, vfg.GeneStateInfo(int(val), active)


class LinearCombination:
    def __init__(self, *, weights=None, gene_names, delay=1):
        self.weights = weights
        self.gene_names = list(gene_names)
        self.delay = delay

    def get_params(self, params_obj):
        target_name = self.target.name
        upstream_names = self.gene_names

        weights = self.weights
        if weights is None:
            weights = [1] * len(upstream_names)

        d = {}
        for n, w in zip(upstream_names, weights):
            param_name = f"{target_name}_w{n}"
            params_obj.add(param_name, value=w, min=-20, max=20, brute_step=0.1)

    def set_params(self, params_obj):
        target_name = self.target.name
        upstream_names = self.gene_names

        weights = []
        d = {}
        for n in upstream_names:
            param_name = f"{target_name}_w{n}"
            val = params_obj[param_name].value
            weights.append(val)

        self.weights = weights

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        if not self.weights:
            raise Exception("need weights")
        assert len(self.weights) == len(self.gene_names)

        delay = self.delay
        input_activity = []
        for name in self.gene_names:
            gene = vfg.get_gene(name)
            gsi = states.get_gene_state_info(
                timepoint=timepoint, delay=delay, gene=gene, tissue=tissue
            )
            input_level = 0
            if gsi and gsi.active:
                input_level = gsi.level
            input_activity.append(input_level)

        output = 0
        for weight, activity in zip(self.weights, input_activity):
            output += weight * activity

        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=delay,
        )
        return self.target, vfg.GeneStateInfo(output, active)


class LogisticActivator:
    "Logistic function: switch on above threshold."

    def __init__(self, *, rate=11, midpoint=50, activator_name, delay=1):
        self.rate = rate
        self.midpoint = midpoint
        self.activator_name = activator_name
        self.delay = delay

    def get_params(self, params_obj):
        target_name = self.target.name
        rate = float(self.rate)
        midpoint = float(self.midpoint)

        param_name = f"{target_name}_rate"
        params_obj.add(param_name, value=self.rate, min=11, max=100, brute_step=1)

        param_name = f"{target_name}_midpoint"
        params_obj.add(param_name, value=midpoint, min=0, max=100, brute_step=1)

    def set_params(self, params_obj):
        target_name = self.target.name

        param_name = f"{target_name}_rate"
        self.rate = params_obj[param_name].value

        param_name = f"{target_name}_midpoint"
        self.midpoint = params_obj[param_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        delay = self.delay

        gene = vfg.get_gene(self.activator_name)
        gsi = states.get_gene_state_info(
            timepoint=timepoint, delay=delay, gene=gene, tissue=tissue
        )

        if gsi is not None and gsi.active:
            input_level = gsi.level
        else:
            input_level = 0

        # calc logistic function, centered at midpoint, with k = log(rate/10)
        level = logistic_output(
            rate=self.rate, input_level=input_level, midpoint=self.midpoint
        )

        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=delay,
        )
        return self.target, vfg.GeneStateInfo(level, active)


Activator = LogisticActivator


class LogisticMultiActivator:
    """Logistic function: combine multiple (weighted) inputs, activate
    if logistic function is above threshold.

    Can be used for both AND and OR, with different weights.
    """

    def __init__(
        self,
        *,
        rate=11,
        midpoint=50,
        activator_names,
        weights=None,
        delay=1,
    ):
        self.rate = rate
        self.midpoint = midpoint
        self.delay = delay
        self.activator_names = activator_names

        # default weights to 1 for each input, if not specified
        if weights is None:
            weights = [1] * len(activator_names)

        assert len(weights) == len(activator_names)
        self.weights = weights

    def get_params(self, params_obj):
        target_name = self.target.name
        rate = float(self.rate)

        param_name = f"{target_name}_rate"
        params_obj.add(param_name, value=self.rate, min=11, max=100, brute_step=1)

        # no need to adjust midpoint here; adjusting the weights can suffice.
        upstream_names = self.activator_names

        weights = self.weights
        if weights is None:
            weights = [1] * len(upstream_names)

        d = {}
        for n, w in zip(upstream_names, weights):
            param_name = f"{target_name}_w{n}"
            params_obj.add(param_name, value=w, min=-20, max=20, brute_step=0.1)

    def set_params(self, params_obj):
        target_name = self.target.name

        param_name = f"{target_name}_rate"
        self.rate = params_obj[param_name].value

        upstream_names = self.activator_names

        weights = []
        d = {}
        for n in upstream_names:
            param_name = f"{target_name}_w{n}"
            val = params_obj[param_name].value
            weights.append(val)
        self.weights = weights

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        delay = self.delay

        # calculate weighted sum of upstream
        activator_sum = 0.0
        for activator, weight in zip(self.activator_names, self.weights):
            r = vfg.get_gene(activator)
            activator_state = states.get_gene_state_info(
                timepoint=timepoint, delay=delay, gene=r, tissue=tissue
            )
            # check active here and elsewhere? @CTB YES FIX/TEST ELSEWHERE :)
            if activator_state is not None and activator_state.active:
                activator_sum += weight * activator_state.level

        activator_output = logistic_output(
            rate=self.rate, midpoint=self.midpoint, input_level=activator_sum
        )

        # are we repressed?
        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=delay,
        )
        return self.target, vfg.GeneStateInfo(activator_output, active)


class LogisticRepressor:
    """Logistic function: activate if activator, unless repressor
    above threshold"""

    def __init__(
        self, *, rate=11, midpoint=50, activator_name, repressor_name, delay=1
    ):
        self.rate = rate
        self.midpoint = midpoint
        self.delay = delay
        self.activator = activator_name
        self.repressor = repressor_name

    def get_params(self, params_obj):
        target_name = self.target.name
        rate = float(self.rate)
        midpoint = float(self.midpoint)

        param_name = f"{target_name}_rate"
        params_obj.add(param_name, value=self.rate, min=11, max=100, brute_step=1)

        param_name = f"{target_name}_midpoint"
        params_obj.add(param_name, value=midpoint, min=0, max=100, brute_step=1)

    def set_params(self, params_obj):
        target_name = self.target.name

        param_name = f"{target_name}_rate"
        self.rate = params_obj[param_name].value

        param_name = f"{target_name}_midpoint"
        self.midpoint = params_obj[param_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        delay = self.delay

        activator = vfg.get_gene(self.activator)
        activator_state = states.get_gene_state_info(
            timepoint=timepoint, delay=delay, gene=activator, tissue=tissue
        )

        # are we activated? if not, then bail early.
        if (
            activator_state is None
            or activator_state.level == 0
            or not activator_state.active
        ):
            return self.target, vfg.GeneStateInfo(0, False)

        # ok, activated - record level and now see if we are repressed...
        activator_level = activator_state.level

        repressor = vfg.get_gene(self.repressor)
        repressor_state = states.get_gene_state_info(
            timepoint=timepoint, delay=delay, gene=repressor, tissue=tissue
        )

        if repressor_state is not None:
            repressor_input = repressor_state.level
        else:
            repressor_input = 0

        # calc logistic function, centered at midpoint, with k = log(rate/10)
        repressor_output = logistic_output(
            rate=self.rate, midpoint=self.midpoint, input_level=repressor_input
        )

        # are we repressed?
        level2 = max(activator_level - repressor_output, 0)
        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=delay,
        )
        return self.target, vfg.GeneStateInfo(level2, active)


class LogisticRepressor2:
    """Logistic function: activate if activator, unless repressor
    above threshold. Both activation and repression are switch-like"""

    def __init__(
        self,
        *,
        activator_rate=11,
        repressor_rate=11,
        activator_midpoint=25,
        activator_name,
        repressor_midpoint=75,
        repressor_name,
        delay=1,
    ):
        self.activator_rate = activator_rate
        self.repressor_rate = repressor_rate
        self.activator_midpoint = activator_midpoint
        self.repressor_midpoint = repressor_midpoint
        self.delay = delay
        self.activator = activator_name
        self.repressor = repressor_name

    def get_params(self, params_obj):
        target_name = self.target.name

        param_name = f"{target_name}_activator_rate"
        params_obj.add(
            param_name, value=float(self.activator_rate), min=11, max=100, brute_step=1
        )

        param_name = f"{target_name}_activator_midpoint"
        params_obj.add(
            param_name,
            value=float(self.activator_midpoint),
            min=0,
            max=100,
            brute_step=1,
        )

        param_name = f"{target_name}_repressor_rate"
        params_obj.add(
            param_name, value=float(self.repressor_rate), min=11, max=100, brute_step=1
        )

        param_name = f"{target_name}_repressor_midpoint"
        params_obj.add(
            param_name,
            value=float(self.repressor_midpoint),
            min=0,
            max=100,
            brute_step=1,
        )

    def set_params(self, params_obj):
        target_name = self.target.name

        param_name = f"{target_name}_activator_rate"
        self.activator_rate = params_obj[param_name].value

        param_name = f"{target_name}_activator_midpoint"
        self.activator_midpoint = params_obj[param_name].value

        param_name = f"{target_name}_repressor_rate"
        self.repressor_rate = params_obj[param_name].value

        param_name = f"{target_name}_repressor_midpoint"
        self.repressor_midpoint = params_obj[param_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        delay = self.delay

        activator = vfg.get_gene(self.activator)
        activator_state = states.get_gene_state_info(
            timepoint=timepoint, delay=delay, gene=activator, tissue=tissue
        )

        # are we activated? if not, then bail early.
        activator_level = 0
        if activator_state:
            activator_level = activator_state.level

        activator_level = logistic_output(
            rate=self.activator_rate,
            input_level=activator_level,
            midpoint=self.activator_midpoint,
        )

        # now get repressor level
        repressor = vfg.get_gene(self.repressor)
        repressor_state = states.get_gene_state_info(
            timepoint=timepoint, delay=delay, gene=repressor, tissue=tissue
        )

        repressor_input = 0
        if repressor_state is not None:
            repressor_input = repressor_state.level

        repressor_output = logistic_output(
            rate=self.repressor_rate,
            midpoint=self.repressor_midpoint,
            input_level=repressor_input,
        )

        # are we repressed?
        level2 = max(activator_level - repressor_output, 0)
        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=delay,
        )
        return self.target, vfg.GeneStateInfo(level2, active)


Repressor = LogisticRepressor2


class LogisticMultiRepressor:
    """Logistic function: activate if activator, unless sum of repressors
    above threshold"""

    def __init__(
        self,
        *,
        rate=11,
        midpoint=50,
        activator_name,
        repressor_names,
        weights=None,
        delay=1,
    ):
        self.rate = rate
        self.midpoint = midpoint
        self.delay = delay
        self.activator = activator_name
        assert not isinstance(repressor_names, str)
        self.repressor_names = repressor_names

        # default weights to 1 for each input, if not specified
        if weights is None:
            weights = [1] * len(repressor_names)

        assert len(weights) == len(repressor_names)
        self.weights = weights

    def get_params(self, params_obj):
        target_name = self.target.name
        rate = float(self.rate)

        param_name = f"{target_name}_rate"
        params_obj.add(param_name, value=self.rate, min=11, max=100, brute_step=1)

        # no need to adjust midpoint here; adjusting the weights can suffice.
        upstream_names = self.repressor_names

        weights = self.weights
        if weights is None:
            weights = [1] * len(upstream_names)

        d = {}
        for n, w in zip(upstream_names, weights):
            param_name = f"{target_name}_w{n}"
            params_obj.add(param_name, value=w, min=-20, max=20, brute_step=0.1)

    def set_params(self, params_obj):
        target_name = self.target.name

        param_name = f"{target_name}_rate"
        self.rate = params_obj[param_name].value

        upstream_names = self.repressor_names

        weights = []
        d = {}
        for n in upstream_names:
            param_name = f"{target_name}_w{n}"
            val = params_obj[param_name].value
            weights.append(val)
        self.weights = weights

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        delay = self.delay

        activator = vfg.get_gene(self.activator)
        activator_state = states.get_gene_state_info(
            timepoint=timepoint, delay=delay, gene=activator, tissue=tissue
        )

        # are we activated? if not, then bail early.
        if (
            activator_state is None
            or activator_state.level == 0
            or not activator_state.active
        ):
            return self.target, vfg.GeneStateInfo(0, False)

        # ok, activated - record level and now see if we are repressed...
        activator_level = activator_state.level

        repressor_sum = 0.0
        for repressor, weight in zip(self.repressor_names, self.weights):
            r = vfg.get_gene(repressor)
            repressor_state = states.get_gene_state_info(
                timepoint=timepoint, delay=delay, gene=r, tissue=tissue
            )
            if repressor_state is not None:
                repressor_sum += weight * repressor_state.level

        repressor_output = logistic_output(
            rate=self.rate, midpoint=self.midpoint, input_level=repressor_sum
        )

        # are we repressed?
        level2 = max(activator_level - repressor_output, 0)
        active = vfg.check_ligand(
            dest=self.target,
            timepoint=timepoint,
            states=states,
            tissue=tissue,
            delay=delay,
        )
        return self.target, vfg.GeneStateInfo(level2, True)


def get_ix2_for_gene_name(gene_name):
    for ix in vfg._rules:
        if ix.dest.name == gene_name:
            if not isinstance(ix, vfg.Interaction_CustomObj):
                raise Exception(f"ix {ix} must be a CustomObj ix")
            yield ix


def run_lmfit(start, stop, *, fit_values, fit_genes, debug=False, method="leastsq"):
    for g in fit_genes:
        assert isinstance(g, vfg.Gene)

    def get_fit_params():
        "Extract initial fit parameters from all genes."
        p = Parameters()
        found = set()

        for fit_gene in fit_genes:
            for ix in get_ix2_for_gene_name(fit_gene.name):
                found.add(fit_gene.name)
                obj = ix.obj
                obj.get_params(p)

        # track genes we were supposed to find (but didn't) & complain.
        if len(found) != len(fit_genes):
            missing = set([g.name for g in fit_genes]) - found
            raise Exception(f"missing: {missing}")
        return p

    def set_fit_params(p):
        "Set fit parameters on all the genes."
        for ix in vfg._rules:
            if ix.dest in fit_genes:
                if not isinstance(ix, vfg.Interaction_CustomObj):
                    raise Exception(f"ix {ix} must be a CustomObj ix")
                obj = ix.obj
                obj.set_params(p)

    tc = dinkum.Timecourse(start=start, stop=stop)

    times = list(range(start, stop + 1))

    def residual(params, xvals, data):
        assert xvals == times
        assert data == fit_values

        if debug:
            print("RUNNING RESIDUAL:", xvals)
            params.pretty_print()

        # now, run the time course!
        set_fit_params(params)
        tc.reset()
        tc.run()

        # retrieve values
        vals = []
        for ga in tc:
            for gene in fit_genes:
                gs = ga.get_by_tissue_name("M").get_gene_state(gene.name)
                vals.append(gs.level)

        residuals = np.array(vals) - np.array(data)
        if debug:
            print("results:", residuals)
        return residuals

    params = get_fit_params()
    res = minimize(residual, params, args=(times, fit_values), method=method)
    set_fit_params(res.params)

    if hasattr(res, "message"):
        print("fit message:", res.message)

    print("fit values:")
    for k in res.init_values:
        print(f"\t{k}: fit={res.params[k].value} (was: {res.init_values[k]})")

    return res


def run_lmfit2(start, stop, *, debug=False, method="leastsq", **kwargs):
    fit_genes = {}
    for gene_name in vfg.get_gene_names():
        val = kwargs.get(gene_name)
        if val is not None:
            fit_genes[gene_name] = val
    if not len(fit_genes):
        raise Exception("error! no genes to fit!?")
    if len(kwargs) != len(fit_genes):
        raise Exception("error! incorrect genes names given??")

    print(f"found {len(fit_genes)} genes to fit: {', '.join(fit_genes)}")
    fit_vals = []

    fit_gobj = []
    for k, timecourse in fit_genes.items():
        n_leading = timecourse.start_time - start
        n_trailing = stop - (len(timecourse.values) + n_leading)
        print(n_leading, n_trailing, len(timecourse.values))
        assert n_leading + n_trailing + len(timecourse.values) == stop - start + 1
        fit_vals.append([0] * n_leading + list(timecourse.values) + [0] * n_trailing)
        fit_gobj.append(vfg.get_gene(k))

    # flatten values
    fit_vals = list(chain.from_iterable(zip(*fit_vals)))
    return run_lmfit(
        start, stop, fit_values=fit_vals, fit_genes=fit_gobj, debug=debug, method=method
    )


def calc_response_1d(
    *,
    timepoint=1,
    target_gene_name,
    variable_gene_name,
    fixed_gene_states={},
    delay=1,
    tissue_name="M",
):
    """
    Return two 101-length arrays, xvals and yvals.

    xvals is range(0, 101).
    yvals is the activity of the target_gene_name as the variable_gene_name
    varies from 0 to 100.

    Optionally fix other genes using fixed_gene_states.

    Run at timepoint with given delay.
    """
    tissue = vfn.get_tissue(tissue_name)

    set_tp = timepoint - delay

    states_d = dinkum.TissueGeneStates()
    for gene_name, gsi in fixed_gene_states.items():
        assert isinstance(gsi, vfg.GeneStateInfo)
        states_d.set_gene_state(
            timepoint=set_tp,
            tissue_name=tissue_name,
            gene_name=gene_name,
            state_info=gsi,
        )

    # get relevant interaction for the given target gene
    ixlist = []
    for ix in get_ix2_for_gene_name(target_gene_name):
        ixlist.append(ix)
    assert len(ixlist) == 1
    ix = ixlist[0]

    xvals = []
    yvals = []
    for activity in range(0, 101):
        in_gsi = vfg.GeneStateInfo(activity, True)
        states_d.set_gene_state(
            timepoint=set_tp,
            tissue_name=tissue_name,
            gene_name=variable_gene_name,
            state_info=in_gsi,
        )
        out_gsi = list(ix.advance(timepoint=timepoint, states=states_d, tissue=tissue))
        assert len(out_gsi) == 1
        _, out_gsi = out_gsi[0]

        xvals.append(activity)
        yvals.append(out_gsi.level)

    return np.array(xvals), np.array(yvals)


def calc_response_2d(
    *,
    timepoint=1,
    target_gene_name,
    x_gene_name,
    y_gene_name,
    fixed_gene_states={},
    delay=1,
    tissue_name="M",
):
    """
    Produce a [101, 101]-dimensional array where the Z values are
    the target_gene_name activity as x_gene_name and y_gene_name
    levels vary between 0 and 100.

    Optionally fix gene states for other genes using 'fixed_gene_states'.

    Run the target_gene_name at the given timepoint with the given delay.
    """

    set_tp = timepoint - delay
    tissue = vfn.get_tissue(tissue_name)

    # initialize states dict, optionally with preset fixed_gene states
    states_d = dinkum.TissueGeneStates()
    for gene_name, gsi in fixed_gene_states.items():
        assert isinstance(gsi, vfg.GeneStateInfo)
        states_d.set_gene_state(
            timepoint=set_tp,
            tissue_name=tissue_name,
            gene_name=gene_name,
            state_info=gsi,
        )

    # find the relevant interaction for the target gene; There Should Only Be One
    ixlist = []
    for ix in get_ix2_for_gene_name(target_gene_name):
        ixlist.append(ix)
    assert len(ixlist) == 1
    ix = ixlist[0]

    # iterate across x and y ranges, setting gene activity
    arr = np.zeros((101, 101))
    for x_activity in range(0, 101):
        x_gsi = vfg.GeneStateInfo(x_activity, True)
        states_d.set_gene_state(
            timepoint=set_tp,
            tissue_name=tissue_name,
            gene_name=x_gene_name,
            state_info=x_gsi,
        )
        for y_activity in range(0, 101):
            y_gsi = vfg.GeneStateInfo(y_activity, True)
            states_d.set_gene_state(
                timepoint=set_tp,
                tissue_name=tissue_name,
                gene_name=y_gene_name,
                state_info=y_gsi,
            )
            out_gsi = list(
                ix.advance(timepoint=timepoint, states=states_d, tissue=tissue)
            )
            assert len(out_gsi) == 1
            _, out_gsi = out_gsi[0]

            arr[y_activity, x_activity] = out_gsi.level

    return arr
