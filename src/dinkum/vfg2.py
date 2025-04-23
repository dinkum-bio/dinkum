import numpy as np
from lmfit import minimize, Parameters
import math

from dinkum import vfg, Timecourse
from dinkum.vfg import GeneStateInfo

# @CTB prevent set_gene from being called multiple times


class Decay:
    def __init__(self, *, start_time, rate, initial_level, tissue):
        self.start_time = start_time
        self.rate = rate
        self.initial_level = initial_level
        self.tissue = tissue
        self.level = None

    def get_params(self, params_obj):
        target_name = self.target.name
        decay_name = f'{target_name}_decay'
        initial_name = f'{target_name}_initial'

        params_obj.add(decay_name, value=self.rate)
        params_obj.add(initial_name, value=self.initial_level)

    def set_params(self, params_obj):
        target_name = self.target.name
        decay_name = f'{target_name}_decay'
        initial_name = f'{target_name}_initial'

        self.rate = params_obj[decay_name].value
        self.initial_level = params_obj[initial_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        # not right tissue, or not started yet? no opinion.
        if tissue != self.tissue or timepoint < self.start_time:
            return None            

        if timepoint == self.start_time:
            # start!!
            self.level = self.initial_level
            return self.target, GeneStateInfo(self.level, True)
        else: 
            self.level /= self.rate
            return self.target, GeneStateInfo(self.level, True)

class Growth:
    def __init__(self, *, start_time, rate, initial_level, tissue):
        self.start_time = start_time
        self.rate = rate
        self.initial_level = initial_level
        self.tissue = tissue
        self.level = None

    def get_params(self, params_obj):
        target_name = self.target.name
        growth_name = f'{target_name}_growth'
        initial_name = f'{target_name}_initial'

        params_obj.add(growth_name, value=self.rate)
        params_obj.add(initial_name, value=self.initial_level)

    def set_params(self, params_obj):
        target_name = self.target.name
        growth_name = f'{target_name}_growth'
        self.rate = params_obj[growth_name].value

        initial_name = f'{target_name}_initial'
        self.initial_level = params_obj[initial_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        # not right tissue, or not started yet? no opinion.
        if tissue != self.tissue or timepoint < self.start_time:
            return None            
        
        if timepoint == self.start_time:
            # start!!
            self.level = self.initial_level
            return self.target, GeneStateInfo(self.level, True)
        else: 
            self.level += int(100 - self.level) * self.rate
            level = min(self.level, 100.0)
            level = max(level, 0)
            return self.target, GeneStateInfo(int(level), True)

class GeneTimecourse:
    def __init__(self, *, start_time, tissue, values, scale=True):
        self.start_time = start_time
        self.tissue = tissue

        values = list(values)
        if scale:
            max_val = max(values) or 1
            values = [ x / max_val * 100 for x in values ]

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

        return self.target, GeneStateInfo(int(val), active)


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
            param_name = f'{target_name}_w{n}'
            params_obj.add(param_name, value=w)

    def set_params(self, params_obj):
        target_name = self.target.name
        upstream_names = self.gene_names

        weights = []
        d = {}
        for n in upstream_names:
            param_name = f'{target_name}_w{n}'
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
            gsi = states.get_gene_state_info(timepoint, delay, gene, tissue)
            input_level = 0
            if gsi and gsi.active:
                input_level = gsi.level
            input_activity.append(input_level)


        output = 0
        for weight, activity in zip(self.weights, input_activity):
            output += weight*activity

        #output = max(output, 0)
        #output = min(output, 100)
        return self.target, GeneStateInfo(output, True)


class LogisticFunction:
    "Logistic function: switch on above threshold."
    def __init__(self, *, rate=11, midpoint=50, upstream_name=None, delay=1):
        assert upstream_name is not None
        self.rate = rate
        self.midpoint = midpoint
        self.upstream_name = upstream_name
        self.delay = delay

    def get_params(self, params_obj):
        target_name = self.target.name
        rate = float(self.rate)
        midpoint = float(self.midpoint)

        param_name = f'{target_name}_rate'
        params_obj.add(param_name, value=self.rate, min=11, max=100,
                       brute_step=1)

        param_name = f'{target_name}_midpoint'
        params_obj.add(param_name, value=midpoint, min=0, max=100,
                       brute_step=1)

    def set_params(self, params_obj):
        target_name = self.target.name

        param_name = f'{target_name}_rate'
        self.rate = params_obj[param_name].value

        param_name = f'{target_name}_midpoint'
        self.midpoint = params_obj[param_name].value

    def set_gene(self, gene):
        self.target = gene

    def advance(self, timepoint, states, tissue):
        delay = self.delay

        gene = vfg.get_gene(self.upstream_name)
        gsi = states.get_gene_state_info(timepoint, delay, gene, tissue)

        if gsi is not None:
            input_level = gsi.level
        else:
            input_level = 0

        # calc logistic function, centered at midpoint, with k = log(rate/10)
        rate = math.log(self.rate / 10)
        expon = -rate * (input_level - self.midpoint)
        denom = 1 + math.exp(expon)
        level = round(100 / denom)

        return self.target, GeneStateInfo(level, True)


def run_lmfit(start, stop, fit_values, fit_genes,
              *, debug=False, method='leastsq'):
    for g in fit_genes:
        assert isinstance(g, vfg.Gene)

    def get_fit_params():
        "Extract initial fit parameters from all genes."
        p = Parameters()
        found = set()
        for ix in vfg._rules:
            if ix.dest in fit_genes:
                if not isinstance(ix, vfg.Interaction_Custom2):
                    raise Exception(f"ix {ix} must be a Custom2 ix")
                found.add(ix.dest.name)
                obj = ix.obj
                obj.get_params(p)

        # track genes we were supposed to find (but didn't) & complain.
        if len(found) != len(fit_genes):
            missing = set([ g.name for g in fit_genes ]) - found
            raise Exception(f'missing: {missing}')
        return p

    def set_fit_params(p):
        "Set fit parameters on all the genes."
        for ix in vfg._rules:
            if ix.dest in fit_genes:
                if not isinstance(ix, vfg.Interaction_Custom2):
                    raise Exception(f"ix {ix} must be a Custom2 ix")
                obj = ix.obj
                obj.set_params(p)

    tc = Timecourse(start=start, stop=stop)

    times = list(range(start, stop + 1))
    def residual(params, xvals, data):
        assert xvals == times
        assert data == fit_values

        if debug:
            print('RUNNING RESIDUAL:', xvals)
            params.pretty_print()

        # now, run the time course!
        set_fit_params(params)
        tc.reset()
        tc.run()

        # retrieve values
        vals = []
        for ga in tc:
            for gene in fit_genes:
                gs = ga.get_by_tissue_name('M').get_gene_state(gene.name)
                vals.append(gs.level)

        residuals = np.array(vals) - np.array(data)
        if debug:
            print('results:', residuals)
        return residuals

    params = get_fit_params()
    res = minimize(residual, params, args=(times, fit_values),
                   method=method)
    set_fit_params(res.params)

    if hasattr(res, 'message'):
        print('fit message:', res.message)

    print('fit values:')
    for k in res.init_values:
        print(f'\t{k}: fit={res.params[k].value} (was: {res.init_values[k]})')

    return res
