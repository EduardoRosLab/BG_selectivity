import time
import nest
import numpy as np

from .basal_ganglia_params import params

nest_seed = 0
np_seed = nest_seed + 1

# Setting up NEST

def restart_kernel(_nest_seed=None, _np_seed=None):
    global nest_seed, np_seed

    nest.ResetKernel()
    nest.SetKernelStatus({'local_num_threads': params.sim.n_threads})
    nest.SetKernelStatus({'resolution': params.sim.step})
    try:
        nest.Install('izhikevich_cond_module')
    except:
        pass

    # Random number generation
    N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
    nest_seed = nest_seed + N_vp + 2 if _nest_seed is None else _nest_seed
    np_seed = nest_seed + N_vp if _np_seed is None else _np_seed

    pyrngs = [np.random.RandomState(s) for s in range(nest_seed, nest_seed + N_vp)]
    nest.SetKernelStatus({'rng_seeds': range(nest_seed, nest_seed + N_vp)})
    nest.SetKernelStatus({'grng_seed': nest_seed + N_vp + 1})
    np.random.seed(np_seed)

restart_kernel()


#TODO: algo más?
