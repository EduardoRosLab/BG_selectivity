import nest
from matplotlib import pyplot as plt
import numpy as np
import itertools

from .utils import raster_plot, verbose, transient_selectivity, steadystate_selectivity#, log_progress
from tqdm import tqdm
from .basal_ganglia import BasalGanglia
from .basal_ganglia_params import params as _params
from . import restart_kernel


class Experiment:

    time_bin = 20.0
    interpolation_kind = ['linear', 'cubic'][1]

    def __init__(self):
        raise NotImplementedError("'Experiment' is a virtual class, you need to inherit it.")

    def run(self):
        raise NotImplementedError("'Experiment' is a virtual class, you need to inherit it.")

    def plot(self):
        raise NotImplementedError("'Experiment' is a virtual class, you need to inherit it.")


    def transient_selectivity(s, nucleus='msn'):
        return transient_selectivity(s.bg.signal(nucleus,0), s.bg.signal(nucleus,1), t_stable=(2000, 2500), t_transit=(1500, 1700))


    def steadystate_selectivity(s, nucleus='msn'):
        return steadystate_selectivity(s.bg.signal(nucleus,1), t_pre=(1250, 1500), t_stable=(2000, 2500))


class ExperimentDAHD(Experiment):

    # Frecuencias de disparo de la neurona cortical
    cortex_low_freq  = 2.2
    cortex_high_freq = 3.6



    def __init__(s, dopamine_level=0.3, huntington_level=0.0, basal_ganglia=None):

        # Initial values
        restart_kernel()
        s.bg = BasalGanglia() if basal_ganglia is None else basal_ganglia

        # Setting dopamine
        s.dopamine_level = dopamine_level
        s.bg.params.neuron.msn_d1.y1 = dopamine_level
        s.bg.params.neuron.msn_d1.y2 = dopamine_level
        s.bg.params.neuron.msn_d2.y1 = dopamine_level
        s.bg.params.neuron.msn_d2.y2 = dopamine_level

        # Setting Huntington
        s.huntington_level = huntington_level
        s.bg.params.syn.ctx_msn.syn_spec_nmda.weight = (3.05 + 3.05*huntington_level) / s.bg.params.neuron.msn_d1.tau_syn[1]
        s.bg.params.net.msn.d2_n = int(s.bg.params.net.msn.n * 0.5 * (1.0 - huntington_level * 0.8))

        # Configurando el modelo
        s.bg.configure()


    @verbose(_params.sim.verbose)
    def run(s, warming_up_time=1500):
        m = s.cortex_low_freq * s.bg.params.net.ctx_afferents_per_neuron
        M = s.cortex_high_freq * s.bg.params.net.ctx_afferents_per_neuron

        # Estado inicial, se mantiene por dos segundos
        for channel in s.bg.poisson:
            nest.SetStatus(channel, {'rate': m})
        nest.SetStatus(s.bg.poisson[0], {'rate': m})
        nest.SetStatus(s.bg.poisson[1], {'rate': m})
        nest.Simulate(warming_up_time)

        # Primera mitad de la rampa, durante 25 ms
        for i in range(25):
            fi = i/50.0
            nest.SetStatus(s.bg.poisson[0], {'rate': m + (M-m)*fi})
            nest.SetStatus(s.bg.poisson[1], {'rate': m + (M-m)*fi})
            nest.Simulate(1.0)

        # Segunda mitad de la rampa, durante 25 ms
        for i in range(25):
            fi = i/50.0
            nest.SetStatus(s.bg.poisson[0], {'rate': (m+M)*0.5 + (M-m)*fi})
            nest.SetStatus(s.bg.poisson[1], {'rate': (m+M)*0.5 - (M-m)*fi})
            nest.Simulate(1.0)

        # Finalización, 1000 ms
        nest.Simulate(1000.0)


    @verbose(_params.sim.verbose)
    def plot(s):
        s._plot_network()
        plt.show()

    @verbose(_params.sim.verbose)
    def _plot_network(s):
        plt.rcParams['figure.figsize'] = [3.0, 8.0]
        for i in range(len(_params.net.p_channels)):
            raster_plot.from_device(s.bg.spikedetector_msn_channels[i], hist_binwidth=s.time_bin, title='Canal {}'.format(i + 1), ylim=(0, 40))



class ExperimentSSTS(Experiment):

    # Frecuencias de disparo de la neurona cortical
    cortex_low_freq  = 2.2
    cortex_high_freq = 3.6
    time_range = (0, 3550)


    def __init__(s):
        resolution=8

        s.da_values = np.linspace(start=0.0, stop=1.0, num=resolution)
        s.hd_values = np.linspace(start=0.0, stop=1.0, num=resolution)
        s.pd_values = np.linspace(start=0.0, stop=1.0, num=resolution)

        s.signal_hd_msn = np.zeros((s.hd_values.size, s.da_values.size, len(_params.net.p_channels), s.time_range[1]-s.time_range[0]))
        s.signal_hd_snr = np.zeros((s.hd_values.size, s.da_values.size, len(_params.net.p_channels), s.time_range[1]-s.time_range[0]))
        s.signal_pd_msn = np.zeros((s.hd_values.size, s.da_values.size, len(_params.net.p_channels), s.time_range[1]-s.time_range[0]))
        s.signal_pd_snr = np.zeros((s.hd_values.size, s.da_values.size, len(_params.net.p_channels), s.time_range[1]-s.time_range[0]))

        s.results_ts_hd_msn = np.zeros((s.hd_values.size, s.da_values.size))
        s.results_ss_hd_msn = np.zeros((s.hd_values.size, s.da_values.size))
        s.results_ts_pd_msn = np.zeros((s.pd_values.size, s.da_values.size))
        s.results_ss_pd_msn = np.zeros((s.pd_values.size, s.da_values.size))
        s.results_ts_hd_snr = np.zeros((s.hd_values.size, s.da_values.size))
        s.results_ss_hd_snr = np.zeros((s.hd_values.size, s.da_values.size))
        s.results_ts_pd_snr = np.zeros((s.pd_values.size, s.da_values.size))
        s.results_ss_pd_snr = np.zeros((s.pd_values.size, s.da_values.size))

        s.bg = None
        s.iterations = None


    def _configure_new_run(s, dopamine_level=0.3, huntington_level=0.0):
        restart_kernel()

        # Initial values
        s.bg = BasalGanglia()

        # Setting the model
        s.bg.params.net.snr.n = 3000 #More neurons in SNr to have more stable measurements

        # Setting dopamine
        s.bg.params.neuron.msn_d1.y1 = dopamine_level
        s.bg.params.neuron.msn_d1.y2 = dopamine_level
        s.bg.params.neuron.msn_d2.y1 = dopamine_level
        s.bg.params.neuron.msn_d2.y2 = dopamine_level

        # Setting Huntington
        s.bg.params.syn.ctx_msn.syn_spec_nmda.weight = (3.05 + 3.05*huntington_level) / s.bg.params.neuron.msn_d1.tau_syn[1]
        s.bg.params.net.msn.d2_n = int(s.bg.params.net.msn.n * 0.5 * (1.0 - huntington_level*0.8))

        # Configurando el modelo
        s.bg.configure()


    def _run_one(s):
        m = s.cortex_low_freq * s.bg.params.net.ctx_afferents_per_neuron
        M = s.cortex_high_freq * s.bg.params.net.ctx_afferents_per_neuron

        # Estado inicial, se mantiene por dos segundos
        for channel in s.bg.poisson:
            nest.SetStatus(channel, {'rate': m})
        nest.SetStatus(s.bg.poisson[0], {'rate': m})
        nest.SetStatus(s.bg.poisson[1], {'rate': m})
        nest.Simulate(2500.0)

        # Primera mitad de la rampa, durante 25 ms
        for i in range(25):
            fi = i/50.0
            nest.SetStatus(s.bg.poisson[0], {'rate': m + (M-m)*fi})
            nest.SetStatus(s.bg.poisson[1], {'rate': m + (M-m)*fi})
            nest.Simulate(1.0)

        # Segunda mitad de la rampa, durante 25 ms
        for i in range(25):
            fi = i/50.0
            nest.SetStatus(s.bg.poisson[0], {'rate': (m+M)*0.5 + (M-m)*fi})
            nest.SetStatus(s.bg.poisson[1], {'rate': (m+M)*0.5 - (M-m)*fi})
            nest.Simulate(1.0)

        # Finalización, 1000 ms
        nest.Simulate(1000.0)


    def transient_selectivity(s, nucleus='msn'):
        return transient_selectivity(s.bg.signal(nucleus,0), s.bg.signal(nucleus,1), t_stable=(3000, 3500), t_transit=(2500, 2700))


    def steadystate_selectivity(s, nucleus='msn'):
        return steadystate_selectivity(s.bg.signal(nucleus,1), t_pre=(2000, 2400), t_stable=(3000, 3500))


    def run(s, times=1, plot=False):
        i_times = 1.0/times
        s.iterations = times if s.iterations==None else s.iterations+times
        conditions_hd = itertools.product(enumerate(s.hd_values), enumerate(s.da_values), range(times))
        conditions_pd = itertools.product(enumerate(s.hd_values), enumerate(s.da_values), range(times))

        for (i, hd), (j, da), _ in tqdm(conditions_hd, total=s.hd_values.size*s.da_values.size*times):
            s._configure_new_run(dopamine_level=da, huntington_level=hd)
            s._run_one()

            for idx_ch, _ in enumerate(s.bg.params.net.p_channels):
                s.signal_hd_msn[j,i,idx_ch,:] += s.bg.signal('msn',idx_ch,time_range=s.time_range)[:] * i_times
                s.signal_hd_snr[j,i,idx_ch,:] += s.bg.signal('snr',idx_ch,time_range=s.time_range)[:] * i_times

            s.results_ts_hd_msn[j,i] += s.transient_selectivity(nucleus='msn') * i_times
            s.results_ts_hd_snr[j,i] += s.transient_selectivity(nucleus='snr') * i_times
            s.results_ss_hd_msn[j,i] += s.steadystate_selectivity(nucleus='msn') * i_times
            s.results_ss_hd_snr[j,i] += s.steadystate_selectivity(nucleus='snr') * i_times

            if plot:
                plt.rcParams['figure.figsize'] = [8.0, 3.0]
                plt.figure()
                s.bg.plot_activity('msn', channels=True, time_bin=20, frequency=True)
                _ = plt.title('DA: {}\nHD: {}'.format(da, hd))

    def plot(s):
        pass
