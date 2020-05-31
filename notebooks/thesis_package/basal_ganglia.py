import nest
import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt
from scipy import signal
from operator import itemgetter

from .basal_ganglia_params import params as default_params
from .utils import verbose, split_in_groups, raster_plot
from .utils.custom_exceptions import ParamError


class BasalGanglia:

    def __init__(s, _params=None):
        s.params = _params if _params is not None else default_params.copy()


    def configure(s):
        s._create_nodes()
        s._create_connections()
        s._create_detectors()


    @verbose(default_params.sim.verbose)
    def _create_nodes(s):
        # Creamos los generadores, las neuronas corticales (parrots) y las neuronas del estriado
        s.poisson = [nest.Create('poisson_generator', 1, {'rate': 0.0}) for _ in range(len(s.params.net.p_channels))]

        # Creamos el estriado
        s.str_channels = [[] for _ in s.params.net.p_channels]
        if s.params.net.msn.is_active:
            s._create_msn()

        if s.params.net.stn.is_active:
            s._create_stn()
            nest.Create('parrot_neuron', 10)
        if s.params.net.gpe.is_active:
            s._create_gpe()
            nest.Create('parrot_neuron', 10)
        if s.params.net.snr.is_active:
            s._create_snr()

        # Humphries 2010, conexiones uno-a-uno, fig 9 pag 12. De la entrada cortical sólo importa el número de
        # spikes por segundo que se recibe.
        if s.params.net.msn.is_active:
            s.parrots_str = [nest.Create('parrot_neuron', len(channel)) for channel in s.str_channels]

        if s.params.net.stn.is_active:
            s.parrots_stn = [nest.Create('parrot_neuron', len(channel)) for channel in s.stn_channels]



    def _create_stn(s):
        # Queremos crear un sustituto de este núcleo que dispare a una frecuencia determinada
        if type(s.params.net.stn.is_active) is not bool:
            s.stn_poisson = nest.Create('poisson_generator', 1, {'rate': s.params.net.stn.is_active})
            s.stn = nest.Create('parrot_neuron', s.params.net.stn_llrs.n + s.params.net.stn_nr.n + s.params.net.stn_rb.n)
            nest.Connect(s.stn_poisson, s.stn)
            s.stn_channels = split_in_groups(s.stn, s.params.net.p_channels)

        # Queremos crear el núcleo normal
        else:
            s.stn = []
            s.stn_channels = [[] for _ in s.params.net.p_channels]

            if s.params.net.stn_rb.is_active:
                s.stn_rb = nest.Create("stn", s.params.net.stn_rb.n)
                nest.SetStatus(s.stn_rb, s.params.neuron.stn_rb)
                s.stn_rb_channels = split_in_groups(s.stn_rb, s.params.net.p_channels)

                s.stn += s.stn_rb
                s.stn_channels = [base + extra for base,extra in zip(s.stn_channels, s.stn_rb_channels)]

            if s.params.net.stn_llrs.is_active:
                s.stn_llrs = nest.Create("stn", s.params.net.stn_llrs.n)
                nest.SetStatus(s.stn_llrs, s.params.neuron.stn_llrs)
                s.stn_llrs_channels = split_in_groups(s.stn_llrs, s.params.net.p_channels)

                s.stn += s.stn_llrs
                s.stn_channels = [base + extra for base, extra in zip(s.stn_channels, s.stn_llrs_channels)]

            if s.params.net.stn_nr.is_active:
                s.stn_nr = nest.Create("stn", s.params.net.stn_nr.n)
                nest.SetStatus(s.stn_nr, s.params.neuron.stn_nr)
                s.stn_nr_channels = split_in_groups(s.stn_nr, s.params.net.p_channels)

                s.stn += s.stn_nr
                s.stn_channels = [base + extra for base, extra in zip(s.stn_channels, s.stn_nr_channels)]


    def _create_gpe(s):

        # Queremos crear un sustituto de este núcleo que dispare a una frecuencia determinada
        if type(s.params.net.gpe.is_active) is not bool:
            s.gpe_poisson = nest.Create('poisson_generator', 1, {'rate': s.params.net.gpe.is_active})
            s.gpe = nest.Create('parrot_neuron', s.params.net.gpe_a.n + s.params.net.gpe_b.n + s.params.net.gpe_c.n)
            nest.Connect(s.gpe_poisson, s.gpe)
            s.gpe_channels = split_in_groups(s.gpe, s.params.net.p_channels)

        # Queremos crear el núcleo normal
        else:
            s.gpe = []
            s.gpe_channels = [[] for _ in s.params.net.p_channels]

            if s.params.net.gpe_a.is_active:
                s.gpe_a = nest.Create("gpe", s.params.net.gpe_a.n)
                nest.SetStatus(s.gpe_a, s.params.neuron.gpe_a)
                s.gpe_a_channels = split_in_groups(s.gpe_a, s.params.net.p_channels)

                s.gpe += s.gpe_a
                for to_channel, from_channel in zip(s.gpe_channels, s.gpe_a_channels):
                    to_channel += from_channel

            if s.params.net.gpe_b.is_active:
                s.gpe_b = nest.Create("gpe", s.params.net.gpe_b.n)
                nest.SetStatus(s.gpe_b, s.params.neuron.gpe_b)
                s.gpe_b_channels = split_in_groups(s.gpe_b, s.params.net.p_channels)

                s.gpe += s.gpe_b
                for to_channel, from_channel in zip(s.gpe_channels, s.gpe_b_channels):
                    to_channel += from_channel

            if s.params.net.gpe_c.is_active:
                s.gpe_c = nest.Create("gpe", s.params.net.gpe_c.n)
                nest.SetStatus(s.gpe_c, s.params.neuron.gpe_c)
                s.gpe_c_channels = split_in_groups(s.gpe_c, s.params.net.p_channels)

                s.gpe += s.gpe_c
                for to_channel, from_channel in zip(s.gpe_channels, s.gpe_c_channels):
                    to_channel += from_channel

    def _create_msn(s):
        # Queremos crear un sustituto de este núcleo que dispare a una frecuencia determinada
        if type(s.params.net.msn.is_active) is not bool:
            s.msn_d1 = nest.Create('parrot_neuron', int(s.params.net.msn.d1_n * s.params.net.msn.reduction_factor))
            s.msn_d2 = nest.Create('parrot_neuron', int(s.params.net.msn.d2_n * s.params.net.msn.reduction_factor))
            s.msn = s.msn_d1 + s.msn_d2
            s.msn_poisson = nest.Create('poisson_generator', 1, {'rate': s.params.net.msn.is_active})
            nest.Connect(s.msn_poisson, s.msn)

            s.msn_d1_channels = split_in_groups(s.msn_d1, s.params.net.p_channels)
            s.msn_d2_channels = split_in_groups(s.msn_d2, s.params.net.p_channels)
            s.msn_channels = [d1+d2 for (d1,d2) in zip(s.msn_d1_channels, s.msn_d2_channels)]
            s.str_channels = deepcopy(s.msn_channels)

        # Queremos crear el núcleo normal
        else:
            s.msn_d1 = nest.Create('msn', s.params.net.msn.d1_n)
            s.msn_d2 = nest.Create('msn', s.params.net.msn.d2_n)
            nest.SetStatus(s.msn_d1, s.params.neuron.msn_d1)
            nest.SetStatus(s.msn_d2, s.params.neuron.msn_d2)
            s.msn = s.msn_d1 + s.msn_d2

            s.msn_d1_channels = split_in_groups(s.msn_d1, s.params.net.p_channels)
            s.msn_d2_channels = split_in_groups(s.msn_d2, s.params.net.p_channels)
            s.msn_channels = [d1+d2 for (d1,d2) in zip(s.msn_d1_channels, s.msn_d2_channels)]

            for i in range(len(s.str_channels)):
                s.str_channels[i] += s.msn_channels[i]

            # Establezco valores iniciales de las neuronas. Pongo el voltaje al mínimo
            dUVms = [{"V_m": s.params.neuron.msn_d1.vp} for _ in s.msn]
            nest.SetStatus(s.msn, dUVms)


    def _create_snr(s):
        s.snr = nest.Create("izhikevich_cond", s.params.net.snr.n)
        nest.SetStatus(s.snr, s.params.neuron.snr)
        s.snr_channels = split_in_groups(s.snr, s.params.net.p_channels)


    @verbose(default_params.sim.verbose)
    def _create_connections(s):
        channels = list(range(len(s.params.net.p_channels)))

        # Cortex "afferents" (poisson)
        for i in channels:
            if s.params.net.msn.is_active: nest.Connect(s.poisson[i], s.parrots_str[i])
            if s.params.net.stn.is_active: nest.Connect(s.poisson[i], s.parrots_stn[i])

        # MSN efferents
        if s.params.net.msn.is_active:

            if type(s.params.net.msn.is_active) is bool:
                # Connect parrots (ctx) to msn
                for i in channels:
                    nest.Connect(s.parrots_str[i][:len(s.msn_channels[i])], s.msn_channels[i],
                                 conn_spec=s.params.syn.ctx_msn.conn_spec, syn_spec=s.params.syn.ctx_msn.syn_spec_ampa)
                    nest.Connect(s.parrots_str[i][:len(s.msn_channels[i])], s.msn_channels[i],
                                 conn_spec=s.params.syn.ctx_msn.conn_spec, syn_spec=s.params.syn.ctx_msn.syn_spec_nmda)


                # Connect msn to itself, depending on the type of connectivity defined (random or physical)
                s._connect_msn_to_msn()

                # Connect MSN to SNr, depending on the type of connectivity (plastic or static)
                if s.params.net.snr.is_active:
                    syn_spec = s.params.syn.msn_snr.syn_spec_static
                    for i in channels:
                        nest.Connect(s.msn_d1_channels[i], s.snr_channels[i], conn_spec=s.params.syn.msn_snr.conn_spec, syn_spec=syn_spec)

                # Connect MSN to GPe, depending on the type of connectivity (plastic or static)
                if type(s.params.net.gpe.is_active) is bool and s.params.net.gpe.is_active == True:
                    syn_spec = s.params.syn.msn_gpe.syn_spec_static

                    for i in channels:
                        nest.Connect(s.msn_d2_channels[i], s.gpe_channels[i], conn_spec=s.params.syn.msn_gpe.conn_spec, syn_spec=syn_spec)

            else:
                # Connect MSN to SNr, depending on the type of connectivity (plastic or static)
                if s.params.net.snr.is_active:
                    syn_spec = s.params.syn.msn_snr.syn_spec_static

                    conn_spec = s.params.syn.msn_snr.conn_spec.copy()
                    conn_spec.p = conn_spec.p / s.params.net.msn.reduction_factor

                    for i in channels:
                        nest.Connect(s.msn_d1_channels[i], s.snr_channels[i], conn_spec=conn_spec, syn_spec=syn_spec)

                # Connect MSN to GPe, depending on the type of connectivity (plastic or static)
                if type(s.params.net.gpe.is_active) is bool and s.params.net.gpe.is_active == True:
                    syn_spec = s.params.syn.msn_gpe.syn_spec_static

                    conn_spec = s.params.syn.msn_gpe.conn_spec.copy()
                    conn_spec.p = conn_spec.p / s.params.net.msn.reduction_factor

                    for i in channels:
                        nest.Connect(s.msn_d2_channels[i], s.gpe_channels[i], conn_spec=conn_spec, syn_spec=syn_spec)

        # STN efferents
        if s.params.net.stn.is_active:
            # Connect parrots (ctx) to stn
            if type(s.params.net.stn.is_active) is bool:
                for i in channels:
                    nest.Connect(s.parrots_stn[i], s.stn_channels[i],
                                 conn_spec=s.params.syn.ctx_stn.conn_spec, syn_spec=s.params.syn.ctx_stn.syn_spec_ampa)
                    nest.Connect(s.parrots_stn[i], s.stn_channels[i],
                                 conn_spec=s.params.syn.ctx_stn.conn_spec, syn_spec=s.params.syn.ctx_stn.syn_spec_nmda)

            # Conexion difusa STN -< GPe
            if type(s.params.net.gpe.is_active) is bool and s.params.net.gpe.is_active:
                nest.Connect(s.stn, s.gpe, conn_spec=s.params.syn.stn_gpe.conn_spec, syn_spec=s.params.syn.stn_gpe.syn_spec_ampa)
                nest.Connect(s.stn, s.gpe, conn_spec=s.params.syn.stn_gpe.conn_spec, syn_spec=s.params.syn.stn_gpe.syn_spec_nmda)

            # Conexión difusa STN -> SNr
            if s.params.net.snr.is_active:
                nest.Connect(s.stn, s.snr, conn_spec=s.params.syn.stn_snr.conn_spec, syn_spec=s.params.syn.stn_snr.syn_spec_static_ampa)
                nest.Connect(s.stn, s.snr, conn_spec=s.params.syn.stn_snr.conn_spec, syn_spec=s.params.syn.stn_snr.syn_spec_static_nmda)

        # GPe efferents
        if s.params.net.gpe.is_active:

            # Solo añadimos la conexión lateral inhibitoria si este núcleo no es un sustituto, sino un nucleo real
            if type(s.params.net.gpe.is_active) is bool and s.params.net.gpe.has_local_collaterals:
                for i in channels:
                    nest.Connect(s.gpe_channels[i], s.gpe_channels[i], conn_spec=s.params.syn.gpe_gpe.conn_spec, syn_spec=s.params.syn.gpe_gpe.syn_spec)

            if type(s.params.net.stn.is_active) is bool and s.params.net.stn.is_active:
                for i in channels:
                    nest.Connect(s.gpe_channels[i], s.stn_channels[i], conn_spec=s.params.syn.gpe_stn.conn_spec, syn_spec=s.params.syn.gpe_stn.syn_spec)

            if s.params.net.snr.is_active:
                for i in channels:
                    nest.Connect(s.gpe_channels[i], s.snr_channels[i], conn_spec=s.params.syn.gpe_snr.conn_spec, syn_spec=s.params.syn.gpe_snr.syn_spec_static)

        # SNr efferents
        if s.params.net.snr.is_active:
            nest.Connect(s.snr, s.snr, conn_spec=s.params.syn.snr_snr.conn_spec, syn_spec=s.params.syn.snr_snr.syn_spec)


    def _connect_msn_to_msn(s):
        nest.Connect(s.msn, s.msn, conn_spec=s.params.syn.msn_msn.conn_spec_random, syn_spec=s.params.syn.msn_msn.syn_spec_random)


    @verbose(default_params.sim.verbose)
    def _create_detectors(s):
        s.spikedetector_global = nest.Create('spike_detector', params={"withgid": True, "withtime": True})

        if s.params.net.msn.is_active:
            s.spikedetector_msn_channels = [nest.Create("spike_detector", params={"withgid": True, "withtime": True}) for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.msn_channels, s.spikedetector_msn_channels):
                nest.Connect(from_channel, to_detector)
            s.spikedetector_msn_d1_channels = [nest.Create("spike_detector", params={"withgid": True, "withtime": True}) for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.msn_d1_channels, s.spikedetector_msn_d1_channels):
                nest.Connect(from_channel, to_detector)
            s.spikedetector_msn_d2_channels = [nest.Create("spike_detector", params={"withgid": True, "withtime": True}) for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.msn_d2_channels, s.spikedetector_msn_d2_channels):
                nest.Connect(from_channel, to_detector)
            nest.Connect(s.msn, s.spikedetector_global)
            s.spikedetector_msn = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
            nest.Connect(s.msn, s.spikedetector_msn)

            if type(s.params.net.msn.is_active) is bool:
                s.multimeter_msn = nest.Create("multimeter")
                nest.SetStatus(s.multimeter_msn, {"withtime":True, "record_from":["V_m", "U_m", "AMPA", "NMDA", "GABA_A"], 'interval': s.params.sim.step})
                nest.Connect(s.multimeter_msn, [s.msn[1]])

        if s.params.net.snr.is_active:
            s.spikedetector_snr_channels = [nest.Create("spike_detector", params={"withgid": True, "withtime": True}) for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.snr_channels, s.spikedetector_snr_channels):
                nest.Connect(from_channel, to_detector)
            nest.Connect(s.snr, s.spikedetector_global)
            s.spikedetector_snr = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
            nest.Connect(s.snr, s.spikedetector_snr)

        if s.params.net.stn.is_active:
            s.spikedetector_stn_channels = [nest.Create("spike_detector", params={"withgid": True, "withtime": True}) for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.stn_channels, s.spikedetector_stn_channels):
                nest.Connect(from_channel, to_detector)
            nest.Connect(s.stn, s.spikedetector_global)
            s.spikedetector_stn = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
            nest.Connect(s.stn, s.spikedetector_stn)

        if s.params.net.gpe.is_active:
            s.spikedetector_gpe_channels = [nest.Create("spike_detector", params={"withgid": True, "withtime": True}) for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.gpe_channels, s.spikedetector_gpe_channels):
                nest.Connect(from_channel, to_detector)
            nest.Connect(s.gpe, s.spikedetector_global)
            s.spikedetector_gpe = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
            nest.Connect(s.gpe, s.spikedetector_gpe)


    def plot_activity(s, nucleus_name, channels=True, frequency=True, time_range=None, time_bin=50, n_bins=None, y_range=None):
        nucleus_name = nucleus_name.lower()
        s.__getattribute__(nucleus_name)

        if type(channels) is bool and channels == False:
            detector = [d for d in s.__getattribute__('spikedetector_' + nucleus_name)]
        elif (type(channels) is bool and channels == True):
            detector = [d[0] for d in s.__getattribute__('spikedetector_' + nucleus_name + '_channels')]
            channels = [i for i, _ in enumerate(detector)]
        elif type(channels) is int:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if channels == idx]
            channels = [channels]
        elif type(channels) is list:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if idx in channels]
        else:
            raise ValueError('Parameter `channel` must be bool, int or list of ints.')

        nucleus_spike_times = [d['times'] for d in nest.GetStatus(detector, 'events')]
        n_senders_channel = [len(nest.GetConnections(target=[d])) for d in detector]

        time_range = (0, nest.GetKernelStatus('time')) if time_range == None else time_range

        if n_bins:
            time_bin = (time_range[1] - time_range[0]) / n_bins
        else:
            n_bins = int((time_range[1] - time_range[0]) / time_bin)

        if y_range: plt.ylim(y_range)

        if frequency:
            for idx, (ch, n_senders) in enumerate(zip(nucleus_spike_times, n_senders_channel)):
                interval, count, __wtf = plt.hist(
                    ch,
                    bins=n_bins,
                    range=time_range,
                    histtype='step',
                    weights=np.ones_like(ch) / (n_senders * (time_bin / 1000.0)),
                    label='ch{}'.format(channels[idx]) if type(channels) is list else None
                )
            plt.ylabel('Hz')

        else:
            for idx, ch in enumerate(nucleus_spike_times):
                interval, count, __wtf = plt.hist(
                    ch,
                    bins=n_bins,
                    range=time_range,
                    histtype='step',
                    label='ch{}'.format(channels[idx]) if type(channels) is list else None
                )
            plt.ylabel('Spikes inside bin')
        plt.xlabel('Time (ms)\nBin size is {} ms'.format(int(time_bin)))
        if type(channels) is list: plt.legend(loc='upper left')


    def get_activity(s, nucleus_name, channels=True, frequency=True, time_range=None, time_bin=50, n_bins=None):
        nucleus_name = nucleus_name.lower()
        s.__getattribute__(nucleus_name)

        if type(channels) is bool and channels == False:
            detector = s.__getattribute__('spikedetector_' + nucleus_name)
        elif (type(channels) is bool and channels == True):
            detector = [d[0] for d in s.__getattribute__('spikedetector_' + nucleus_name + '_channels')]
        elif type(channels) is int:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if channels == idx]
        elif type(channels) is list:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if idx in channels]
        else:
            raise ValueError('Parameter `channel` must be bool, int or list of ints.')

        nucleus_spike_times = [d['times'] for d in nest.GetStatus(detector, 'events')]
        n_senders_channel = [len(nest.GetConnections(target=[d])) for d in detector]

        time_range = (0, nest.GetKernelStatus('time')) if time_range == None else time_range

        if n_bins:
            time_bin = (time_range[1] - time_range[0]) / n_bins
        else:
            n_bins = int((time_range[1] - time_range[0]) / time_bin)

        count_channel = []

        if frequency:
            for ch, n_senders in zip(nucleus_spike_times, n_senders_channel):
                count, interval = np.histogram(
                    ch,
                    bins=n_bins,
                    range=time_range,
                    weights=np.ones_like(ch) / (n_senders * (time_bin / 1000.0))
                )
                count_channel.append(count)
                intervals = interval

        else:
            for ch in nucleus_spike_times:
                count, interval = np.histogram(
                    ch,
                    bins=n_bins,
                    range=time_range,
                )
                count_channel.append(count)
                intervals = interval

        return intervals, count_channel


    def plot_raster(s, nucleus_name, channels=False, time_range=None, time_bin=50, n_bins=None, y_range=None):
        nucleus_name = nucleus_name.lower()
        s.__getattribute__(nucleus_name)

        if type(channels) is bool and channels == False:
            detectors = [s.__getattribute__('spikedetector_' + nucleus_name)]
        elif (type(channels) is bool and channels == True):
            detectors = [d[0] for d in s.__getattribute__('spikedetector_' + nucleus_name + '_channels')]
        elif type(channels) is int:
            detectors = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if channels == idx]
        elif type(channels) is list:
            detectors = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if idx in channels]
        else:
            raise ValueError('Parameter `channel` must be bool, int or list of ints.')

        time_range = (0, nest.GetKernelStatus('time')) if time_range == None else time_range

        if n_bins:
            time_bin = (time_range[1] - time_range[0]) / n_bins
        else:
            n_bins = int((time_range[1] - time_range[0]) / time_bin)

        for d in detectors:
            raster_plot.from_device(d, hist_binwidth=time_bin, ylim=y_range)


    def _compute_signals(s, nucleus='msn', time_range=None, mode='gaussian_convolution'):
        time_range = (0, 2500) if time_range is None else time_range

        if mode == 'gaussian_convolution':
            times, spikes_channels = s.get_activity(nucleus, time_bin=1)
            win = signal.gaussian(M=50, std=7.5)
            sum_win = sum(win)
            s._signals[nucleus] = [signal.convolve(s, win, mode='same')[time_range[0]:time_range[1]]/sum_win for s in spikes_channels]


    def signal(s, nucleus, channels=True, *args, **kwargs):

        try:
            s.__getattribute__('_signals')
        except AttributeError:
            s._signals = {}

        if not nucleus in s._signals:
            s._compute_signals(nucleus, *args, **kwargs)

        if type(channels) is bool:
            if channels==True:
                channels = list(range(len(s.params.net.p_channels)))
            else:
                return None #Why do this? Maybe to make sure that values are precomputed?
        elif type(channels) is int:
            channels = [channels]

        return itemgetter(*channels)(s._signals[nucleus])
