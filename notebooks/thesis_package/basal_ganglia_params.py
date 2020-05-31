import numpy as np
from math import floor, ceil

from .utils import dotdict, split_in_groups


# Modifiable params

is_verbose = [False, True][0]

is_msn_active = [False, True][1]
is_snr_active = [False, True][1]
is_stn_active = [False, True, 10.0][1]
is_stn_rb_active =   [False, True][1]
is_stn_llrs_active = [False, True][1]
is_stn_nr_active =   [False, True][1]
is_gpe_active = [False, True, 30.0][1]
is_gpe_a_active = [False, True][1]
is_gpe_b_active = [False, True][1]
is_gpe_c_active = [False, True][1]


# Definition of parameters

params = dotdict()


# SIMULATION params ############################################################

params.sim = dotdict()
params.sim.step = 0.01 #Humphries 2009
params.sim.n_threads = 6
params.sim.verbose = is_verbose


# NETWORK params ###############################################################

params.net = dotdict()
params.net.p_channels = [0.4, 0.4, 0.2] #Tomkins 2014, 3 canales de 40%, 40% y 20% de neuronas, pag 5
params.net.ctx_afferents_per_neuron = 250.0 #Humphries 2010

params.net.n_msn_per_mm3 = 84900 #Humphries 2010, pag 2
params.net.box_side_mm = 0.3 #Tomkins 2014, pag 4
params.net.box_side_um = params.net.box_side_mm * 1000.0
params.net.box_volume_mm3 = params.net.box_side_mm ** 3

params.net.msn = dotdict()
params.net.msn.n = int(params.net.box_volume_mm3 * params.net.n_msn_per_mm3)
params.net.msn.d1_n = int(params.net.msn.n * 0.5)
params.net.msn.d2_n = int(params.net.msn.n * 0.5)
params.net.msn.is_active = is_msn_active
params.net.msn.reduction_factor = 1.0 / 10.0 #Cuando se usa una población poisson en lugar del msn, reducimos el número de neuronas para acelerar los calculos.

params.net.stn = dotdict()
params.net.stn.is_active = is_stn_active
params.net.stn.n = ceil(13600/300) #Fountas

params.net.stn_rb = dotdict()
params.net.stn_rb.is_active = is_stn_rb_active
params.net.stn_rb.n = ceil(params.net.stn.n * 0.6)
params.net.stn_llrs = dotdict()
params.net.stn_llrs.is_active = is_stn_llrs_active
params.net.stn_llrs.n = ceil(params.net.stn.n * 0.25)
params.net.stn_nr = dotdict()
params.net.stn_nr.is_active = is_stn_nr_active
params.net.stn_nr.n = ceil(params.net.stn.n * 0.15)

params.net.gpe = dotdict()
params.net.gpe.is_active = is_gpe_active
params.net.gpe.n = ceil(46000.0/300.0) #Fountas
params.net.gpe.has_local_collaterals = True

params.net.gpe_a = dotdict()
params.net.gpe_a.is_active = is_gpe_a_active
params.net.gpe_a.n = ceil(params.net.gpe.n * 0.0405)
params.net.gpe_b = dotdict()
params.net.gpe_b.is_active = is_gpe_b_active
params.net.gpe_b.n = ceil(params.net.gpe.n * 0.85)
params.net.gpe_c = dotdict()
params.net.gpe_c.is_active = is_gpe_c_active
params.net.gpe_c.n = ceil(params.net.gpe.n * 0.1095)

params.net.snr = dotdict()
params.net.snr.n = int(26300/300) #Fountas, Oorschot 1996
params.net.snr.is_active = is_snr_active


# NEURONS params ###############################################################

params.neuron = dotdict()
params.neuron.msn_d1 = dotdict()
params.neuron.msn_d1.C_m = 15.0
params.neuron.msn_d1.E_rev = (0.0, 0.0, -60.0) #TODO: falta fuente
params.neuron.msn_d1.tau_syn = (6.0, 160.0, 11.0) #(6.0, 160.0, 4.0) #GABA_A Modificado según Götz1997
params.neuron.msn_d1.Kplus = 1.0
params.neuron.msn_d1.V_T = -30.0
params.neuron.msn_d1.vp = -80.0 #-76.5
params.neuron.msn_d1.V_th = 40.0
params.neuron.msn_d1.a = 0.01
params.neuron.msn_d1.b = -20.0
params.neuron.msn_d1.c = -55.0
params.neuron.msn_d1.d = 91.0
params.neuron.msn_d1.y1 = 0.3
params.neuron.msn_d1.y2 = 0.3
params.neuron.msn_d1.alpha_1 = 6.3 #Humphries
params.neuron.msn_d1.alpha_2 = 0.215 #Humphries
params.neuron.msn_d1.alpha = 0.0 #D1 type
params.neuron.msn_d1.c_1 = 0.0289 #D1 type
params.neuron.msn_d1.c_2 = 0.331 #D1 type


params.neuron.msn_d2 = dotdict()
params.neuron.msn_d2.C_m = 15.0
params.neuron.msn_d2.E_rev = (0.0, 0.0, -60.0) #TODO: falta fuente
params.neuron.msn_d2.tau_syn = (6.0, 160.0, 11.0) #(6.0, 160.0, 4.0) #GABA_A Modificado según Götz1997
params.neuron.msn_d2.Kplus = 1.0
params.neuron.msn_d2.V_T = -30.0
params.neuron.msn_d2.vp = -80.0 #-76.5
params.neuron.msn_d2.V_th = 40.0
params.neuron.msn_d2.a = 0.01
params.neuron.msn_d2.b = -20.0
params.neuron.msn_d2.c = -55.0
params.neuron.msn_d2.d = 91.0
params.neuron.msn_d2.y1 = 0.3
params.neuron.msn_d2.y2 = 0.3
params.neuron.msn_d2.alpha_1 = 6.3 #Humphries
params.neuron.msn_d2.alpha_2 = 0.215 #Humphries
params.neuron.msn_d2.alpha = 0.032 #D2 type
params.neuron.msn_d2.c_1 = 0.0 #D2 type
params.neuron.msn_d2.c_2 = 0.0 #D2 type

params.neuron.snr = dotdict()
params.neuron.snr.a = 0.113 #Fountas #TODO: revisar si la fuente es fiable
params.neuron.snr.b = 11.057 #Fountas
params.neuron.snr.c = -62.7 #Fountas
params.neuron.snr.d = 138.4 #Fountas
params.neuron.snr.V_m = -64.58 # -55.8, #Fountas
params.neuron.snr.V_T = -51.8 # -55.2, #Fountas
params.neuron.snr.V_th = 9.8 # 20.0, #Fountas
params.neuron.snr.C_m = 172.1 # 80.0, #Fountas #200.0 C_sim
params.neuron.snr.Kplus = 0.7836 # 1.731, #Fountas
params.neuron.snr.tau_syn = [2.0, 100.0, 5.2, 2.1, 3.0] #Fountas #TODO: Revisar tau canales
params.neuron.snr.E_rev = [0.0, 0.0, -80.0, -80.0, -80.0] #Fountas #TODO: Revisar E_rev canales
params.neuron.snr.I_e = 690.4457293525152 #Búsqueda local siguiendo procedimiento de Fountas #235.0  # Vivo:235.0, Vitro:150.0 #Fountas


params.neuron.gpe_a = dotdict()
params.neuron.gpe_b = dotdict()
params.neuron.gpe_c = dotdict()

params.neuron.gpe_a.tau_syn = [2.0, 100.0, 6.0, 5.0]
params.neuron.gpe_a.E_rev = [0.0, 0.0, -65.0, -65.0]
params.neuron.gpe_a.V_m = -50.7
params.neuron.gpe_a.V_T = -42.0
params.neuron.gpe_a.V_th = 38.0
params.neuron.gpe_a.C_m = 55.0
params.neuron.gpe_a.a = 0.29
params.neuron.gpe_a.b = 4.26
params.neuron.gpe_a.c = -57.4
params.neuron.gpe_a.d = 110.0
params.neuron.gpe_a.Kplus = 0.06
params.neuron.gpe_a.I_e = 167.0  # Vivo:167.0, Vitro:107.0

params.neuron.gpe_b.tau_syn = [2.0, 100.0, 6.0, 5.0]
params.neuron.gpe_b.E_rev = [0.0, 0.0, -65.0, -65.0]
params.neuron.gpe_b.V_m = -53.0
params.neuron.gpe_b.V_T = -44.0
params.neuron.gpe_b.V_th = 25.0
params.neuron.gpe_b.C_m = 68.0
params.neuron.gpe_b.a = 0.0045
params.neuron.gpe_b.b = 3.895
params.neuron.gpe_b.c = -58.36
params.neuron.gpe_b.d = 0.353
params.neuron.gpe_b.Kplus = 0.943
params.neuron.gpe_b.I_e = 64.0 # Vivo:64.0, Vitro:52

params.neuron.gpe_c.tau_syn = [2.0, 100.0, 6.0, 5.0]
params.neuron.gpe_c.E_rev = [0.0, 0.0, -65.0, -65.0]
params.neuron.gpe_c.V_m = -54.0
params.neuron.gpe_c.V_T = -43.0
params.neuron.gpe_c.V_th = 34.5
params.neuron.gpe_c.C_m = 57.0
params.neuron.gpe_c.a = 0.42
params.neuron.gpe_c.b = 7.0
params.neuron.gpe_c.c = -52.0
params.neuron.gpe_c.d = 166.0
params.neuron.gpe_c.Kplus = 0.099
params.neuron.gpe_c.I_e = 237.5  # Vivo:237.5, Vitro:187.5


params.neuron.stn_rb = dotdict()
params.neuron.stn_llrs = dotdict()
params.neuron.stn_nr = dotdict()

params.neuron.stn_rb.tau_syn = [2.0, 100.0, 8.0]
params.neuron.stn_rb.E_rev = [0.0, 0.0, -84.0]
params.neuron.stn_rb.V_m = -56.2
params.neuron.stn_rb.V_T = -41.4
params.neuron.stn_rb.V_th = 15.4
params.neuron.stn_rb.C_m = 23.0
params.neuron.stn_rb.a = 0.021
params.neuron.stn_rb.b = 4.0
params.neuron.stn_rb.c = -47.7
params.neuron.stn_rb.d = 17.1
params.neuron.stn_rb.alpha = 0.123
params.neuron.stn_rb.beta = 0.015
params.neuron.stn_rb.delta_P = -68.4
params.neuron.stn_rb.V_epsp = -60.0
params.neuron.stn_rb.Kplus = 0.439
params.neuron.stn_rb.Wmax = 0.1
params.neuron.stn_rb.Wmin = 0.0
params.neuron.stn_rb.type_id = False
params.neuron.stn_rb.I_e = 56.1  # Vivo:56.1, Vitro:56.1

params.neuron.stn_llrs.tau_syn = [2.0, 100.0, 8.0]
params.neuron.stn_llrs.E_rev = [0.0, 0.0, -84.0]
params.neuron.stn_llrs.V_m = -56.2
params.neuron.stn_llrs.V_T = -50.0
params.neuron.stn_llrs.V_th = 15.4
params.neuron.stn_llrs.C_m = 40.0
params.neuron.stn_llrs.a = 0.05
params.neuron.stn_llrs.b = 0.2
params.neuron.stn_llrs.c = -60.0
params.neuron.stn_llrs.d = 1.0
params.neuron.stn_llrs.alpha = 0.001
params.neuron.stn_llrs.beta = 0.3
params.neuron.stn_llrs.delta_P = 10.0
params.neuron.stn_llrs.V_epsp = -60.0
params.neuron.stn_llrs.Kplus = 0.3
params.neuron.stn_llrs.Wmax = 0.01
params.neuron.stn_llrs.Wmin = 0.0
params.neuron.stn_llrs.type_id = False
params.neuron.stn_llrs.I_e = 8.0  # Vivo:8.0, Vitro:25.0

params.neuron.stn_nr.tau_syn = [2.0, 100.0, 8.0]
params.neuron.stn_nr.E_rev = [0.0, 0.0, -84.0]
params.neuron.stn_nr.V_m = -58.5
params.neuron.stn_nr.V_T = -43.75
params.neuron.stn_nr.V_th = 15.4
params.neuron.stn_nr.C_m = 30.0
params.neuron.stn_nr.a = 0.44
params.neuron.stn_nr.b = -1.35
params.neuron.stn_nr.c = -52.34
params.neuron.stn_nr.d = 17.65
params.neuron.stn_nr.alpha = 0.32
params.neuron.stn_nr.beta = 3.13
params.neuron.stn_nr.delta_P = 92.0
params.neuron.stn_nr.V_epsp = -43.2
params.neuron.stn_nr.Kplus = 0.105
params.neuron.stn_nr.Wmax = 0.001
params.neuron.stn_nr.Wmin = 1.0
params.neuron.stn_nr.type_id = True
params.neuron.stn_nr.I_e = -18.0 # Vivo:-18.0, Vitro:-1.0


# SYNAPSES params ###########################################################

params.syn = dotdict()

params.syn.ctx_msn = dotdict()
params.syn.ctx_msn.conn_spec = dotdict()
params.syn.ctx_msn.conn_spec.rule = 'one_to_one'
params.syn.ctx_msn.syn_spec_ampa = dotdict()
params.syn.ctx_msn.syn_spec_ampa.model = 'static_synapse'
params.syn.ctx_msn.syn_spec_ampa.delay = 10.0 #A mano, no importa
params.syn.ctx_msn.syn_spec_ampa.weight = 6.1 / params.neuron.msn_d1.tau_syn[0] #Humphries 2009b, pag 1178, tabla 2
params.syn.ctx_msn.syn_spec_ampa.receptor_type = 1
params.syn.ctx_msn.syn_spec_nmda = dotdict()
params.syn.ctx_msn.syn_spec_nmda.model = 'static_synapse'
params.syn.ctx_msn.syn_spec_nmda.delay = 10.0 #A mano, no importa
params.syn.ctx_msn.syn_spec_nmda.weight = 3.05 / params.neuron.msn_d1.tau_syn[1] #Humphries 2009b, pag 1178, tabla 2
params.syn.ctx_msn.syn_spec_nmda.receptor_type = 2


params.syn.msn_msn = dotdict()
params.syn.msn_msn.conn_spec_random = dotdict()
params.syn.msn_msn.conn_spec_random.rule = 'pairwise_bernoulli'
params.syn.msn_msn.conn_spec_random.p = 728.0 / params.net.msn.n #Tomkins 2014, pag 5
params.syn.msn_msn.conn_spec_random.autapses = False
params.syn.msn_msn.syn_spec_random = dotdict()
params.syn.msn_msn.syn_spec_random.model = 'static_synapse'
params.syn.msn_msn.syn_spec_random.delay = {"distribution": "uniform", "low": 1.0, "high": 2.0} #Humphries, código
params.syn.msn_msn.syn_spec_random.weight = 0.25 #A mano
params.syn.msn_msn.syn_spec_random.receptor_type = 3

params.syn.msn_snr = dotdict()
params.syn.msn_snr.conn_spec = dotdict()
params.syn.msn_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.msn_snr.conn_spec.p = 0.033 #Fountas
params.syn.msn_snr.syn_spec_static = dotdict()
params.syn.msn_snr.syn_spec_static.model = 'static_synapse'
params.syn.msn_snr.syn_spec_static.delay = 1.0 #Fountas
params.syn.msn_snr.syn_spec_static.weight = 57.654700925145136 # Búsqueda local siguiendo procedimiento de Fountas #4.5 #Fountas
params.syn.msn_snr.syn_spec_static.receptor_type = 3

params.syn.msn_snr.syn_spec_plastic = dotdict()
params.syn.msn_snr.syn_spec_plastic.model = 'tsodyks2_synapse'
params.syn.msn_snr.syn_spec_plastic.U = 0.0192 #Fountas
params.syn.msn_snr.syn_spec_plastic.u = 1.0#0.0192 #Fountas
params.syn.msn_snr.syn_spec_plastic.x = 0.0192#1.0
params.syn.msn_snr.syn_spec_plastic.tau_rec = 623.0 #Fountas
params.syn.msn_snr.syn_spec_plastic.tau_fac = 559.0 #Fountas
params.syn.msn_snr.syn_spec_plastic.delay = 1.0 #Fountas
params.syn.msn_snr.syn_spec_plastic.weight = 156.3 #Fountas
params.syn.msn_snr.syn_spec_plastic.receptor_type = 3

params.syn.snr_snr = dotdict()
params.syn.snr_snr.conn_spec = dotdict()
params.syn.snr_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.snr_snr.conn_spec.p = 0.1 #Fountas
params.syn.snr_snr.syn_spec = dotdict()
params.syn.snr_snr.syn_spec.model = 'static_synapse'
params.syn.snr_snr.syn_spec.delay = 1.0 #Fountas
params.syn.snr_snr.syn_spec.weight = 0.32539679447089986 # Búsqueda local siguiendo procedimiento de Fountas #0.2 #Fountas
params.syn.snr_snr.syn_spec.receptor_type = 5


params.syn.ctx_stn = dotdict()
params.syn.ctx_stn.conn_spec = dotdict()
params.syn.ctx_stn.conn_spec.rule = 'one_to_one'

params.syn.ctx_stn.syn_spec_ampa = dotdict()
params.syn.ctx_stn.syn_spec_ampa.model = 'static_synapse'
params.syn.ctx_stn.syn_spec_ampa.delay = 2.5
params.syn.ctx_stn.syn_spec_ampa.weight = 0.0215 #Ajuste manual siguiendo procedimiento de Fountas #0.01 #Fountas
params.syn.ctx_stn.syn_spec_ampa.receptor_type = 1

params.syn.ctx_stn.syn_spec_nmda = params.syn.ctx_stn.syn_spec_ampa.copy()
params.syn.ctx_stn.syn_spec_nmda.weight *= 0.6 #Fountas
params.syn.ctx_stn.syn_spec_nmda.receptor_type = 2


params.syn.stn_gpe = dotdict()
params.syn.stn_gpe.conn_spec = dotdict()
params.syn.stn_gpe.conn_spec.rule = 'pairwise_bernoulli'
params.syn.stn_gpe.conn_spec.p = 0.3 #Fountas

params.syn.stn_gpe.syn_spec_ampa = dotdict()
params.syn.stn_gpe.syn_spec_ampa.model = 'static_synapse'
params.syn.stn_gpe.syn_spec_ampa.delay = 2.0
params.syn.stn_gpe.syn_spec_ampa.weight = 0.3 #Ajuste manual siguiendo procedimiento de Fountas  #1.447 #Fountas
params.syn.stn_gpe.syn_spec_ampa.receptor_type = 1

params.syn.stn_gpe.syn_spec_nmda = params.syn.stn_gpe.syn_spec_ampa.copy()
params.syn.stn_gpe.syn_spec_nmda.weight *= 0.36
params.syn.stn_gpe.syn_spec_nmda.receptor_type = 2


params.syn.gpe_stn = dotdict()
params.syn.gpe_stn.conn_spec = dotdict()
params.syn.gpe_stn.conn_spec.rule = 'pairwise_bernoulli'
params.syn.gpe_stn.conn_spec.p = 0.1

params.syn.gpe_stn.syn_spec = dotdict()
params.syn.gpe_stn.syn_spec.model = 'static_synapse'
params.syn.gpe_stn.syn_spec.delay = 4.0
params.syn.gpe_stn.syn_spec.weight = 0.518
params.syn.gpe_stn.syn_spec.receptor_type = 3


params.syn.gpe_gpe = dotdict()
params.syn.gpe_gpe.conn_spec = dotdict()
params.syn.gpe_gpe.conn_spec.rule = 'pairwise_bernoulli'
params.syn.gpe_gpe.conn_spec.p = 0.1

params.syn.gpe_gpe.syn_spec = dotdict()
params.syn.gpe_gpe.syn_spec.model = 'static_synapse'
params.syn.gpe_gpe.syn_spec.delay = 1.0
params.syn.gpe_gpe.syn_spec.weight = 0.765
params.syn.gpe_gpe.syn_spec.receptor_type = 4


params.syn.msn_gpe = dotdict()
params.syn.msn_gpe.conn_spec = dotdict()
params.syn.msn_gpe.conn_spec.rule = 'pairwise_bernoulli'
params.syn.msn_gpe.conn_spec.p = 0.033 #Fountas

params.syn.msn_gpe.syn_spec_plastic = dotdict()
params.syn.msn_gpe.syn_spec_plastic.model = 'tsodyks2_synapse'
params.syn.msn_gpe.syn_spec_plastic.U = 0.24 #Fountas
params.syn.msn_gpe.syn_spec_plastic.u = 1.0#0.24 #Fountas
params.syn.msn_gpe.syn_spec_plastic.x = 0.24#1.0
params.syn.msn_gpe.syn_spec_plastic.tau_rec = 11.0 #Fountas
params.syn.msn_gpe.syn_spec_plastic.tau_fac = 73.0 #Fountas
params.syn.msn_gpe.syn_spec_plastic.delay = 5.0 #Fountas
params.syn.msn_gpe.syn_spec_plastic.weight = 21.6 #Fountas
params.syn.msn_gpe.syn_spec_plastic.receptor_type = 3

params.syn.msn_gpe.syn_spec_static = dotdict()
params.syn.msn_gpe.syn_spec_static.model = 'static_synapse'
params.syn.msn_gpe.syn_spec_static.delay = 5.0 #Fountas
params.syn.msn_gpe.syn_spec_static.weight = 10.0 # Ajuste manual siguiendo procedimiento de Fountas #5.435 #Fountas
params.syn.msn_gpe.syn_spec_static.receptor_type = 3


params.syn.stn_snr = dotdict()
params.syn.stn_snr.conn_spec = dotdict()
params.syn.stn_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.stn_snr.conn_spec.p = 0.3

params.syn.stn_snr.syn_spec_plastic_ampa = dotdict()
params.syn.stn_snr.syn_spec_plastic_ampa.model = 'tsodyks2_synapse'
params.syn.stn_snr.syn_spec_plastic_ampa.U = 0.35
params.syn.stn_snr.syn_spec_plastic_ampa.u = 1.0#0.35
params.syn.stn_snr.syn_spec_plastic_ampa.x = 0.35#1.0
params.syn.stn_snr.syn_spec_plastic_ampa.tau_rec = 800.0
params.syn.stn_snr.syn_spec_plastic_ampa.tau_fac = 0.0
params.syn.stn_snr.syn_spec_plastic_ampa.delay = 1.5
params.syn.stn_snr.syn_spec_plastic_ampa.weight = 49.5
params.syn.stn_snr.syn_spec_plastic_ampa.receptor_type = 1

params.syn.stn_snr.syn_spec_plastic_nmda = dotdict()
params.syn.stn_snr.syn_spec_plastic_nmda.model = 'tsodyks2_synapse'
params.syn.stn_snr.syn_spec_plastic_nmda.U = 0.35
params.syn.stn_snr.syn_spec_plastic_nmda.u = 1.0#0.35
params.syn.stn_snr.syn_spec_plastic_nmda.x = 0.35#1.0
params.syn.stn_snr.syn_spec_plastic_nmda.tau_rec = 800.0
params.syn.stn_snr.syn_spec_plastic_nmda.tau_fac = 0.0
params.syn.stn_snr.syn_spec_plastic_nmda.delay = 1.5
params.syn.stn_snr.syn_spec_plastic_nmda.weight = 49.5 * 0.2
params.syn.stn_snr.syn_spec_plastic_nmda.receptor_type = 2

params.syn.stn_snr.syn_spec_static_ampa = dotdict()
params.syn.stn_snr.syn_spec_static_ampa.model = 'static_synapse'
params.syn.stn_snr.syn_spec_static_ampa.delay = 1.5
params.syn.stn_snr.syn_spec_static_ampa.weight = 3.39200118729286 #Búsqueda local siguiendo procedimiento de Fountas #20.8 Fountas
params.syn.stn_snr.syn_spec_static_ampa.receptor_type = 1

params.syn.stn_snr.syn_spec_static_nmda = dotdict()
params.syn.stn_snr.syn_spec_static_nmda.model = 'static_synapse'
params.syn.stn_snr.syn_spec_static_nmda.delay = 1.5
params.syn.stn_snr.syn_spec_static_nmda.weight = 0.2 * params.syn.stn_snr.syn_spec_static_ampa.weight
params.syn.stn_snr.syn_spec_static_nmda.receptor_type = 2


params.syn.gpe_snr = dotdict()
params.syn.gpe_snr.conn_spec = dotdict()
params.syn.gpe_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.gpe_snr.conn_spec.p = 0.1066

params.syn.gpe_snr.syn_spec_static = dotdict()
params.syn.gpe_snr.syn_spec_static.model = 'static_synapse'
params.syn.gpe_snr.syn_spec_static.delay = 3.0
params.syn.gpe_snr.syn_spec_static.weight = 59.67154181984671 #Búsqueda local siguiendo procedimiento de Fountas #93.0 #Fountas
params.syn.gpe_snr.syn_spec_static.receptor_type = 4

params.syn.gpe_snr.syn_spec_plastic = dotdict()
params.syn.gpe_snr.syn_spec_plastic.model = 'tsodyks2_synapse'
params.syn.gpe_snr.syn_spec_plastic.U = 0.196
params.syn.gpe_snr.syn_spec_plastic.u = 1.0#0.196
params.syn.gpe_snr.syn_spec_plastic.x = 0.196#1.0
params.syn.gpe_snr.syn_spec_plastic.tau_rec = 969.0
params.syn.gpe_snr.syn_spec_plastic.tau_fac = 0.0
params.syn.gpe_snr.syn_spec_plastic.delay = 3.0
params.syn.gpe_snr.syn_spec_plastic.weight = 603.9
params.syn.gpe_snr.syn_spec_plastic.receptor_type = 4
