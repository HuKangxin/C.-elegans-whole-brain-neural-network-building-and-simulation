# Created by Kangxin Hu at Beihang University, on Oct.2023.

from __future__ import division
from neuron import h, nrn, gui

import matplotlib.pyplot as plt
h.nrn_load_dll('/mod/nrnmech.dll')
plt.ion()
#################################################### Build Model #######################################################
h.celsius = 30
soma = h.Section(name="soma")
soma.L = 5  # um
soma.diam = 5  # um
# h("forall { nseg = int((L/(0.1*lambda_f(100))+0.9)/2)*2 + 1  }")
soma.nseg = 1
soma.cm = 3.1
# ##
# dend = h.Section(name="dend")
# dend.L = 20  # um
# dend.diam = 1  # um
# h("forall { nseg = int((L/(0.1*lambda_f(100))+0.9)/2)*2 + 1  }")
# print(dend.nseg)
# dend.cm = 3.1
# ##
# dend.connect(soma)
#################################################### Ion Channel #######################################################
soma.insert('shl_pas')  # add potassium channel
soma.gshl1bar_shl_pas = 2.9
soma.insert('SHK1_pas')
soma.gshk1bar_SHK1_pas = 0.1
soma.insert('KVS1_pas')
soma.gkvs1bar_KVS1_pas = 0.8
soma.insert('KQT3_pas')   # question
soma.gkqt3bar_KQT3_pas = 0.55
soma.insert('EGL2_pas')
soma.gegl2bar_EGL2_pas = 0.85
soma.insert('EGL36_pas')
soma.gegl36bar_EGL36_pas = 0
soma.insert('IRK_pas')
soma.girkbar_IRK_pas = 0.65

soma.insert('EGL19_ca')  # add calcium channel
soma.gegl19bar_EGL19_ca = 1.55
soma.insert('UNC2_ca')
soma.gunc2bar_UNC2_ca = 1
soma.insert('CCA1_ca')
soma.gcca1bar_CCA1_ca = 0.7

soma.insert('SLO1_EGL19_cak')  # add potassium channel
soma.gslo1_egl19bar_SLO1_EGL19_cak = 0.11
h.setpointer(soma(0.5)._ref_minf_EGL19_ca, 'minf_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_mtau_EGL19_ca, 'mtau_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_m_EGL19_ca, 'm_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_h_EGL19_ca, 'h_egl19', soma(0.5).SLO1_EGL19_cak)
soma.insert('SLO1_UNC2_cak')
soma.gslo1_unc2bar_SLO1_UNC2_cak = 0.11
h.setpointer(soma(0.5)._ref_minf_UNC2_ca, 'minf_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_mtau_UNC2_ca, 'mtau_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_m_UNC2_ca, 'm_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_h_UNC2_ca, 'h_unc2', soma(0.5).SLO1_UNC2_cak)
soma.insert('SLO2_EGL19_cak')
soma.gslo2_egl19bar_SLO2_EGL19_cak = 0.1
h.setpointer(soma(0.5)._ref_minf_EGL19_ca, 'minf_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_mtau_EGL19_ca, 'mtau_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_m_EGL19_ca, 'm_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_h_EGL19_ca, 'h_egl19', soma(0.5).SLO2_EGL19_cak)
soma.insert('SLO2_UNC2_cak')
soma.gslo2_unc2bar_SLO2_UNC2_cak = 0.1
h.setpointer(soma(0.5)._ref_minf_UNC2_ca, 'minf_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_mtau_UNC2_ca, 'mtau_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_m_UNC2_ca, 'm_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_h_UNC2_ca, 'h_unc2', soma(0.5).SLO2_UNC2_cak)

soma.insert('Intercell_CaSK')  # add potassium channel
soma.insert('KCNL_CaSK')
soma.gkcnlbar_KCNL_CaSK = 1.8

soma.insert('na_leak')  # add sodium and leak channel
soma.gnabar_na_leak = 0.06
soma.glbar_na_leak = 0.38
soma(0.5).el_na_leak = -84  # -80 for RMD
soma.ena = 55

stim3 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
stim3.delay = 10  # ms
stim3.dur = 50  # ms
stim3.amp = 40  # nA
'''
syn_ = h.Exc2Syn(soma(0.5)) # assumes synapse is in middle of dend
# NetStim that drives the first activation of syn
ns1 = h.NetStim(soma(0.5))
# statements that specify number and timing of ns1's events
# for example, to make ns1 generate one event at t = 10 ms
ns1.number = 1
ns1.start = 10 # time of event in ms
# NetCon that conveys ns1's events to syn
nc1 = h.NetCon(ns1, syn_)
nc1.delay=10
nc1.weight[0] = 5
'''
#################################################### Set up experiment #################################################
# stim = h.VClamp(soma(0.5))  # add a voltage clamp at the middle of the soma
# stim.dur[0] = 200
# stim.amp[0] = 110
soma_ik = h.Vector()
soma_ik.record(soma(0.5)._ref_ik)
soma_ica = h.Vector()
soma_ica.record(soma(0.5)._ref_ica)
soma_ina = h.Vector()
soma_ina.record(soma(0.5)._ref_ina)
soma_il = h.Vector()
soma_il.record(soma(0.5)._ref_il_na_leak)
soma_i = h.Vector()
soma_i.record(soma(0.5)._ref_itotal_na_leak)
#


# stim3 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim3.delay = 600  # ms
# stim3.dur = 400  # ms
# stim3.amp = 0  # pA
#
# stim4 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim4.delay = 1000  # ms
# stim4.dur = 500  # ms
# stim4.amp = 20  # pA
#
# stim5 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim5.delay = 1500  # ms
# stim5.dur = 500  # ms
# stim5.amp = 0  # pA
#
# stim6 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim6.delay = 2000  # ms
# stim6.dur = 500  # ms
# stim6.amp = 20  # pA
#
# stim7 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim7.delay = 2500  # ms
# stim7.dur = 500  # ms
# stim7.amp = 0  # pA
#
soma_v = h.Vector()
soma_v.record(soma(0.5)._ref_v)


# stim1 = h.NetStim() # Make a new stimulator
# syn_ = h.ExpSyn(dend(0.1))
# stim1.number = 10
# stim1.start = 9
# ncstim = h.NetCon(stim1, syn_)
# ncstim.delay = 10
# ncstim.weight[0] = 1 # NetCon weight is a vector.
#
print(soma.psection())
# print(h.units('cai'))

t = h.Vector()
t.record(h._ref_t)  # record time
h.v_init = -65 # set starting voltage
h.tstop = 300  # set simulation time
h.run()  # run simulation
#################################################### Plot Results ######################################################
plt.figure(figsize=(8, 5))
# plt.plot(t, soma_ik, color='b', label='ik')
# plt.plot(t, soma_ica, color='g', label='ica')
# plt.plot(t, soma_ina, color='r', label='ina')
plt.plot(t, soma_i, color='c', label='itotal')
plt.plot(t, soma_v, color='b', label='V')
plt.xlim(0, h.tstop)
plt.xlabel('Time (ms)', fontsize=15)
plt.ylabel('I (nA)', fontsize=15)
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.ioff()
plt.show()
