# Created by Kangxin Hu at Beihang University, on Oct.2023.

from __future__ import division
from neuron import h, nrn, gui

h.load_file('stdrun.hoc')
h.nrn_load_dll('/mod/nrnmech.dll')

import matplotlib.pyplot as plt
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
soma.gshl1bar_shl_pas = 4.4627
soma.insert('SHK1_pas')
soma.gshk1bar_SHK1_pas = 1.10337
soma.insert('KVS1_pas')
soma.gkvs1bar_KVS1_pas = 0.78251
soma.insert('KQT3_pas')   # question
soma.gkqt3bar_KQT3_pas = 0.00405
soma.insert('EGL2_pas')
soma.gegl2bar_EGL2_pas = 0.11278
soma.insert('EGL36_pas')
soma.gegl36bar_EGL36_pas = 0
soma.insert('IRK_pas')
soma.girkbar_IRK_pas = 0.17725

soma.insert('EGL19_ca')  # add calcium channel
soma.gegl19bar_EGL19_ca = 1.46841
soma.insert('UNC2_ca')
soma.gunc2bar_UNC2_ca = 0.98413
soma.insert('CCA1_ca')
soma.gcca1bar_CCA1_ca = 7.48234

soma.insert('SLO1_EGL19_cak')  # add potassium channel
soma.gslo1_egl19bar_SLO1_EGL19_cak = 0.18972
h.setpointer(soma(0.5)._ref_minf_EGL19_ca, 'minf_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_mtau_EGL19_ca, 'mtau_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_m_EGL19_ca, 'm_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_h_EGL19_ca, 'h_egl19', soma(0.5).SLO1_EGL19_cak)
soma.insert('SLO1_UNC2_cak')
soma.gslo1_unc2bar_SLO1_UNC2_cak = 0.18972
h.setpointer(soma(0.5)._ref_minf_UNC2_ca, 'minf_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_mtau_UNC2_ca, 'mtau_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_m_UNC2_ca, 'm_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_h_UNC2_ca, 'h_unc2', soma(0.5).SLO1_UNC2_cak)
soma.insert('SLO2_EGL19_cak')
soma.gslo2_egl19bar_SLO2_EGL19_cak = 0.12285
h.setpointer(soma(0.5)._ref_minf_EGL19_ca, 'minf_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_mtau_EGL19_ca, 'mtau_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_m_EGL19_ca, 'm_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_h_EGL19_ca, 'h_egl19', soma(0.5).SLO2_EGL19_cak)
soma.insert('SLO2_UNC2_cak')
soma.gslo2_unc2bar_SLO2_UNC2_cak = 0.12285
h.setpointer(soma(0.5)._ref_minf_UNC2_ca, 'minf_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_mtau_UNC2_ca, 'mtau_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_m_UNC2_ca, 'm_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_h_UNC2_ca, 'h_unc2', soma(0.5).SLO2_UNC2_cak)

soma.insert('Intercell_CaSK')  # add potassium channel
soma.insert('KCNL_CaSK')
soma.gkcnlbar_KCNL_CaSK = 0.06

soma.insert('leak')  # add sodium and leak channel
soma.glbar_leak = 0.35
soma(0.5).el_leak = -80  # -80 for RMD
#################################################### Set up experiment #################################################
# stim = h.VClamp(soma(0.5))  # add a voltage clamp at the middle of the soma
# stim.dur[0] = 200
# stim.amp[0] = 110

stim = h.IClamp(soma(0.5))  # add a current clamp on the middle of the neuron[0].soma #ASH
stim.delay = 10  # ms
stim.dur = 400  # ms
stim.amp = 150  # pA

soma_ik = h.Vector()
soma_ik.record(soma(0.5)._ref_ik)
soma_ica = h.Vector()
soma_ica.record(soma(0.5)._ref_ica)
soma_il = h.Vector()
soma_il.record(soma(0.5)._ref_il_leak)
soma_i = h.Vector()
soma_i.record(soma(0.5)._ref_i)
#
# stim1 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim1.delay = 0  # ms
# stim1.dur = 100  # ms
# stim1.amp = 0  # pA
#
# stim2 = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
# stim2.delay = 100  # ms
# stim2.dur = 500  # ms
# stim2.amp = 20  # pA
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
h.v_init = -70 # set starting voltage
h.tstop = 400  # set simulation time
h.run()  # run simulation
#################################################### Plot Results ######################################################
plt.figure(figsize=(8, 5))
# plt.plot(t, soma_ik, color='b', label='ik')
# plt.plot(t, soma_ica, color='g', label='ica')
# plt.plot(t, soma_ina, color='r', label='ina')
# plt.plot(t, soma_il, color='c', label='il')
plt.plot(t, soma_v, color='b', label='V')
plt.xlim(0, h.tstop)
plt.xlabel('Time (ms)', fontsize=15)
plt.ylabel('v (mv)', fontsize=15)
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.ioff()
plt.show()