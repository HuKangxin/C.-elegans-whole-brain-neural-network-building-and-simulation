# Created by Kangxin Hu at Beihang University, on Oct.2023.

from __future__ import division
from neuron import h, nrn, gui

import matplotlib.pyplot as plt
plt.ion()
#################################################### Build Model #######################################################
h.celsius = 30
soma = h.Section(name="soma")
soma.L = 10  # um
soma.diam = 10  # um
# h("forall { nseg = int((L/(0.1*lambda_f(100))+0.9)/2)*2 + 1  }")
soma.nseg = 1
soma.cm = 6.1
soma.Ra = 1

h.load_file('stdrun.hoc')
h.nrn_load_dll('/mod/nrnmech.dll')

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
'''
soma.insert('shl_pas')  # add potassium channel
soma.gshl1bar_shl_pas = 4.52005
soma.insert('SHK1_pas')
soma.gshk1bar_SHK1_pas = 1.10234
soma.insert('KVS1_pas')
soma.gkvs1bar_KVS1_pas = 1.88709
soma.insert('KQT3_pas')   # question
soma.gkqt3bar_KQT3_pas = 0.2
soma.insert('EGL2_pas')
soma.gegl2bar_EGL2_pas = 0.02235
soma.insert('EGL36_pas')
soma.gegl36bar_EGL36_pas = 0.3

soma.insert('IRK_pas')
soma.girkbar_IRK_pas = 0.28083
'''
'''
soma.insert('EGL19_ca')  # add calcium channel
soma.gegl19bar_EGL19_ca = 1.59354
soma.insert('UNC2_ca')
soma.gunc2bar_UNC2_ca = 2.95138
soma.insert('CCA1_ca')
soma.gcca1bar_CCA1_ca = 3.97335  # 3.97335

soma.insert('SLO1_EGL19_cak')  # add potassium channel
soma.gslo1_egl19bar_SLO1_EGL19_cak = 0.0408
h.setpointer(soma(0.5)._ref_minf_EGL19_ca, 'minf_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_mtau_EGL19_ca, 'mtau_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_m_EGL19_ca, 'm_egl19', soma(0.5).SLO1_EGL19_cak)
h.setpointer(soma(0.5)._ref_h_EGL19_ca, 'h_egl19', soma(0.5).SLO1_EGL19_cak)
soma.insert('SLO1_UNC2_cak')
soma.gslo1_unc2bar_SLO1_UNC2_cak = 0.0408
h.setpointer(soma(0.5)._ref_minf_UNC2_ca, 'minf_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_mtau_UNC2_ca, 'mtau_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_m_UNC2_ca, 'm_unc2', soma(0.5).SLO1_UNC2_cak)
h.setpointer(soma(0.5)._ref_h_UNC2_ca, 'h_unc2', soma(0.5).SLO1_UNC2_cak)
soma.insert('SLO2_EGL19_cak')
soma.gslo2_egl19bar_SLO2_EGL19_cak = 0.31636
h.setpointer(soma(0.5)._ref_minf_EGL19_ca, 'minf_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_mtau_EGL19_ca, 'mtau_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_m_EGL19_ca, 'm_egl19', soma(0.5).SLO2_EGL19_cak)
h.setpointer(soma(0.5)._ref_h_EGL19_ca, 'h_egl19', soma(0.5).SLO2_EGL19_cak)
soma.insert('SLO2_UNC2_cak')
soma.gslo2_unc2bar_SLO2_UNC2_cak = 0.31636
h.setpointer(soma(0.5)._ref_minf_UNC2_ca, 'minf_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_mtau_UNC2_ca, 'mtau_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_m_UNC2_ca, 'm_unc2', soma(0.5).SLO2_UNC2_cak)
h.setpointer(soma(0.5)._ref_h_UNC2_ca, 'h_unc2', soma(0.5).SLO2_UNC2_cak)
'''
'''
soma.insert('Intercell_CaSK')  # add potassium channel
soma.insert('KCNL_CaSK')
soma.gkcnlbar_KCNL_CaSK = 0.5

soma.insert('leak')  # add sodium and leak channel
soma.glbar_leak = 0.3
soma(0.5).el_leak = -80  # -80 for RMD
'''
soma.insert('hh')
#################################################### Set up experiment #################################################

#
stim = h.IClamp(soma(0.5))  # add a current clamp at the middle of the soma
stim.delay = 5
stim.dur = 0.17  # ms
stim.amp = 14

soma_v = h.Vector()
soma_v.record(soma(0.5)._ref_v)

print(soma.psection())
# print(h.units('cai'))

t = h.Vector()
t.record(h._ref_t)  # record time
h.v_init = -70 # set starting voltage
h.tstop = 30  # set simulation time
h.dt = 0.0001
h.run()  # run simulation
#################################################### Plot Results ######################################################
plt.figure(figsize=(8, 5))
plt.plot(t, soma_v, color='b', label='V')
plt.xlim(0, h.tstop)
plt.xlabel('Time (ms)', fontsize=15)
plt.ylabel('V (mV)', fontsize=15)
plt.legend(fontsize=14, frameon=False)
plt.tight_layout()
plt.ioff()
plt.show()