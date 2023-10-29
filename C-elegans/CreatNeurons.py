# Created by Kangxin Hu at Beihang University, on Oct.2023.

from neuron import h, nrn, gui, rxd


class Cell:
    def __init__(self, shapePara, RPara, CPara, Para):  # grid is the number of neurons
        self._setup_morphology(shapePara)
        self.all = self.soma.wholetree()
        self._setup_biophysics(RPara, CPara, Para)
        h.define_shape()

    def _setup_morphology(self, shapePara):
        self.soma = h.Section(name='soma', cell=self)
        self.soma.L = self.soma.diam = shapePara  # 12.6157

    def _setup_biophysics(self, RPara, CPara, Para):
        for sec in self.all:
            sec.Ra = RPara
            sec.cm = CPara

        self.soma.insert('shl_pas')  # add potassium channel
        self.soma.insert('SHK1_pas')
        self.soma.insert('KVS1_pas')
        self.soma.insert('KQT3_pas')
        self.soma.insert('EGL2_pas')
        self.soma.insert('EGL36_pas')
        self.soma.insert('IRK_pas')

        self.soma.insert('Intercell_CaSK')  # add potassium channel
        self.soma.insert('KCNL_CaSK')
        self.soma.insert('leak')  # add sodium and leak channel
        self.soma.insert('ITR1_ca')
        self.soma.insert('MEC10_na')
        self.soma.insert('OSM9_trpv')
        self.soma.insert('TRP4_trpv')

        for seg in self.soma:
            seg.gshl1bar_shl_pas = Para[0]*5
            seg.gshk1bar_SHK1_pas = Para[1]*5
            seg.gkvs1bar_KVS1_pas = Para[2]*5
            seg.gkqt3bar_KQT3_pas = Para[3]*5
            seg.gegl2bar_EGL2_pas = Para[4]*5
            seg.gegl36bar_EGL36_pas = Para[5]*5
            seg.girkbar_IRK_pas = Para[6]*5

            seg.gkcnlbar_KCNL_CaSK = Para[14]
            seg.glbar_leak = Para[16]
            seg.el_leak = Para[17]  # -80 for RMD
            seg.ena = Para[18]  ##前面的这个for又是在干啥呀，看不懂
            seg.gitr1bar_ITR1_ca = Para[19]
            seg.gmec10bar_MEC10_na = Para[20]
            seg.gosm9bar_OSM9_trpv = Para[21]
            seg.gtrp4bar_TRP4_trpv = Para[22]
            seg.ek=-75
