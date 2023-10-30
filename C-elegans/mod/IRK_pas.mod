TITLE IRK_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. girkbar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX IRK_pas
        USEION k READ ek WRITE ik    : ek is reversal voltage
        RANGE  girkbar, girk                : girkbar is max girk
        GLOBAL minf, m_V, k_a, mtau, mtau_a, mtau_b, mtau_c, mtau_d, mtau_e, mtau_f
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        girkbar  = 0.25 (nS/cm2)   : 0.25 for AWC_on, 0.2 for RMD
    : minf
	m_V    = -82 (mV)
        k_a    = 13 (mV)
    : mtau
        mtau_a    = 17.1 (ms)
        mtau_b    = -17.8 (mV)
        mtau_c    = 20.3 (mV)
        mtau_d    = -43.4 (mV)
        mtau_e    = 11.2 (mV)
        mtau_f    = 3.8 (ms)
}

STATE {
        m
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        girk (nS/cm2)
	ik (pA/cm2)
        minf
	    mtau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        girk = girkbar*m
	ik = girk*(v -ek)
}

INITIAL {
	rates(v)
	m = 0
}

? states
DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10
        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)
        :TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" IRK type potassium activation system
		minf = 1/(1+exp((v-m_V)/k_a))
        mtau = mtau_a/(exp((mtau_b-v)/mtau_c)+exp((v-mtau_d)/mtau_e))+mtau_f   : multiplied by q10, if consider temperature
}

UNITSON
