TITLE SHK_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gshk1bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX SHK1_pas
        USEION k READ ek WRITE ik    : ek is reversal voltage
        RANGE  gshk1bar, gshk1                : gshk1bar is max gshk1
        GLOBAL minf, mtau, hinf, htau, m_V, k_a, h_V, k_i, m_a, m_b, m_c, m_d, m_e, m_f, h_a
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (nS) = (nanosiemens)
}

PARAMETER {
        gshk1bar  = 0.1 (nS/cm2)   : 0.1 for AWC_ON, 1.1 for RMD
    : minf
		m_V = 20.4 (mV)
        k_a    = 7.7 (mV)
    : hinf
        h_V = -7.0 (mV)
        k_i    = 5.8 (mV)
    : mtau
        m_a    = 26.6 (ms)
        m_b    = -33.7 (mV)
        m_c    = 15.8 (mV)
        m_d    = -33.7 (mV)
        m_e    = 15.4 (mV)
        m_f    = 2.0 (ms)
    : htau
        h_a    = 1400 (ms)
}

STATE {
        m h
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gshk1 (nS/cm2)
		ik (pA/cm2)
        minf hinf
	    mtau (ms) htau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gshk1 = gshk1bar*m*h
	    ik = gshk1*(v-ek)
}

INITIAL {
	rates(v)
	m = 0
	h = 1
}

? states
DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10
        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)
        :TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" SHK-1 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = m_a/(exp((m_b-v)/m_c)+exp((v-m_d)/m_e))+m_f   : multiplied by q10, if consider temperature

                :"h" SHK-1 type potassium inactivation system
		hinf = 1/(1+exp((v-h_V)/k_i))
		htau = h_a                                           : multiplied by q10, if consider temperature
}

UNITSON