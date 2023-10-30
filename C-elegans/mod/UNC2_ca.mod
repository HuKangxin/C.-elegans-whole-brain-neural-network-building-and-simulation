TITLE UNC2_ca.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gunc2bar, eca, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX UNC2_ca
        USEION ca READ eca WRITE ica    	  : eca is reversal voltage
        RANGE  gunc2bar, gunc2, minf, mtau, hinf, htau
		GLOBAL m_V, k_a, h_V, k_i, m_a, m_b, m_c, m_d, m_e, h_a, h_b, h_c, h_d, h_e, h_f
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (nS) = (nanosiemens)
}

PARAMETER {
        gunc2bar  = 1.0 (nS/cm2)   : 1.0 for AWC_ON, 0.9 for RMD
    : minf
		m_V    = -37.2 (mV)        : -37.2 as an exchange
        k_a    = 4.0 (mV)
    : hinf
        h_V    = -77.5 (mV)        : -77.5 as an exchange
        k_i    = 5.6 (mV)
    : mtau
        m_a    = 1.5 (ms)
        m_b    = -38.2 (mV)       : 3.0 as an exchange
        m_c    = 9.1 (mV)
        m_d    = 15.4 (mV)         : 0.1 as an exchange
        m_e    = 0.1 (ms)         : 0.1 as an exchange
    : htau
        h_a    = 142.5 (ms)        : 142.5 as an exchange
        h_b    = 22.9 (mV)        : 22.9 as an exchange
        h_c    = -3.5 (mV)
        h_d    = 122.6 (ms)        : 122.6 as an exchange
        h_e    = -6.1 (mV)        : -6.1 as an exchange
        h_f    = -3.6  (mV)
}

STATE {
        m h
}

ASSIGNED {
        v (mV)
        celsius (degC)
        eca (mV)
        gunc2 (nS/cm2)
		ica (pA/cm2)
        minf hinf
	    mtau (ms) htau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gunc2 = gunc2bar*m*h
	    ica = gunc2*(v - 45)
}

INITIAL {
	rates(v)
	m = 0  : minf
	h = 1  : hinf
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
        : TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" UNC-2 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = m_a/(exp((m_b-v)/m_c)+exp((v-m_b)/m_d))+m_e   : multiplied by q10, if consider temperature

                :"h" UNC-2 type potassium inactivation system
		hinf = 1/(1+exp((v-h_V)/k_i))
		htau = h_a/(1+exp((h_b-v)/h_c))+h_d/(1+exp((v-h_e)/h_f))   : multiplied by q10, if consider temperature
}

UNITSON