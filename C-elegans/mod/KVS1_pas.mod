TITLE KVS1_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gkvs1bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX KVS1_pas
        USEION k READ ek WRITE ik    : ek is reversal voltage
        RANGE  gkvs1bar, gkvs1                : gkvs1bar is max gkvs1
        GLOBAL minf, mtau, hinf, htau, m_V, k_a, h_V, k_i, m_a, m_b, m_c, m_d, h_a, h_b, h_c, h_d
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (nS) = (nanosiemens)
}

PARAMETER {
        gkvs1bar  = 0.8 (nS/cm2)   : 0.8 for AWC_ON
    : minf
		m_V = 27.1 (mV)        : 27.1 as an exchange
        k_a    = 25.0 (mV)
    : hinf
        h_V = 17.3 (mV)        : 17.3 as an exchange
        k_i    = 11.1 (mV)
    : mtau
        m_a    = 3.0 (ms)        : 3.0 as an exchange
        m_b    = 18.1 (mV)
        m_c    = 20 (mV)
        m_d    = 0.1 (mV)         : 0.1 as an exchange
    : htau
        h_a    = 8.9 (ms)        : 8.9 as an exchange
        h_b    = 50.0 (mV)
        h_c    = 15.0 (mV)
        h_d    = 5.3 (ms)        : 5.3 as an exchange
}

STATE {
        m h
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gkvs1 (nS/cm2)
		ik (pA/cm2)
        minf hinf
	    mtau (ms) htau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkvs1 = gkvs1bar*m*h
	    ik = gkvs1*(v -ek)
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
        : TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" KVS-1 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = m_a/(1+exp((v-m_b)/m_c))+m_d   : multiplied by q10, if consider temperature

                :"h" KVS-1 type potassium inactivation system
		hinf = 1/(1+exp((v-h_V)/k_i))
		htau = h_a/(1+exp((v-h_b)/h_c))+h_d   : multiplied by q10, if consider temperature
}

UNITSON