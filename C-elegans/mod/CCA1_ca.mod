TITLE CCA1_ca.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gcca1bar, eca, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX CCA1_ca
        USEION ca READ eca WRITE ica    : eca is reversal voltage
        RANGE  gcca1bar, gcca1                : gcca1bar is max gcca1
        GLOBAL minf, mtau, hinf, htau, m_V, k_a, h_V, k_i, m_a, m_b, m_c, m_d, h_a, h_b, h_c, h_d
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        gcca1bar  = 0.7 (nS/cm2)   : 0.7 for AWC_ON, 3.1 for RMD
    : minf
	m_V = -57.7 (mV)        : -57.7 as an exchange
        k_a =2.4 (mV)           : 2.4 as an exchange
    : hinf
        h_V = -73.0 (mV)         : -73.0 as an exchange
        k_i = 8.1 (mV)           : 8.1 as an exchange
    : mtau
        m_a = 20.0 (ms)        : 20.0 as an exchange
        m_b = -92.5 (mV)       : -92.5 as an exchange
        m_c = 21.1 (mV)       : 21.1 as an exchange
        m_d = 0.4 (ms)         : 0.4 as an exchange
    : htau
        h_a = 22.4 (ms)        : 22.4 as an exchange
        h_b = -75.7 (mV)      : -75.7 as an exchange
        h_c = 9.4 (mV)        : 9.4 as an exchange
        h_d = 1.6 (ms)       : 1.6 as an exchange
}

STATE {
        m h
}

ASSIGNED {
        v (mV)
        celsius (degC)
        eca (mV)
        gcca1 (nS/cm2)
	ica (pA/cm2)
        minf 
        hinf
	mtau (ms)
        htau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gcca1 = gcca1bar*m*m*h
	ica = gcca1*(v - 45)
}

INITIAL {
	rates(v)
	m = 0  :minf
	h = 1  :hinf
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

                :"m" CCA-1 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = m_a/(1+exp((m_b-v)/m_c))+m_d   : multiplied by q10, if consider temperature

                :"h" CCA-1 type potassium inactivation system
		hinf = 1/(1+exp((v-h_V)/k_i))
		htau = h_a/(1+exp((v-h_b)/h_c))+h_d   : multiplied by q10, if consider temperature
UNITSON
}

