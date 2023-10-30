TITLE shl_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gshl1bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX shl_pas
        USEION k READ ek WRITE ik    : ek is reversal voltage
        RANGE  gshl1bar, gshl1                : gshl1bar is max gshl1
        GLOBAL minf, mtau, hinf, h_ftau, h_stau, m_V, k_a, h_V, k_i, m_a, m_b, m_c, m_d, m_e, m_f, h_f_a, h_f_b, h_f_c, h_f_d, h_s_a, h_s_b, h_s_c, h_s_d
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (nS) = (nanosiemens)
}

PARAMETER {
        gshl1bar  = 2.9 (nS/cm2)   : 2.9 for AWC_ON, 2.5 for RMD
    : minf
		m_V = -6.8 (mV)        : -6.8 as an exchange
        k_a    = 14.1 (mV)
    : hinf
        h_V = -33.1 (mV)
        k_i    = 8.3 (mV)
    : mtau
        m_a    = 1.4 (ms)        : 1.4 as an exchange
        m_b    = -17.5 (mV)
        m_c    = 12.9 (mV)
        m_d    = -3.7 (mV)
        m_e    = 6.5 (mV)
        m_f    = 0.2 (ms)         : 0.2 as an exchange
    : htao_fast
        h_f_a  = 53.9 (ms)       : 53.9 as an exchange
        h_f_b  = -28.2 (mV)
        h_f_c  = 4.9 (mV)
        h_f_d  = 2.7 (ms)
    : htao_slow
        h_s_a  = 842.2 (ms)      : 842.2 as an exchange
        h_s_b  = -37.7 (mV)
        h_s_c  = 6.4 (mV)
        h_s_d  = 11.9 (ms)       : 11.9 as an exchange
}

STATE {
        m h_f h_s
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gshl1 (nS/cm2)
		ik (pA/cm2)
        minf hinf
	    mtau (ms)
		h_ftau (ms)
		h_stau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gshl1 = gshl1bar*m*m*m*(0.7*h_f+0.3*h_s)
	    ik = gshl1*(v -ek)
}

INITIAL {
	rates(v)
	m = 0
	h_f = 1
	h_s = 1
}

? states
DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
        h_f' = (hinf-h_f)/h_ftau
		h_s' = (hinf-h_s)/h_stau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10
        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)
        : TABLE minf, mtau, hinf, h_ftau, h_stau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" SHL-1 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = m_a/(exp((m_b-v)/m_c)+exp((v-m_d)/m_e))+m_f   : multiplied by q10, if consider temperature

                :"h_fast" SHL-1 type potassium inactivation system
		hinf = 1/(1+exp((v-h_V)/k_i))
		h_ftau =h_f_a/(1+exp((v-h_f_b)/h_f_c))+h_f_d         : multiplied by q10, if consider temperature

				:"h_fast" SHL-1 type potassium inactivation system
		h_stau =h_s_a/(1+exp((v-h_s_b)/h_s_c))+h_s_d         : multiplied by q10, if consider temperature
}

UNITSON
