TITLE KQT3_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gkqt3bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX KQT3_pas
        USEION k READ ek WRITE ik    : ek is reversal voltage
        RANGE  gkqt3bar, gkqt3       : gkqt3bar is max gkqt3
        GLOBAL minf, m_V, mk_a, m_ftau, m_ftau_a, m_ftau_b, m_ftau_c, m_stau, m_stau_a, m_stau_b, m_stau_c, m_stau_d, m_stau_e, m_stau_f, m_stau_g, winf, w_V, wk_i, winf_a, winf_b, sinf, s_V, sk_i, sinf_a, sinf_b, wtau, wtau_a, wtau_b, wtau_c, wtau_d, stau, stau_a
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (nS) = (nanosiemens)
}

PARAMETER {
        gkqt3bar  =0.55 (nS/cm2)   : 0.55 for AWC_ON
    : minf
		m_V = 7.7 (mV)        : 7.7 as an exchange
        mk_a = 15.8 (mV)
    : m_ftau
        m_ftau_a   = 39.5 (ms)        : 39.5 as an exchange
        m_ftau_b   = 38.1 (mV)
        m_ftau_c   = 33.6 (mV)
	: m_stau
        m_stau_a   = 550.3 (ms)        : 550.3 as an exchange
        m_stau_b   = -534.5 (ms)        : 554.5 as an exchange
        m_stau_c   = 35.3357 (mV)      : 1/x
        m_stau_d   = -23.9 (mV)         : 39.5 as an exchange
        m_stau_e   = -459.1 (ms)          : 459.1 as an exchange
        m_stau_f   = 28.0112 (mV)      : 1/x
        m_stau_g   = 14.2 (mV)
    : winf
        w_V = -1.1 (mV)
        wk_i   = 28.8 (mV)
        winf_a = 0.5
        winf_b = 0.5
    : sinf
        s_V = -45.3 (mV)
        sk_i    = 12.3 (mV)
        sinf_a = 0.3
        sinf_b = 0.7
    : wtau
		wtau_a = 0.5 (ms)
		wtau_b = 2.9 (ms)
		wtau_c = -48.1 (mV)
		wtau_d = 48.8 (mV)
    : stau
		stau_a = 500 (ms)
}

STATE {
        m_f m_s w s
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gkqt3 (nS/cm2)
		ik (mA/cm2)
        minf  winf  sinf
	    m_ftau (ms) m_stau (ms) wtau (ms) stau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkqt3 = gkqt3bar*(0.3*m_f+0.7*m_s)*w*s
	    ik = gkqt3*(v -ek)
}

INITIAL {
	rates(v)
	m_f = 0
	m_s = 0
	  w = 0
      s = 0
}

? states
DERIVATIVE states {
        rates(v)
        m_f' = (minf-m_f)/m_ftau
        m_s' = (minf-m_s)/m_stau
		w' = (winf-w)/wtau
		s' = (sinf-s)/stau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)
                :"m_fast" KQT-3 type potassium activation system
		minf = 1/(1+exp((m_V-v)/mk_a))
        m_ftau = m_ftau_a/(1+((v+m_ftau_b)/m_ftau_c)^2)

                :"m_slow" KQT-3 type potassium activation system
        m_stau = m_stau_a+m_stau_b/(1+10^(-1*(m_stau_d-v)/m_stau_c))+m_stau_e/(1+10^(-1*(v+m_stau_g)/m_stau_f)) : multiplied by q10, if consider temperature

                :"w" KQT-3 type potassium inactivation system
		winf = winf_a+winf_b/(1+exp((v-w_V)/wk_i))
		wtau = wtau_a+wtau_b/(1+((v-wtau_c)/wtau_d)*((v-wtau_c)/wtau_d))         : multiplied by q10, if consider temperature

				:"s" KQT-3 type potassium inactivation system
		sinf = sinf_a+sinf_b/(1+exp((v-s_V)/sk_i))
		stau = stau_a         : multiplied by q10, if consider temperature
}

UNITSON