TITLE EGL36_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gegl36bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX EGL36_pas
        USEION k READ ek WRITE ik    : ek is reversal voltage
        RANGE  gegl36bar, gegl36     : gegl36bar is max gegl36
        GLOBAL minf, m_V, k_a, m_ftau, m_mtau, m_stau, m_ftau_a, m_mtau_a, m_stau_a
	THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        gegl36bar  = 1.3 (nS/cm2)   : 1.3 for RMD
    : minf
	m_V = 63.0 (mV)
        k_a = 28.5 (mV)
    : mtau
        m_stau_a    =  355.0 (ms)
        m_mtau_a    =  63.0 (ms)
        m_ftau_a    =  13.0 (ms)
}

STATE {
        m_f m_m  m_s
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gegl36 (nS/cm2)
	ik (pA/cm2)
        minf
	m_ftau (ms)
	m_mtau (ms)
	m_stau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gegl36 = gegl36bar*(0.33*m_f+0.36*m_m+0.39*m_s)
	ik = gegl36*(v -ek)
}

INITIAL {
	rates(v)
	m_f = 0
	m_m = 0
	m_s = 0
}

? states
DERIVATIVE states {
        rates(v)
        m_f' = (minf-m_f)/m_ftau
        m_m' = (minf-m_m)/m_mtau
        m_s' = (minf-m_s)/m_stau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10
        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)
        : TABLE minf, m_ftau, m_mtau, m_stau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" EGL-36 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        m_ftau = m_ftau_a   : multiplied by q10, if consider temperature
		m_mtau = m_mtau_a   : multiplied by q10, if consider temperature
		m_stau = m_stau_a   : multiplied by q10, if consider temperature
}

UNITSON
