TITLE EGL2_pas.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gegl2bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX EGL2_pas
        USEION k READ ek WRITE ik     : ek is reversal voltage
        RANGE  gegl2bar, gegl2        : gegl2bar is max gegl2
        GLOBAL minf, mtau, m_V, k_a, m_a, m_b, m_c, m_d
	THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        gegl2bar = 0.85 (nS/cm2)
    : minf
	m_V = -6.9 (mV)
        k_a = 14.9 (mV)
    : mtau
        m_a = 8.39 (ms)        : 8.4 as an exchange
        m_b = -122.6 (mV)
        m_c = 13.8 (mV)
        m_d = 4.04 (ms)         : 4.1 as an exchange
}

STATE {
        m
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gegl2 (nS/cm2)
	ik (pA/cm2)
        minf
	mtau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gegl2 = gegl2bar*m
	ik = gegl2*(v-ek)
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

		:"m" EGL-2 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = m_a/(1+exp((v-m_b)/m_c))+m_d   : multiplied by q10, if consider temperature
}

UNITSON