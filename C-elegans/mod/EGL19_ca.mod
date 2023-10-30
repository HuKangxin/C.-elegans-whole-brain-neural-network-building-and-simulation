TITLE EGL19_ca.mod Calcium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gegl19bar, eca, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX EGL19_ca
        USEION ca READ eca WRITE ica    : eca is reversal voltage
        RANGE  gegl19bar, gegl19, minf, mtau, hinf, htau
	GLOBAL m_V, k_a,  h_V, k_i, k_i_b, h_V_b, hinf_a, hinf_b,hinf_c, hinf_d, mtau_a, mtau_b, mtau_c, mtau_d, mtau_e, mtau_f, mtau_g, htau_a, htau_b, htau_c, htau_d, htau_e, htau_f, htau_g, htau_h
	THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        gegl19bar  = 1.55 (nS/cm2)   : 1.55 for AWC_ON, 0.99 for RMD
    : minf
	m_V = -4.4 (mV)        : -4.4 as an exchange
        k_a = 7.5 (mV)
    : hinf
        h_V = 14.9 (mV)       : 14.9 as an exchange
        k_i = 12.0 (mV)
        h_V_b = -20.5 (mV)      : -20.5 as an exchange
        k_i_b = 8.1 (mV)
        hinf_a = 1.43
        hinf_b = 0.14
        hinf_c = 5.96
        hinf_d = 0.6
    : mtau
        mtau_a    = 2.9 (ms)
        mtau_b    = -4.8 (mV)    : -4.8 as an exchange
        mtau_c    = 6.0 (mV)
        mtau_d    = 1.9 (ms)
        mtau_e    = -8.6 (mV)    : -8.6 as an exchange
        mtau_f    = 30.0 (mV)
        mtau_g    = 2.3 (ms)
	: htau
        htau_a    = 0.4
        htau_b    = 44.6 (ms)
        htau_c    = -33.0 (mV)   : -33.0 as an exchange
        htau_d    = 5.0 (mV)
        htau_e    = 36.4 (ms)
        htau_f    = 18.7 (mV)   : 18.7 as an exchange
        htau_g    = 3.7 (mV)
        htau_h    = 43.1 (ms)
}

STATE {
        m h
}

ASSIGNED {
        v (mV)
        celsius (degC)
        eca (mV)
        gegl19 (nS/cm2)
	ica (pA/cm2)
        minf
        hinf
	mtau (ms)
        htau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gegl19 = gegl19bar*m*h
	ica = gegl19*(v - 45)
}

INITIAL {
	rates(v)
	m =  0  : minf
	h =  1  : hinf
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

                :"m" EGL-19 type potassium activation system
		minf = 1/(1+exp((m_V-v)/k_a))
        mtau = mtau_a*exp(-1*((v-mtau_b)/mtau_c)*((v-mtau_b)/mtau_c))+mtau_d*exp(-1*((v-mtau_e)/mtau_f)*((v-mtau_e)/mtau_f))+mtau_g  : multiplied by q10, if consider temperature

                :"h" EGL-19 type potassium inactivation system
		hinf = (hinf_a/(1+exp((h_V-v)/k_i))+hinf_b)*(hinf_c/(1+exp((v-h_V_b)/k_i_b))+hinf_d)
		htau = htau_a*(htau_b/(1+exp((v-htau_c)/htau_d))+htau_e/(1+exp((v-htau_f)/htau_g))+htau_h)      : multiplied by q10, if consider temperature
}

UNITSON