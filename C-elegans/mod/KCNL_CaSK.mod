TITLE KCNL_CaSK.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gkcnlbar, ekcnl, celsius
ENDCOMMENT

NEURON	{
	SUFFIX KCNL_CaSK
    USEION ca READ cai
    USEION k READ ek WRITE ik
	RANGE gkcnlbar, gkcnl
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (nS) = (nanosiemens)
}

PARAMETER {
        gkcnlbar  = 0.06 (nS/cm2)
}

STATE {
        m
}

ASSIGNED {
        cai (mM)
        v (mV)
		ek (mV)
        celsius (degC)
        ekcnl (mV)
        gkcnl (pS/cm2)
		ik (pA/cm2)
        minf
	    mtau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkcnl = gkcnlbar*m
	    ik = gkcnl*(v -ek)
}

INITIAL {
	:rates(v)
	m = 0.13563
}

? states
DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
}

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10
        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" EGL-2 type potassium activation system
		minf = cai/(0.00033+cai)
        mtau = 6.3   : multiplied by q10, if consider temperature
}

UNITSON
