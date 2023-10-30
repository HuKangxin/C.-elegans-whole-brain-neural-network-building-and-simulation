TITLE na_leak.mod na ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gegl2bar, eegl2, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX leak
	NONSPECIFIC_CURRENT il
        RANGE  glbar, el, il
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        glbar  = 2700 (nS/cm2)    :0.04
}

STATE {
        m
}

ASSIGNED {
        v (mV)
        celsius (degC)
        el (mV)
	il (pA/cm2)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
	il= glbar*(v - el)
}

INITIAL {
	m = 1
}

DERIVATIVE states {
        m' = 0
}