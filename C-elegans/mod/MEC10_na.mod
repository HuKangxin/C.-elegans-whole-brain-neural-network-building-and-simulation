TITLE MEC10_na.mod Sodium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gmec10bar, eca, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX MEC10_na
        USEION na READ ena WRITE ina    : eca is reversal voltage
        RANGE  gmec10bar, gmec10                : gcca1bar is max gcca1
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
        gmec10bar  = 1.7 (nS/cm2)   : 0.7 for AWC_ON, 3.1 for RMD
}

STATE {
        m h
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        gmec10 (nS/cm2)
	ina (pA/cm2)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gmec10 = gmec10bar*m*m*m*h
	ina = gmec10*(v - 50)
}

INITIAL {
	m = alpha(v)/(alpha(v)+beta(v))
	h = 0.5  :hinf
}

? states
DERIVATIVE states {
        m' = (1-m)*alpha(v)-m*beta(v)
        h' = (1-h)*alpha(v)-h*beta(v)
}


FUNCTION alpha(v(mV))(/ms) {  
        LOCAL x
        UNITSOFF
        x=(v+55)/10
        if (fabs(x)>1e-6) {
                alpha=0.01*(v+55)/(1-exp(-(v+55)/10))
        } else {
                alpha=0.1/(1-0.5*x)
        }
        UNITSON
}

FUNCTION beta(v(mV))(/ms) {
        UNITSOFF
        beta=0.125*exp(-(v+55)/80)
        UNITSON
}