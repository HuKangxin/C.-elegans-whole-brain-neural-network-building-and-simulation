TITLE SLO2_EGL19_cak.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gslo2_egl19bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX SLO2_EGL19_cak
        USEION k READ ek WRITE ik
		POINTER minf_egl19, mtau_egl19, m_egl19, h_egl19
        RANGE  gslo2_egl19bar, gslo2_egl19                     : gslo2_egl19bar is max gslo2_egl19
        GLOBAL minf, mtau
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (pS) = (picosiemens)
}

PARAMETER {
        gslo2_egl19bar  = 0.1 (pS/cm2)    : 0.1 AWC; 0.3 RMD
}

STATE {
        m
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gslo2_egl19 (pS/cm2)
		ik (pA/cm2)
		minf
	    mtau (ms)
        minf_egl19
		mtau_egl19  (ms)
		m_egl19
		h_egl19
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gslo2_egl19 = gslo2_egl19bar*m*h_egl19
	    ik = gslo2_egl19*(v -ek)
}

INITIAL {
	:rates(v)
	m = 0
}

? states
DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants   at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, sgn, Ca_oBK, k1, k2, k3, alpha, beta
        : Ca_oBK is nanoscale

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" SLO1_EGL19 type potassium activation system
		if ((v-60)> 0) { sgn = 1 } else { sgn = -1}
		Ca_oBK = (10^6)*(10^-3)*40*(10^-15)*sgn*(v-60)/(8*3.1415926535897932*13*(10^-9)*250*(10^-12)*96485)*exp(-0.100697567)+0.05
		k1= 0.027*exp(0.024*v)/(1+(93.45/Ca_oBK)^1.84)   :k_o_+
		k2= 0.9*exp(-0.019*v)/(1+(0.05/3294.55)^0.00001)   :k_c_-
        k3= 0.87*exp(-0.019*v)/(1+(Ca_oBK/3294.55)^0.00001)    :k_o_-
        alpha = minf_egl19/mtau_egl19
        beta = 1/mtau_egl19-alpha
		minf = m_egl19*k1*(alpha+beta+k2)/((k1+k3)*(k2+alpha)+beta*k2)
		mtau = (alpha+beta+k2)/((k1+k3)*(k2+alpha)+beta*k2)  : multiplied by q10, if consider temperatur
}

UNITSON
