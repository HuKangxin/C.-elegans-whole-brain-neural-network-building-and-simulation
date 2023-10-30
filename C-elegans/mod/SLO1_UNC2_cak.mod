TITLE SLO1_UNC2_cak.mod Potassium ion channel of HH model in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.
3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.
4. gslo1_unc2bar, ek, celsius
ENDCOMMENT

? interface
NEURON {
        SUFFIX SLO1_UNC2_cak
        USEION k READ ek WRITE ik
		POINTER minf_unc2, mtau_unc2, m_unc2, h_unc2
        RANGE  gslo1_unc2bar, gslo1_unc2                     : gslo1_unc2bar is max gslo1_unc2
        GLOBAL minf, mtau
		THREADSAFE : assigned GLOBALs will be per thread
}

UNITS {
        (pA) = (picoamp)
        (mV) = (millivolt)
	    (pS) = (picosiemens)
}

PARAMETER {
        gslo1_unc2bar  = 0.11 (pS/cm2)    : 0.11 AWC; 0.3 RMD
}

STATE {
        m
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gslo1_unc2 (pS/cm2)
		ik (pA/cm2)
		minf
	    mtau (ms)
        minf_unc2
		mtau_unc2  (ms)
		m_unc2
		h_unc2
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gslo1_unc2 = gslo1_unc2bar*m*h_unc2
	    ik = gslo1_unc2*(v -ek)
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
PROCEDURE rates(v(mV)) {  :Computes rate and other constants   at current v.
                          :Call once from HOC to initialize inf at resting v.
        LOCAL  q10,  sgn, Ca_oBK, k1, k2, k3, alpha, beta

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

                :"m" SLO1_UNC2 type potassium activation system
		if ((v-60)> 0) { sgn = 1 } else { sgn = -1}
		Ca_oBK = (10^6)*(10^-3)*40*(10^-15)*sgn*(v-60)/(8*3.1415926535897932*13*(10^-9)*250*(10^-12)*96485)*exp(-0.100697567)+0.05
		k1= 0.156217*exp(0.028*v)/(1+(55.73/Ca_oBK)^1.30)   :k_o_+
		k2= 3.15*exp(-0.013*v)/(1+(0.05/34.34)^0.0001)   :k_c_-
        k3= 3.15*exp(-0.013*v)/(1+(Ca_oBK/34.34)^0.0001)    :k_o_-
        alpha = minf_unc2/mtau_unc2
        beta = 1/mtau_unc2-alpha
		minf = m_unc2*k1*(alpha+beta+k2)/((k1+k3)*(k2+alpha)+beta*k2)
		mtau = (alpha+beta+k2)/((k1+k3)*(k2+alpha)+beta*k2)  : multiplied by q10, if consider temperatur
}

UNITSON