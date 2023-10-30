TITLE Intercell_CaSK.mod Intracellular calcium ion concentration of neurons in C.elegans

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.
2. Refer to "Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD", Plos One, 14(7), Jul. 2019.

ENDCOMMENT

NEURON	{
	SUFFIX Intercell_CaSK
	USEION ca  READ ica WRITE cai
    GLOBAL cellvolum, H
}

UNITS	{
	(mV) = (millivolt)
	(pA) = (picoamp)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um)	= (micron)
}

PARAMETER	{
	cellvolum = 31.16    :unite
}

ASSIGNED	{
		ica (pA/cm2)
		H
}

STATE	{
	cai (mM)
}

BREAKPOINT	{ SOLVE states METHOD cnexp
}

INITIAL {
	cai=0.00005
}

DERIVATIVE states	{
UNITSOFF
	if (ica >= 0) { H = 1 } else {H=0}
	cai' = (1-H)*(-0.001*ica/(2*cellvolum*96.485)-(cai-0.00005)/50)+H*((0.00005-cai)/50)  : frady constant is 96.485
}
UNITSON