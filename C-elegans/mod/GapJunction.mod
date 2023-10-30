NEURON {
	POINT_PROCESS GapJunction
	RANGE i, gmax, s
	POINTER Vpre
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
    gmax = 100 (pS)
}

ASSIGNED {
	v (mV)
	i (nA)
    k
	Vpre (mV)
}

STATE {
	s
}

INITIAL {
	s=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = gmax*(v-Vpre)
}

DERIVATIVE state {
    s' = 0
}