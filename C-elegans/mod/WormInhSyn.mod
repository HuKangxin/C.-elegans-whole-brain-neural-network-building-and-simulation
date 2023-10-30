NEURON {
	POINT_PROCESS WormInhSyn
	RANGE tau, e, i, gmax, g,s, k
	POINTER Vpre
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = -80	(mV)
	gmax = 500 (pS)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (pS)
	k
	Vpre (mV)
}

STATE {
	s
}

INITIAL {
	s=0
    net_send(0, 1)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gmax*s*k
	i = g*(v-e)
    :printf("t = %g,  g = %g,  syni = %g,  synv = %g\n", t, g, i, v)
}

DERIVATIVE state {
    :s' = 5*(1-s)/(1+exp(-0.1*(Vpre-95)))-0.15*s
	if (k==0){
		s' = 0
	} else {
		s' = 5*(1-s)/(1+exp(-0.1*(Vpre-95)))-0.15*s
	}
}

NET_RECEIVE(weight (uS)) {
  if (flag == 0) { : presynaptic spike (after last post so depress)
        :s = s + w
        k = 1
  }else if (flag == 2) { : postsynaptic spike (after last pre so potentiate)
        k = 0
  } else { : flag == 1 from INITIAL block
        WATCH (g < 0.01) 2
  }
}