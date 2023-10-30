NEURON {
	POINT_PROCESS WormExcSyn
	RANGE tau, e, i, gmax, g, s, k
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
	e = 0(mV)
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
	:printf("t= %g, k = %g, s=%g, Vpre=%g \n", t, k, s, Vpre)
    :printf("t = %g, k = %g, synv = %g\n", t, k, v)
}

DERIVATIVE state {
    :s' = 3*(1-s)/(1+exp(-0.1*(Vpre-95)))-0.1*s
	if (k==0){
		s' = 0
	} else {
		s' = 3*(1-s)/(1+exp(-0.1*(Vpre-95)))-0.1*s
	}
}

NET_RECEIVE(weight (uS)) {
  if (flag == 0) { : presynaptic spike (after last post so depress)
        :s = s + weight
        k = 1
		:printf("t= %g, k = %g, s=%g, Vpre=%g \n", t, k, s, Vpre)
  }else if (flag == 2) { : postsynaptic spike (after last pre so potentiate)
        k = 0
  } else { : flag == 1 from INITIAL block
        WATCH (g < 0.001) 2
  }
}