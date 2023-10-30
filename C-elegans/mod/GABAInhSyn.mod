TITLE Syn.mod

COMMENT
1. Created by Yu Zhang at Tsinghua University, on Otc. 2021.
2. Refer to "Multistability and Long-Timescale Transients Encoded by Network Structure in a Model of C. elegans Connectome Dynamics", James M. Kunert-Graf, et.al, 11, Jun. 2017.
3. Refer to "Mathematical foundation of neuroscience", G. Bard. Ermentrout and David H. Terman, pp. 146-152.
4. Attached to the presynapse
ENDCOMMENT

NEURON {
  POINT_PROCESS GABAInhSyn
  RANGE tau, vgaba, i, d, p, taud, taup
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  tau = 0.1 (ms) <1e-9,1e9> :0.1
  vgaba = -45 (mV)
  d = 0.0053 <0,1>: depression factor
  p = 0.0096 <0, 1e9>: potentiation factor
  taud = 34 (ms) : depression effectiveness time constant
  taup = 16.8 (ms) : Bi & Poo (1998, 2001)
}

ASSIGNED {
  v (mV)
  i (nA)
  tpost (ms)
}

STATE {
  g (uS)
}

INITIAL {
  g=0
  tpost = -1e9
  net_send(0, 1)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  i = g*(v - vgaba)
  : printf(" syn i=%g\n", i)
  : printf("volatage = %g,  time = %g, g = %g,  current = %g \n", v, t, g, i)
}

DERIVATIVE state {
    g' = -g/tau
}

FUNCTION factor(Dt (ms)) { : Dt is interval between most recent presynaptic spike
    : and most recent postsynaptic spike
    : calculated as tpost - tpre (i.e. > 0 if pre happens before post)
  : the following rule is the one described by Bi & Poo
  if (Dt>0) {
    factor = 1 + p*exp(-Dt/taup) : potentiation
  } else if (Dt<0) {
    factor = 1 - d*exp(Dt/taud) : depression
  } else {
    factor = 1 : no change if pre and post are simultaneous
  }
}

: w    intrinsic synaptic weight
: k    plasticity factor (k in [0,1) for depression, k>1 for potentiation)
: tpre time of previous postsynaptic spike
NET_RECEIVE(w (uS), k, tpre (ms)) {  :from Netcon
  INITIAL { k = 1  tpre = -1e9 }
  if (flag == 0) { : presynaptic spike (after last post so depress)
        :printf("Presyn spike--entry flag=%g t=%g w=%g k=%g tpre=%g tpost=%g g=%g\n", flag, t, w, k, tpre, tpost, g)
        g = g + w*k
        tpre = t
        k = k * factor(tpost - t)
: printf("  new k %g, tpre %g\n", k, tpre)
  }else if (flag == 2) { : postsynaptic spike (after last pre so potentiate)
        :printf("Postsyn spike--entry flag=%g t=%g tpost=%g g=%g\n", flag, t, tpost, g)
        tpost = t
        FOR_NETCONS(w1, k1, tp) { : also can hide NET_RECEIVE args
: printf("entry FOR_NETCONS w1=%g k1=%g tp=%g\n", w1, k1, tp)
        k1 = k1*factor(t - tp)
: printf("  new k1 %g\n", k1)
    }
  } else { : flag == 1 from INITIAL block
        :printf("entry flag=%g t=%g g=%g\n", flag, t, g)
        WATCH (v > -20) 2
  }
}