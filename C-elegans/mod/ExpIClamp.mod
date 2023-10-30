COMMENT
An event-driven current source that delivers an electrode current which decays monoexponentially with time.
Arrival of an event with weight w instantaneously perturbs the current by w nA.
For example, a single event with weight 0.1 will elicit a single current pulse with peak amplitude 0.1 nA.

Developed by NTC based on ExpSyn.
ENDCOMMENT

NEURON {
        POINT_PROCESS ExpIClamp
        RANGE tau_1,tau_2, i
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
        tau_1 = 0.2 (ms) <1e-9,1e9>
        tau_2 = 0.1 (ms) <1e-9,1e9>
}

ASSIGNED {
        i (nA)
}

STATE {
        i_1 (nA)
        i_2 (nA)
}

INITIAL {
        i_1=0
        i_2=0
}

BREAKPOINT {
        SOLVE state METHOD cnexp
        if (t>=0(ms) && t<100(ms)) {
                i = i_1 + 150
        }   else   {
                i=i_2+20
        }
}

DERIVATIVE state {
        i_1' = -i_1/tau_1
        i_2' = -i_2/tau_2
}
