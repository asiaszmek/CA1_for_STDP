: Tsodyks-Markram model

TITLE TM depressing synapse

:use and t_rec taken from John F. Wesseling and Donald C. Lo
:other parameters are fitted for consistency with the traditional
:NMDA model
NEURON {
	POINT_PROCESS depletion
        RANGE R, T
	RANGE t_rec, use, ase, factor

	RANGE tau1, tau2
	NONSPECIFIC_CURRENT i
}

UNITS {
	(mM) = (milli/liter)			
        (nA) = (nanoamp)
}

PARAMETER {
	t_rec = 4200 (ms)	: recovery of synaptic vesicles
	use = 0.044 : probability of release
	ase = 36	(mM) : availability of resources
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 1 (ms) <1e-9,1e9>
        gmax = 1 (nA/mM)

}

INITIAL {
     LOCAL tp
     if (tau1/tau2 > .9999) {
	tau1 = .9999*tau2
     }
     A = 0
     B = 0
     tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
     factor = -exp(-tp/tau1) + exp(-tp/tau2)
     factor = 1/factor
     R = 1
     T = 0
}

ASSIGNED {
	i (nA)
	factor
}

STATE {
        A
        B
        R (1)
	T (ms)
}

BREAKPOINT {
     SOLVE state METHOD cnexp
     T = ase*(B-A)
     i = T*gmax
}

DERIVATIVE state {
	R' = (1-R)/t_rec
	A' = -A/tau1
	B' = -B/tau2

}

NET_RECEIVE (weight (uS)){
	state_discontinuity(R, (1-use)*R)
	state_discontinuity(A, A + (use*R)*factor)
	state_discontinuity(B, B + (use*R)*factor)


}
