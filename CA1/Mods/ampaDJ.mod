TITLE detailed model of glutamate AMPA receptors

COMMENT
-----------------------------------------------------------------------------

	Kinetic model of AMPA receptors
	===============================

	7-state gating model:
	taken from
 	Diamond and Jahn (1997)
  
	C0 ---- CA1 -- CA2 -- O
	         |     |      |
      	        DA1 -- DA2   DA2b

Implemented by: Asia Jedrzejewska-Szmek
-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS AMPADJ
	POINTER T
	RANGE C0, CA1, CA2, DA1, DA2, DA2b, O, glu
	RANGE g, gmax, rb
	RANGE Erev
	RANGE Ka		    
	GLOBAL ka, k_a, k1, k_1, k2, k_2, k3, k_3
	GLOBAL alpha, beta
	GLOBAL vmin, vmax
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
        (pS) = (picosiemens)
	(uS) = (microsiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= 0    (mV)	: reversal potential
	gmax	= 5e-4  (uS)	: maximal conductance
	vmin = -120	(mV)
	vmax = 100	(mV)
	
: Rates

	ka	= 13.3   (/mM-ms): binding 
	k_a     = 6.4	(/ms)   : diffusion limited (DO NOT ADJUST)
	k1	= 0.361  (/ms)	: desensitization
	k_1	= 0.02  (/ms)	: resensitization
	k2	= 0.65   (/ms)	: desensitization
	k_2	= 0.018 (/ms)	: resensitization 
	k3	= 1    (/ms)	: desensitization
	k_3	= 3    (/ms)	: resensitization
	alpha	= 1.1    (/ms)	: opening
	beta	= 5.7    (/ms)	: closing
}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
											T (mM) : transmitter       
	glu 		(mM)		: pointer to glutamate concentration

	Ka		(/ms)    : binding
}

STATE {
	: Channel states (all fractions)
	C0		: unbound
	CA1		: single glu bound
	CA2		: double glu bound
 	DA1		: single glu bound, desensitized
 	DA2		: double glu bound, desensitized
	DA2b            : desensitized
	O		: open state 2
}

INITIAL {
	C0=1
	CA1=0
	CA2=0
	DA1=0
        DA2=0
	DA2b=0
	O=0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = (1e+06)*(gmax * O)
	i =  (1e-06)*(g * (v - Erev))
}

KINETIC kstates {
     if (T<0) {glu = 0}
	else {glu = T}
	Ka = ka * glu 

	~ C0  <-> CA1	(2*Ka, k_a)
	~ CA1 <-> CA2	(Ka, 2*k_a)
	~ CA1 <-> DA1	(k1, k_1)
	~ CA2 <-> DA2	(k2, k_2)
	~ DA1 <-> DA2     (Ka, k_a)	       
	~ CA2 <-> O	(alpha, beta)
	~ O <-> DA2b    (k3, k_3)
    
	CONSERVE C0+CA1+CA2+DA1+DA2+DA2b+O = 1
}

