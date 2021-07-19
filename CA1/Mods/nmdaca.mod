: "Standard" NMDA model
: Dual exponential, with voltage-dependent max conductance
: A fraction of the current is carried by calcium
: BPG 22-8-11

TITLE nmda synapse 

NEURON {
	POINT_PROCESS NMDAca
	NONSPECIFIC_CURRENT i
	USEION ca READ cai, cao WRITE  ica
        RANGE g,a,b,gNMDAmax,tcon,tcoff, mgconc, eta, gamma, ggk
        RANGE fCa,i,ica,cai, Erev, c_2
}

UNITS {
        (uS) = (microsiemens)
        (nA) = (nanoamp)
        (mV) = (millivolt)
	(mM) = (milli/liter)			
	FARADAY = 96485 (coul)
	R = 8.3134 (joule/degC)
}

PARAMETER {
	fCa = 1.0	: fraction of current is calcium (BPG)
	tcon = 3 (ms)
	tcoff = 40 (ms): 150 (ms)
        gNMDAmax = 1	(uS)
	Erev	= 0    (mV)	: reversal potential					
	mgconc = 1	(mM)	: magnesium concentration
	eta = 0.33	(/mM)
	gamma = 0.06	(/mV)
	cai = 100.e-6 (mM)
        cao = 2 (mM)
	c_2 = 2 (mM)
	celsius = 34	(degC)
}

ASSIGNED {
	v 	(mV)
	i	(nA)
	ica	(nA)
	g       
	factor
}

INITIAL { 
   LOCAL tp
   a=0  
   b=0 
:   factor=tcon*tcoff/(tcoff-tcon)
   tp = (tcon*tcoff)/(tcoff - tcon) * log(tcoff/tcon)
   factor = -exp(-tp/tcon) + exp(-tp/tcoff)
   factor = 1/factor
}

STATE {
      a
      b
}

BREAKPOINT {
	LOCAL s
	SOLVE states METHOD derivimplicit
	s = 1.0/(1+eta*mgconc*exp(-gamma*v))
        g = b-a
	i = (1-fCa)*gNMDAmax*g*s*(v-Erev)
	ica = fCa*gNMDAmax*g*s*ghk(v, cai, cao)*cao/c_2
}

DERIVATIVE states {
	a' = -a/tcon
	b' = -b/tcoff
}

NET_RECEIVE(wgt) {
        LOCAL x
	x=wgt*factor
	state_discontinuity(a,a+x)
	state_discontinuity(b,b+x)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f
        f = KTF(celsius)
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
     KTF = (1e+3)*(R*(celsius + 273.15)/FARADAY/2)
     :KTF = RT/(zF)
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}


