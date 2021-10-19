TITLE K-DR channel
: from Klee Ficker and Heinemann
: modified to account for Dax et al.
: M.Migliore 1997

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = 96485 (coul)
	R = 8.3134 (joule/degC)
}

PARAMETER {
	v (mV)
        ek (mV)		: must be explicitely def. in hoc
	celsius		(degC)
	gbar=.003 (mho/cm2)
        vhalfn=13   (mV)
        a0n=0.02      (/ms)
        zetan=-1.5    
        gmn=0.7  (1)
	nmax=2  (ms)
	q10=1
	z (/mV)
}


NEURON {
	SUFFIX kdr
	USEION k READ ek WRITE ik
        RANGE gkdr, gbar
	RANGE ninf,taun
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkdr (mho/cm2)
        taun (ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gbar*n
	ik = gkdr*(v-ek)

}

INITIAL {
        z = (0.001)*2*FARADAY/(R*(celsius+273.15 (degC)))
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(zetan*(v-vhalfn)*z)
}

FUNCTION betn(v(mV)) {
  betn = exp(zetan*gmn*(v-vhalfn)*z)
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt = q10^((celsius-24 (degC))/10 (degC))
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
	if (taun<nmax) {taun=nmax}
}
