TITLE T-calcium channel
: T-type calcium channel


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER {
	v (mV)
	celsius = 25	(degC)
	gbar=.003 (mho/cm2)
	cai = 50.e-6 (mM)
	cao = 2 (mM)
	q10 = 5
	mmin=0.2 (ms)
	hmin=10 (ms)
	a0h =0.015
	zetah = 3.5 (/mV)
	vhalfh = -75 (mV)
	gmh=0.6	
	a0m =0.04
	zetam = 2 (/mV)
	vhalfm = -28 (mV)
        gmm=0.1
	tshift = 25 (degC)
        tslope = 10 (degC)
	A = 0.2 (/mV)
}


NEURON {
	SUFFIX cat
	USEION ca READ cai,cao WRITE ica
        RANGE gbar, ica, gcat
        RANGE hinf,minf,mtau,htau
}

STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        gcat (mho/cm2)
	hinf
	htau (ms)
	minf
	mtau (ms)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcat = gbar*m*m*h
	ica = gcat*ghk(v,cai,cao)

}

DERIVATIVE states {	: exact when v held constant
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
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

FUNCTION alph(v(mV)) {
  alph = exp(0.0378*zetah*(v-vhalfh)) 
}

FUNCTION beth(v(mV)) (ms) {
  beth = exp(0.0378*zetah*gmh*(v-vhalfh)) *1 (ms)
}

FUNCTION alpmt(v(mV)) {
  alpmt = exp(0.0378*zetam*(v-vhalfm)) 
}

FUNCTION betmt(v(mV)) (ms){
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm)) *1 (ms)
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a,b, qt
        qt=q10^((celsius-tshift)/tslope)

	a = A*(-1.0*v+19.26)/(exp((-1.0*v+19.26 (mV))/10.0 (mV))-1.0)
	b = 0.009*exp(-v/22.03 (mV))
	minf = a/(a+b)
	mtau = betmt(v)/(qt*a0m*(1+alpmt(v)))
	if (mtau<mmin) {mtau=mmin}

	a = 1.e-6*exp(-v/16.26 (mV))
	b = 1/(exp((-v+29.79 (mV))/10. (mV))+1.) 
	hinf = a/(a+b)
	htau = beth(v)/(a0h*(1+alph(v)))
	if (htau<hmin) {htau=hmin}
}

