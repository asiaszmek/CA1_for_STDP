TITLE l-calcium channel
: l-type calcium channel


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	FARADAY = 96485 (coul)
	R = 8.3134 (joule/degC)
}

PARAMETER {
	v (mV)
	celsius= 34	(degC)
	gbar=0 (mho/cm2)
	ki=.001 (mM)
	cai = 100.e-6 (mM)
        cao = 2 (mM)
	c_2 = 2 (mM)
     	tfa = 5 
        alpha_th = -27.01 (mV)
	alpha_slope = 3.8 (mV)
	alpha_mult = 0.055 (1/ms-mV)
        beta_th = -63.01 (mV)
	beta_slope = 17 (mV)
	beta_mult = 0.94 (1/ms)
}


NEURON {
	SUFFIX cal
	USEION ca READ cai,cao WRITE ica
        RANGE gbar,cai, ica, gcal, ggk
        RANGE minf,taum
}

STATE {
	m
}

ASSIGNED {
	ica (mA/cm2)
        gcal (mho/cm2)
        minf
        taum  (ms)
}

INITIAL {
	rate(v)
	m = minf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gcal = gbar*m*h2(cai)    
	ica = gcal*ghk(v, cai, cao)*cao/c_2

}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
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




FUNCTION alpm(v(mV)) (1/ms) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = alpha_mult*(alpha_th - v)/(exp((alpha_th-v)/alpha_slope) - 1)
}


FUNCTION betm(v(mV)) (1/ms){
        TABLE FROM -150 TO 150 WITH 200
        betm = beta_mult*exp((beta_th-v)/beta_slope)
}



DERIVATIVE state {  
        rate(v)
        m' = (minf - m)/taum
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a, b, qt
        a = alpm(v)
        taum = 1/(tfa*(a+betm(v))) : estimation of activation tau
        minf = a/(a+betm(v))       : estimation of activation steady state value

	
}
