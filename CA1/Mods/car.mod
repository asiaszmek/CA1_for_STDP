TITLE Ca R-type channel with medium threshold for activation
: used in distal dendritic regions, together with calH.mod, to help
: the generation of Ca++ spikes in these regions
: uses channel conductance (not permeability)
: written by Yiota Poirazi on 11/13/00 poirazi@LNC.usc.edu
:
: updated to use CVode by Carl Gold 08/10/03
:  Updated by Maria Markaki  03/12/03

NEURON {
	SUFFIX car
	USEION ca READ cai, cao WRITE ica
        RANGE gbar, m, h, ica
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = 96485 (coul)
	R = 8.3134 (joule/degC)
}


ASSIGNED {               : parameters needed to solve DE
	ica (mA/cm2)
:	iCa (mA/cm2)
        inf[2]
	tau[2]		(ms)
        v               (mV)
        
	   

}


PARAMETER {              : parameters that can be entered when function is called in cell-setup
        gbar = 0      (mho/cm2) : initialized conductance
	cai = 1.e-4   (mM)       : initial internal Ca++ concentration
	cao = 2       (mM)       : initial external Ca++ concentration
        celsius = 34 	(degC)
	t_act = 50 (ms)
        t_inact = 5 (ms)
        th_act = -48.5 (mV)
	slope_act = -3 (mV)
        th_inact = -53 (mV)
	slope_inact = 1 (mV)
	
}  

STATE {	
	m 
	h 
}            : unknown activation and inactivation parameters to be solved in the DEs  


INITIAL {
	rates(v)
        m = 0    : initial activation parameter value
	h = 1    : initial inactivation parameter value
}

BREAKPOINT {
     SOLVE states METHOD cnexp
	   
	ica = gbar*m*m*m*h*ghk(v, cai, cao)

}


DERIVATIVE states {
	rates(v)
	m' = (inf[0]-m)/t_act
	h' = (inf[1]-h)/t_inact
}

PROCEDURE rates(v(mV)) {LOCAL a, b :rest = -70
	FROM i=0 TO 1 {
	inf[i] = varss(v,i)
	}
}




FUNCTION varss(v(mV), i) {
	if (i==0) {
	    varss = 1 / (1 + exp((v - th_act)/slope_act)) : Ca activation
	}
	else if (i==1) {
             varss = 1/ (1 + exp((v - th_inact)/slope_inact)) : Ca inactivation
	}
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














