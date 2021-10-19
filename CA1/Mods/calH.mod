TITLE Ca L-type channel with high treshold of activation
: inserted in distal dendrites to account for distally
: restricted initiation of Ca++ spikes
: uses channel conductance (not permeability)
: written by Yiota Poirazi, 1/8/00 poirazi@LNC.usc.edu

NEURON {
	SUFFIX calH
	USEION ca READ cai, cao WRITE ica
        RANGE gbar, cai, ica, m, h
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)			
	FARADAY = 96485 (coul)
	R = 8.3134 (joule/degC)
		
}



PARAMETER {          : parameters that can be entered when function is called in cell-setup
        v               (mV)
        celsius = 34	(degC)
	dt              (ms)
        gbar = 0     (mho/cm2) : initialized conductance
	cai = 50.e-6 (mM)
	cao = 2 (mM)				
	t_act = 3.6 (ms)
        t_inact = 29 (ms)
	th_act = -37 (mV)
	slope_act = -1 (mV)
	th_inact = -41 (mV)
	slope_inact = 0.5 (mV)
        }

STATE {	m h }                     : unknown activation and inactivation parameters to be solved in the DEs  

ASSIGNED {
	ica (mA/cm2)
        inf[2]
	tau[2] (ms)

        
}


INITIAL {
      m = 0    : initial activation parameter value
	h = 1    : initial inactivation parameter value
	rate(v)
	
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ica = gbar*m*m*m*h*ghk(v, cai, cao)

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


DERIVATIVE state {  
        rate(v)
        m' = (inf[0]-m)/t_act
        h' = (inf[1]-h)/t_inact

}

PROCEDURE rate(v (mV)) { :callable from hoc
       FROM i=0 TO 1 {
		inf[i] = varss(v,i)
	}
   

}


FUNCTION varss(v (mV), i) {
	if (i==0) { 
             varss = 1 / (1 + exp((v-th_act)/slope_act))  : Ca activation 
	}
	else if (i==1) { 
             varss = 1 / (1 + exp((v-th_inact)/slope_inact)) : Ca inactivation 
	}
}


















