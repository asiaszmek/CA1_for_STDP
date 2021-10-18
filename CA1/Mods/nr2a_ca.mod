TITLE detailed kinetic model of NMDA receptors
:Ausra Saudargiene 2020

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of NMDA receptors
    ===============================
  
                 RA2d1 
                 |
    R  -- RA  -- RA2  --  RA2f, RA2s --O 
                 |
                 RA2d2



                 RA2d1Mg 
                 |
    RMg  -- RAMg  -- RA2Mg  --  RA2fMg, RA2sMg --OMg 
                 |
                 RA2d2Mg


    Voltage dependence of Mg2+ block:
    Ascher and Nowak 1988
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------

M. Vargas Caballero /  H. Robinson 2002

-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS NR2A_CA
    POINTER T 
    NONSPECIFIC_CURRENT i
    USEION ca READ cai, cao WRITE ica
    RANGE state_R, state_RA, state_RA2, state_RA2d1,state_RA2d2,state_RA2f, state_RA2s, O
    RANGE state_RMg, state_RAMg, state_RA2Mg, state_RA2d1Mg,state_RA2d2Mg,state_RA2fMg, state_RA2sMg, OMg, glu
    RANGE gmax, mg_on, mg_off, rb, rb_Mg, fCa
    RANGE kd1_plus, kd1_minus, kd2_plus, kd2_minus, kf_plus, ks_plus0, ks_plus, kf_minus, ks_minus, kon, koff, fCa, cai, ggk, Erev, c_2
    GLOBAL mg
    RANGE vmin, vmax

}

UNITS {
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (mM) = (milli/liter)
    (nA) = (nanoamp)
    (uM) = (micro/liter)
    FARADAY = 96485 (coul)
    R = 8.3134 (joule/degC)
}

PARAMETER {
    fCa = 0.06	: fraction of current is calcium (BPG)
    Erev    = -3    (mV) : reversal potential
    mg  = 1     (mM)     : external magnesium concentration
    vmin = -180 (mV)
    vmax = 100  (mV)
    : gmax=0.05 (mho/cm2)  : conductance of NMDA channels
    gmax=5 (uS)  : conductance of NMDA channels
    cai = 100.e-6 (mM)   : initial internal Ca++ concentration
    cao = 2       (mM)
    c_2 = 2 (mM)
    
   
    :Erreger et al, 2005
    :NR2A_CA subunit
    kon  = 31.6e-3    (/uM /ms)    : binding       
    koff  = 1010e-3   (/ms)      : unbinding     
    
    kd1_plus=85.1e-3   (/ms)
    kd1_minus=29.7e-3 (/ms)
    kd2_plus=230e-3   (/ms)
    kd2_minus=1.01e-3 (/ms)

    ks_plus0=230e-3     (/ms) :Table 2
    ks_minus=178e-3   (/ms)
    kf_plus=3140e-3   (/ms)  
    kf_minus=174e-3   (/ms)
   celsius = 34	(degC)
   mg_on_rate = 610 (/ms)
    mg_off_rate = 5 (/ms)
}

ASSIGNED { 
    v       (mV)        : postsynaptic voltage
    rb      (/ms)       : binding
    rb_mg   (/ms)   : binding of glutamate to blocked channels
    mg_on   (/ms)       : blocking rate
    mg_off  (/ms)       : unblocking rate
    glu (mM)   : glutamate from release
    T (mM)   : glutamate from release
 
   i     (nA) :(mA/cm2)  : Na current density
    ica     (nA) :(mA/cm2)  : Ca current density
    ks_plus (/ms)     : rate 
 
}


STATE {
    : Channel states (all fractions)

   state_R
   state_RA
   state_RA2
   state_RA2d1
   state_RA2d2
   state_RA2f
   state_RA2s
   O
   state_RMg
   state_RAMg
   state_RA2Mg
   state_RA2d1Mg
   state_RA2d2Mg
   state_RA2fMg
   state_RA2sMg
   OMg
   

}

INITIAL {
    rates(v)
:   state_C0 = 1
:   Obidos: state_C0 - Mg free unbound
    state_RMg=1
    mg = 1 (mM)
}

BREAKPOINT {
     :transmitter()
    SOLVE kstates METHOD sparse
    i = (1-fCa)* gmax*O * (v - Erev)   : current in mA/cm2 
    ica = fCa* gmax*O * ghk(v, cai, cao) *cao/c_2  : current in mA/cm2
}


KINETIC kstates {
     rates(v)
     if (T < 0) {glu = 0}
     else {glu = T}	  
	  
     
    rb = (1e+3)*(kon * glu) 
    rb_mg =(1e+3)*(kon * glu)
    ks_plus=ks_plus0*exp((v+100 (mV))/175 (mV)) :Vm - in mV, Vdep=175mV

   ~ state_R   <-> state_RA     (2*rb,koff)
   ~ state_RA  <-> state_RA2    (rb,2*koff)
   ~ state_RA2 <-> state_RA2d1  (kd1_minus,kd1_plus)
   ~ state_RA2 <-> state_RA2d2  (kd2_plus,kd2_minus)
   ~ state_RA2f <-> state_RA2   (kf_minus,kf_plus)
   ~ state_RA2s <-> state_RA2   (ks_minus,ks_plus)
   ~ O <-> state_RA2f   (ks_minus,ks_plus)
   ~ O <-> state_RA2s   (kf_minus,kf_plus)
   ~ O <-> OMg   (mg_on, mg_off) 
   ~ state_RMg  <-> state_RAMg    (2*rb,koff)
   ~ state_RAMg  <-> state_RA2Mg    (rb,2*koff)
   ~ state_RA2Mg <-> state_RA2d1Mg  (kd1_minus,kd1_plus)
   ~ state_RA2Mg <-> state_RA2d2Mg  (kd2_plus,kd2_minus)
   ~ state_RA2fMg <-> state_RA2Mg   (kf_minus,kf_plus)
   ~ state_RA2sMg <-> state_RA2Mg   (ks_minus,ks_plus)
   ~ OMg <-> state_RA2fMg   (ks_minus,ks_plus)
   ~ OMg <-> state_RA2sMg   (kf_minus,kf_plus)

    
    CONSERVE state_R+state_RA+state_RA2+state_RA2d1+state_RA2d2+state_RA2f+ state_RA2s+ O+state_RMg+state_RAMg+state_RA2Mg+state_RA2d1Mg+state_RA2d2Mg+state_RA2fMg+ state_RA2sMg+ OMg = 1
}


FUNCTION KTF(celsius (degC)) (mV) {
     KTF = (1e+3)*(R*(celsius + 273.15)/FARADAY/2)
     :KTF = RT/(zF)
}
 
FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f
        f = KTF(celsius)
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}



FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}




PROCEDURE rates(v(mV)) {
    :Antonov and Johnson, 2008
:   mg_on = 1.1e3*exp(-v/55) :  (/mM /ms)    
 :  mg_off = 110*exp(v/52.7): (/ms)



    : from Ascher and Nowak - mg replaced with 1
    mg_on = mg_on_rate*exp(-v/17(mV))
    mg_off = mg_off_rate*exp(v/47(mV))
}

: Use procedure only when working in non-stationary conditions, 
: when in use, uncomment the call for this
: procedure at the breakpoint

: a brief square pulse may reproduce the synaptic activation
: of the receptor



