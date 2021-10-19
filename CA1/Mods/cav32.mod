TITLE T-type calcium current (Cav3.2)

COMMENT
This is translated from Blackwell's moose_nerp, taking into accout that moose's
units are V and sec
    
ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = 96485 (coul)
    R = 8.3134 (joule/degC)
}

NEURON {
    THREADSAFE
    SUFFIX cav32
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE gbar, ica
    RANGE qfactCaT
}

PARAMETER {
    gbar = 0.0 (cm/s)
    a = 0.17
    Z (/mV)
    mmin = 0.0
    mmax = 1
    mvhalf = -52.3 (mV)
    mvslope = -6.8 (mV)
    mtmin = 0.916 (ms)
    mtmax = 13.3 (ms)
    mtvhalf = -53.6 (mV)
    mtvslope = 10.4 (mV)
    hmin = 0.
    hmax = 1.
    hvhalf = -73.8 (mV)
    hvslope = 2.8 (mV)
    htmin = 20.4 (ms)
    htmax = 45.6 (ms)
    htvhalf = -45.5 (mV)
    htvslope = 3 (mV)
    qfactCaT = 1
} 

ASSIGNED { 
    v (mV)
    ica (mA/cm2)
    celsius (degC)
    cai (mM)
    cao (mM)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = gbar*m^3*h*ghk(v, cai, cao)
}

INITIAL {
    Z = (0.001)*2*FARADAY/(R*(celsius+273.15 (degC)))
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states { 
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
}

PROCEDURE rates(v (mV)) {
    minf = mmin + mmax / (1 + exp((v + mvhalf) / mvslope))
    mtau = (mtmin + mtmax / (1 + exp((v + mtvhalf) / mtvslope)))/qfactCaT
    hinf = hmin + hmax / (1 + exp((v + hvhalf) / hvslope))
    htau = (htmin + htmax / (1 + exp((v + htvhalf) / htvslope)))/qfactCaT

}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = Z*v
    if(z == 0) {
        z = z+1e-6
    }
    eco = co*(z)/(exp(z)-1)
    eci = ci*(-z)/(exp(-z)-1)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}

