TITLE HVA L-type calcium current (Cav1.2)

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
    SUFFIX cal12
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE gbar, ica
    RANGE m_A, m_B
    RANGE qfactCaL
}

PARAMETER {
    gbar = 0.0 (cm/s)
    a = 0.17
    Z (/mV)
    hmin = 0.83
    hmax = 0.17
    hvhalf = -55 (mV)
    hvslope = 8 (mV)
    A_A = -0.880066 (/ms)
    A_B = -0.220 (/ms-mV)
    A_C = -1
    A_vhalf = 4 (mV)
    A_vslope = -8 (mV)
    B_A = -0.2840213 (/ms)
    B_B = 71e-3 (/mV-ms)
    B_C = -1
    B_vhalf = -4 (mV)
    B_vslope = 5 (mV)

    htau_in = 44.3 (ms)
    qfactCaL = 2
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
    m_A (/ms)
    m_B (/ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = gbar*m*h*ghk(v, cai, cao)
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
    hinf = hmin + hmax / (1 + exp((v + hvhalf) / hvslope))
    htau = htau_in/qfactCaL

    m_A = (A_A+A_B*v)/(A_C+exp((v+A_vhalf)/A_vslope))*qfactCaL
    m_B = (B_A+B_B*v)/(B_C+exp((v+B_vhalf)/B_vslope))*qfactCaL
    minf = m_A/(m_A + m_B)
    mtau = 1/(m_A + m_B)

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

