TITLE R-type calcium current (CaR)

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
    SUFFIX car
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE gbar, ica
    RANGE qfactCaR
    RANGE m_A, m_B
    RANGE h_A, h_B
}

PARAMETER {
    gbar = 0.0 (cm/s)
    a = 0.17
    Z (/mV)
    mA_A = 0.24 (/ms)
    mA_B =  0 (/mV-ms)
    mA_C = 0 
    mA_vhalf = 0 (mV)
    mA_vslope = -28 (mV)
    mB_A = 1264(/ms)
    mB_B =  8(/mV-ms)
    mB_C = -1
    mB_vhalf = 158 (mV)
    mB_vslope = 13.6 (mV)
    hA_A = 1.10 (/ms)
    hA_B =  0.01(/mV-ms)
    hA_C = -1
    hA_vhalf = 110 (mV)
    hA_vslope = 17 (mV)
    hB_A = 20e-3(/ms)
    hB_B =  0(/mV-ms)
    hB_C = 0
    hB_vhalf = 0 (mV)
    hB_vslope = -30 (mV)
    qfactCaR = 2

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
    h_A (/ms)
    h_B (/ms)
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
    m_A = (mA_A + mA_B*v)/(mA_C+exp((v+mA_vhalf)/mA_vslope))*qfactCaR
    m_B = (mB_A + mB_B*v)/(mB_C+exp((v+mB_vhalf)/mB_vslope))*qfactCaR
    h_A = (hA_A + hA_B*v)/(hA_C+exp((v+hA_vhalf)/hA_vslope))
    h_B = (hB_A + hB_B*v)/(hB_C+exp((v+hB_vhalf)/hB_vslope))
				
    minf = m_A/(m_A+m_B)
    mtau = 1/(m_A+m_B)
    hinf = h_A/(h_A+h_B)
    htau = 1/(h_A+h_B)

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

