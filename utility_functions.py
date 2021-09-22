from neuron import h
import numpy as np

def record_specie_vec(specie, region, record_dt, spine=False, special_name=""):
    '''use vec.record to record cyt calcium concentration. 

    vectors, one vector for each node'''
    Recordl = []
    header = ""
    for mynode in specie[region].nodes:
        if spine:
            if "head" in mynode.sec.name():
                if not special_name:
                    cytcavec = h.Vector()
                    cytcavec.record(mynode._ref_concentration, record_dt)
                    Recordl.append(cytcavec)
                    header += " " + "%s_%s" %(specie.name, mynode.segment)
                elif special_name in mynode.sec.name():
                    cytcavec = h.Vector()
                    cytcavec.record(mynode._ref_concentration, record_dt)
                    Recordl.append(cytcavec)
                    header += " " + "%s_%s" %(specie.name, mynode.segment)
        else:
            if "head" not in mynode.sec.name():
                if not special_name:
                    cytcavec = h.Vector()
                    cytcavec.record(mynode._ref_concentration, record_dt)
                    Recordl.append(cytcavec)
                    header += " " + "%s_%s_%s" %(specie.name, region.name,
                                                 mynode.segment)
                elif special_name in mynode.sec.name():
                    cytcavec = h.Vector()
                    cytcavec.record(mynode._ref_concentration, record_dt)
                    Recordl.append(cytcavec)
                    header += " " + "%s_%s_%s" %(specie.name, region.name,
                                                 mynode.segment)
    return Recordl, header



def add_bath(conc, t_start, t_stop, v_init): #in sec!
    pre = h.Section(name="soma")
    pre.L = 10
    pre.diam = 10
    pre.insert("pas")
    pre.g_pas = 1/5000
    pre.e_pas = v_init
    pre.insert("bath")
    pre.conc_bath = conc
    pre.v_init_bath = v_init
    iclamp = h.IClamp(pre(0.5))
    iclamp.delay = t_start*1000
    iclamp.dur = (t_stop-t_start)*1000
    iclamp.amp = 0.00063
    
    return pre, iclamp
    
def spine_average(indicator, spine):
    volumes = []
    for seg in spine:
        volumes.append(seg.volume())
    volumes = np.array(volumes)
    out = (indicator*volumes[:, None]).sum(axis=0)/volumes.sum()
    return out
    
def dend_average(indicator, dend, factors):
    volumes = np.zeros((len(factors), dend.nseg))
    seg_len = dend.L/dend.nseg
    for i, seg in enumerate(dend):
        radius = seg.diam/2
        prev_vol = seg.volume()
        for j, factor in enumerate(factors):
            radius = radius - factor*seg.diam/2
            volumes[j, i] = prev_vol - np.pi*seg_len*radius**2
            prev_vol = np.pi*seg_len*radius**2
    
    return (indicator*volumes[:, :, None]).sum(axis=0)/volumes.sum(axis=0)[:, None]

def calculate_T(t, T0, t_1, t_2, t_0):
    return T0*(-np.exp(-(t-t_0)/t_1)+np.exp(-(t-t_0)/t_2))

def calculate_R(t, R0, t_rec, t0):
    return 1 - (1 - R0)*np.exp(-(t-t0)/t_rec)

def add_transmitter(no_pulses, freq, stim_start, t_stop, dt=0.1,
                    t_rec=4167,
                    use=0.0044, t_in=0.5, t_off=0.1, ase=18):
    t_peak = t_off*t_in/(t_in - t_off) * np.log(t_in/t_off)
    factor = -np.exp(-t_peak/t_off) + np.exp(-t_peak/t_in)
    isi = 1000/freq
    t = np.arange(0, t_stop+dt, dt)
    T = np.zeros(len(t))
    R = np.zeros(len(t))
    for i in range(no_pulses):
        tstim = (stim_start + i*isi)
        tstamp = int((stim_start + i*isi)/dt)
        R_ = R[tstamp-1]
        R_new = (1-use)*R_
        T_new = ase*(use*R_)/factor
        R[tstamp:] = calculate_R(t[tstamp:], R_new, t_rec, tstim)
        T[tstamp:] = calculate_T(t[tstamp:], T_new, t_off,
                                 t_in, tstim)
    
    return t, T

def add_stim(syn, pairings, freq, start, w):
    stim = h.NetStim()
    stim.number = pairings
    stim.interval = 1000.0 / freq
    stim.start = start
    stim.noise = 0
    netcon = h.NetCon(stim, syn, 0, 0, w)
    return netcon, stim
