import numpy as np
import matplotlib.pyplot as plt
import neuron
    
if __name__ == "__main__":
    neuron.h.nrn_load_dll('./Mods/nrnmech.dll')
    neuron.h.xopen('pyramidal_cell_weak_bAP_original.hoc')
    cell = neuron.h.CA1_PC_Tomko()
    neuron.h.CVode()
    neuron.h.CVode().active(True)
    fig, ax = plt.subplots(1, 1)
    injections = [-.2, 0, .200, .400, 0.5]
    neuron.h.celsius = 34
    soma = cell.soma[0]
    print(soma)
    for inj in injections:
        print(inj)
        ic = neuron.h.IClamp(soma(0.5))
        rec_v = neuron.h.Vector().record(soma(0.5)._ref_v)
        time = neuron.h.Vector().record(neuron.h._ref_t)
        ic.delay = 300
        ic.dur = 1000
        ic.amp = inj
        neuron.h.finitialize(-70)
        neuron.h.fcurrent()
        neuron.run(1500)
        ax.plot(time, rec_v, label="%4.3fnA" % inj)


 
    ax.set_xlabel("time (s)")
    ax.set_ylabel("V (mV)")
    ax.legend()
    plt.show()
