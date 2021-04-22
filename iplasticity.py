import numpy as np
import matplotlib.pyplot as plt
from CA1 import CA1_PC
import neuron
    
if __name__ == "__main__":
    where_ca = []
    cell = CA1_PC(where_ca=where_ca)
    neuron.h.CVode()
    neuron.h.CVode().active(True)
    fig, ax = plt.subplots(1, 1)
    injections = [-.2, 0, .200, .400, 0.5]
    neuron.h.celsius = 35
    soma = cell.soma
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


    f = open("python_model_no_calc", "w")
    for sec in cell.apical:
        f.write("%s\n" % sec.name())
        for key in sorted(sec.psection()):
            if isinstance(sec.psection()[key], dict):
                for key2 in sec.psection()[key]:
                    if isinstance(sec.psection()[key][key2], dict):
                        for key3 in sec.psection()[key][key2].keys():
                            f.write("%s, %s\n" % (key3,
                                                  sec.psection()[key][key2][key3]))
                    else:
                        f.write("%s, %s\n" %(key2, sec.psection()[key][key2]))
            elif isinstance(sec.psection()[key], list):
                for item in sec.psection()[key]:
                    f.write("%s, %s\n" % (key, item))
            elif isinstance(sec.psection()[key], set):
                for item in list(sec.psection()[key]):
                    f.write("%s, %s\n" % (key, item))
            else:
                f.write("%s, %s\n" % (key, sec.psection()[key]))
    
    ax.set_xlabel("time (s)")
    ax.set_ylabel("V (mV)")
    ax.legend()
    plt.show()
