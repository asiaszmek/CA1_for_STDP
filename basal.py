import time
import numpy as np
import matplotlib.pyplot as plt

from CA1 import CA1_PC
import utility_functions as uf
from neuron import h

if __name__ == "__main__":
    add_ER = False
    where_spines = []#["radTprox1"]
    where_ca = ["soma", "apical"]
    t_stop = 5000
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines)
    h.CVode().re_init()
    section_order = [sec.name() for sec in cell.sections_rxd
                     if "head" not in sec.name()]
    first_shells = [cell.shells[name][0] for name in section_order]

    t = h.Vector().record(h._ref_t, 100)
    ca_apic = {}
    for idx, sec_name in enumerate(section_order):
        shell = first_shells[idx]
        ca_apic[sec_name], header = uf.record_specie_vec(cell.ca,
                                                         shell, 100)
        
        
    start = time.time()
    cell.make_a_run(t_stop)
    print(time.time() - start)
    ca_new = {}
    for sec_name in ca_apic.keys():
        ca_new[sec_name] = np.array([c.as_numpy() for c in ca_apic[sec_name]])

    for sec in ca_new.keys():
        fig7, ax7 = plt.subplots(1, 1)
        for i in range(ca_new[sec].shape[0]):
            ax7.plot(t.as_numpy(), 1e6*ca_new[sec][i, :],
                     label="segment %d" %i)
        ax7.legend()
        ax7.set_xlabel("time (ms)")
        ax7.set_ylabel("Ca in the dendite (nM)")
        ax7.set_title(sec)
        fig7.savefig("%s_Ca_cytosol_dendrite.png" % sec)

    plt.show()
