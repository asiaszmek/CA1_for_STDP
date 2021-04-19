import time
import numpy as np
import matplotlib.pyplot as plt

from CA1 import CA1_PC
import utility_functions as uf
from neuron import h

def advance_a_bit(t_stop):
    for i in range(t_stop):
        h.fadvance()


if __name__ == "__main__":
    add_ER = False
    where_spines = []  # ["CA1_PC_Tomko[0].radTprox1"]
    where_ca = ["apical"]
    t_stop = 10000
    h.CVode()
    h.CVode().active(True)
    h.finitialize(-70)
    h.fcurrent()
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines)
    h.CVode().re_init()
    h.celsius = 34
    h.finitialize(-70)
    h.fcurrent()
    t = h.Vector().record(h._ref_t, 10)
    sec_name = "CA1_PC_Tomko[0].lm_medium2"
    print(sec_name)
    ca_apic = {}
    for sec in cell.apical:
        ca_apic[sec.name()] = []
        for shell in cell.shells[sec.name()]:
            ca_apic[sec.name()].append(uf.record_specie_vec(cell.ca,
                                                            shell, 10))
    
    
    ca_ecs = uf.record_specie_vec(cell.ca, cell.ECS, 10)
    start = time.time()
    h.tstop = t_stop
    h.run(t_stop)
    print(time.time() - start)
    ca_new = {}
    for sec_name in ca_apic.keys():
        ca_new[sec_name] = []
        for shell in ca_apic[sec_name]:
            ca_new[sec_name].append(np.array([c.as_numpy() for c in shell]))

    ca_ecs_new = np.array([c.as_numpy()  for c in ca_ecs])
    for sec in ca_new.keys():
                                    
        fig7, ax7 = plt.subplots(1, 1)
        for i, ca_shell in enumerate(ca_new[sec]):
            ax7.plot(t.as_numpy(), 1e6*ca_shell.mean(axis=0),
                     label="Shell %d" %i)
        ax7.legend()
        ax7.set_xlabel("time (ms)")
        ax7.set_ylabel("Ca in the dendite (nM)")
        ax7.set_title(sec)
        fig7.savefig("%s_Ca_cytosol_dendrite.png" % sec)

    fig8, ax8 = plt.subplots(1, 1)
    show8 = ax8.imshow(ca_ecs_new, extent=[0, t_stop, 0, 19], origin="lower",
                       interpolation="none", aspect="auto")
    fig8.colorbar(show8)
    ax8.set_title("extracellular Ca")
    ax8.set_xlabel("time (ms)")
    ax8.set_ylabel("dendrite (mm)")
    fig8.savefig("extracellular_Ca_dendrite.png")
    
    plt.show()
