import time
import numpy as np
import matplotlib.pyplot as plt

from CA1 import CA1_PC
import utility_functions as uf
from neuron import h

if __name__ == "__main__":
    add_ER = False
    where_spines = []  # ["CA1_PC_Tomko[0].radTprox1"]
    where_ca = ["apical"]
    t_stop = 10000
    h.load_file("stdrun.hoc")
    h.CVode()
    h.CVode().active(True)
    h.CVode().re_init()
    h.celsius = 34
    h.finitialize(-65)
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines)
    h.CVode().re_init()
    
    #f = open("python_model_3", "w")
    # for sec in cell.sections:
    #     f.write("%s\n" % sec.name())
    #     for key in sorted(sec.psection()):
    #         if isinstance(sec.psection()[key], dict):
    #             for key2 in sec.psection()[key]:
    #                 if isinstance(sec.psection()[key][key2], dict):
    #                     for key3 in sec.psection()[key][key2].keys():
    #                         f.write("%s, %s\n" % (key3,
    #                                               sec.psection()[key][key2][key3]))
    #                 else:
    #                     f.write("%s, %s\n" %(key2, sec.psection()[key][key2]))
    #         elif isinstance(sec.psection()[key], list):
    #             for item in sec.psection()[key]:
    #                 f.write("%s, %s\n" % (key, item))
    #         elif isinstance(sec.psection()[key], set):
    #             for item in list(sec.psection()[key]):
    #                 f.write("%s, %s\n" % (key, item))
    #         else:
    #             f.write("%s, %s\n" % (key, sec.psection()[key]))
    
    
    
    t = h.Vector().record(h._ref_t, 10)
    ca_apic = {}
    for sec in cell.sections:
        sec_name = sec.name()
        if sec_name in cell.shells:
            ca_apic[sec_name] = []
            for shell in cell.shells[sec_name]:
                ca_apic[sec_name].append(uf.record_specie_vec(cell.ca,
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
