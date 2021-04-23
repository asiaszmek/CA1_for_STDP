import time
import numpy as np
import matplotlib.pyplot as plt

from CA1 import CA1_PC
import utility_functions as uf
from neuron import h

if __name__ == "__main__":
    add_ER = False
    where_spines = ["radTprox1"]
    where_ca = ["apical"]
    t_stop = 100
    h.load_file("stdrun.hoc")
    h.CVode()
    h.CVode().active(True)
    h.CVode().re_init()
    h.celsius = 34
    h.finitialize(-65)
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines)
    h.CVode().re_init()
    section_order = [sec.name() for sec in cell.sections_rxd]
    first_shells = [cell.shells[name][0] for name in section_order
                    if "head" not in name]
    camn_apic = uf.record_specie_vec(cell.camn, first_shells, 10)
    camc_apic = uf.record_specie_vec(cell.camc, first_shells, 10)
    calbca_apic = uf.record_specie_vec(cell.calbca, first_shells, 10)
    fixed_apic = uf.record_specie_vec(cell.fixedca, first_shells, 10)
    pmcaca_apic = uf.record_specie_vec(cell.pmca, cell.membrane_list, 10)
    ncxca_apic = uf.record_specie_vec(cell.ncx, cell.membrane_list, 10)
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
    camn_new = np.array([c.as_numpy()  for c in camn_apic])
    camc_new = np.array([c.as_numpy()  for c in camc_apic])
    calbca_new = np.array([c.as_numpy()  for c in calbca_apic])
    fixed_new = np.array([c.as_numpy()  for c in fixed_apic])

    pmcaca_new = np.array([c.as_numpy()  for c in pmcaca_apic])
    ncxca_new = np.array([c.as_numpy()  for c in ncxca_apic])

    fig1, ax1 = plt.subplots(1, 1)
    show1 = ax1.imshow(pmcaca_new, extent=[0, t_stop, 0, 1000],
                       origin="lower", interpolation="none", aspect="auto")
    fig1.colorbar(show1)
    ax1.set_title("pmcaCa membrane")
    ax1.set_xlabel("time (ms)")
    ax1.set_ylabel("dendrite (mm)")
    fig1.savefig("pmcaca_membrane_dend.png")
    fig3, ax3 = plt.subplots(1, 1)
    show3 = ax3.imshow(ncxca_new, extent=[0, t_stop, 0, 1000], origin="lower",
                       interpolation="none", aspect="auto")
    fig3.colorbar(show3)
    ax3.set_title("ncxCa membrane")
    ax3.set_xlabel("time (ms)")
    ax3.set_ylabel("dendrite (mm)")
    fig3.savefig("ncxca_membrane.png")

    fig5, ax5 = plt.subplots(1, 1)
    show5 = ax5.imshow(fixed_new, extent=[0, t_stop, 0, 1000], origin="lower",
                       interpolation="none", aspect="auto")
    fig5.colorbar(show5)
    ax5.set_title("Ca bound fixed")
    ax5.set_xlabel("time (ms)")
    ax5.set_ylabel("dendrite (mm)")
    fig5.savefig("fixed_bound_dendrite.png")

    fig9, ax9 = plt.subplots(1, 1)
    show9 = ax9.imshow(camn_new, extent=[0, t_stop, 0, 1000], origin="lower",
                       interpolation="none", aspect="auto")
    fig9.colorbar(show9)
    ax9.set_title("CalmodulinNCa cytosol")
    ax9.set_xlabel("time (ms)")
    ax9.set_ylabel("dendrite (mm)")
    fig9.savefig("CaMN.png")


    fig10, ax10 = plt.subplots(1, 1)
    show10 = ax10.imshow(camc_new, extent=[0, t_stop, 0, 1000], origin="lower",
                       interpolation="none", aspect="auto")
    fig10.colorbar(show10)
    ax10.set_title("CalmodulinCCa cytosol")
    ax10.set_xlabel("time (ms)")
    ax10.set_ylabel("dendrite (mm)")
    fig10.savefig("CaMC.png")
    
    fig11, ax11 = plt.subplots(1, 1)
    show11 = ax11.imshow(calbca_new, extent=[0, t_stop, 0, 1000],
                         origin="lower", interpolation="none", aspect="auto")
    fig11.colorbar(show11)
    ax11.set_title("Ca bound Calbindin cytosol")
    ax11.set_xlabel("time (ms)")
    ax11.set_ylabel("dendrite (mm)")
    fig11.savefig("CalbindingCa.png")

    plt.show()
