import time
import numpy as np
import matplotlib.pyplot as plt

from CA1 import CA1_PC
import utility_functions as uf
from neuron import h

if __name__ == "__main__":
    add_ER = False
    where_spines = []#["radTprox1"]
    where_ca = ["radTprox1", "radTprox2"]
    t_stop = 50000
    h.load_file("stdrun.hoc")
    h.CVode()
    h.CVode().active(True)
    h.CVode().re_init()
    h.celsius = 34
    h.finitialize(-65)
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines)
    h.CVode().re_init()
    section_order = [sec.name() for sec in cell.sections_rxd
                     if "head" not in sec.name()]
    first_shells = [cell.shells[name][0] for name in section_order]
    camn_apic = uf.record_specie_vec(cell.camn, first_shells, 100)
    camc_apic = uf.record_specie_vec(cell.camc, first_shells, 100)
    calbca_apic = uf.record_specie_vec(cell.calbca, first_shells, 100)
    fixed_apic = uf.record_specie_vec(cell.fixedca, first_shells, 100)
    t = h.Vector().record(h._ref_t, 100)
    ca_apic = {}
    for idx, sec_name in enumerate(section_order):
        shell = first_shells[idx]
        ca_apic[sec_name] = uf.record_specie_vec(cell.ca,
                                                 shell, 100)
    start = time.time()
    h.tstop = t_stop
    h.run(t_stop)
    print(time.time() - start)
    ca_new = {}
    for sec_name in ca_apic.keys():
        ca_new[sec_name] = np.array([c.as_numpy() for c in ca_apic[sec_name]])

    for sec in ca_new.keys():
        fig7, ax7 = plt.subplots(1, 1)
        print(sec, ca_new[sec].shape, ca_new[sec][:, -1])
        for i in range(ca_new[sec].shape[0]):
            ax7.plot(t.as_numpy(), 1e6*ca_new[sec][i, :],
                     label="segment %d" %i)
        ax7.legend()
        ax7.set_xlabel("time (ms)")
        ax7.set_ylabel("Ca in the dendite (nM)")
        ax7.set_title(sec)
        fig7.savefig("%s_Ca_cytosol_dendrite.png" % sec)

    camn_new = np.array([c.as_numpy()  for c in camn_apic])
    camc_new = np.array([c.as_numpy()  for c in camc_apic])
    calbca_new = np.array([c.as_numpy()  for c in calbca_apic])
    fixed_new = np.array([c.as_numpy()  for c in fixed_apic])



    fig5, ax5 = plt.subplots(1, 1)
    show5 = ax5.imshow(fixed_new, extent=[0, t_stop, 0, 10000], origin="lower",
                       interpolation="none", aspect="auto")
    fig5.colorbar(show5)
    ax5.set_title("Ca bound fixed")
    ax5.set_xlabel("time (ms)")
    ax5.set_ylabel("dendrite (mm)")
    fig5.savefig("fixed_bound_dendrite.png")

    fig9, ax9 = plt.subplots(1, 1)
    show9 = ax9.imshow(camn_new, extent=[0, t_stop, 0, 10000], origin="lower",
                       interpolation="none", aspect="auto")
    fig9.colorbar(show9)
    ax9.set_title("CalmodulinNCa cytosol")
    ax9.set_xlabel("time (ms)")
    ax9.set_ylabel("dendrite (mm)")
    fig9.savefig("CaMN.png")


    fig10, ax10 = plt.subplots(1, 1)
    show10 = ax10.imshow(camc_new, extent=[0, t_stop, 0, 10000], origin="lower",
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
