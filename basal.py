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
    where_spines = False
    where_ca = ["lm_medium"]
    t_stop = 10000
    h.CVode()
    h.CVode().active(True)
    h.finitialize(-70)
    h.fcurrent()
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=[])
    h.CVode().re_init()
    h.celsius = 34
    h.finitialize(-70)
    h.fcurrent()
    t = h.Vector().record(h._ref_t, 10)
    sec_name = list(cell.shells.keys())[0]
    rec_v_soma = h.Vector().record(cell.soma[0](0.5)._ref_v, 10)
    ca_apic_0 = uf.record_specie_vec(cell.ca, cell.shells[sec_name][0], 10)
    ca_apic_1 = uf.record_specie_vec(cell.ca, cell.shells[sec_name][1], 10)
    ca_apic_2 = uf.record_specie_vec(cell.ca, cell.shells[sec_name][2], 10)
    ca_ecs = uf.record_specie_vec(cell.ca, cell.ECS, 10)
    camn_apic = uf.record_specie_vec(cell.camn, cell.shells[sec_name][0], 10)
    camc_apic = uf.record_specie_vec(cell.camc, cell.shells[sec_name][0], 10)
    calbca_apic = uf.record_specie_vec(cell.calbca, cell.shells[sec_name][0], 10)
    fixed_apic = uf.record_specie_vec(cell.fixedca, cell.shells[sec_name][0], 10)
    pmcaca_apic = uf.record_specie_vec(cell.pmcaca, cell.membrane[sec_name], 10)
    ncxca_apic = uf.record_specie_vec(cell.ncxca, cell.membrane[sec_name], 10)

    
    start = time.time()
    h.tstop = t_stop
    h.run(t_stop)
    print(time.time() - start)

    ca_new_0 = np.array([c.as_numpy()  for c in ca_apic_0])
    ca_new_1 = np.array([c.as_numpy()  for c in ca_apic_1])
    ca_new_2 = np.array([c.as_numpy()  for c in ca_apic_2])

    ca_ecs_new = np.array([c.as_numpy()  for c in ca_ecs])
    camn_new = np.array([c.as_numpy()  for c in camn_apic])
    camc_new = np.array([c.as_numpy()  for c in camc_apic])
    calbca_new = np.array([c.as_numpy()  for c in calbca_apic])
    fixed_new = np.array([c.as_numpy()  for c in fixed_apic])

    pmcaca_new = np.array([c.as_numpy()  for c in pmcaca_apic])
    ncxca_new = np.array([c.as_numpy()  for c in ncxca_apic])


    ca_new_2 = np.array([c.as_numpy()  for c in ca_apic_2])
    
    fig7, ax7 = plt.subplots(1, 1)
    ax7.plot(t.as_numpy(), 1e6*ca_new_0.mean(axis=0), label="Shell 0")
    ax7.plot(t.as_numpy(), 1e6*ca_new_1.mean(axis=0), label="Shell 1")
    ax7.plot(t.as_numpy(), 1e6*ca_new_2.mean(axis=0), label="Shell 2")
    ax7.legend()
    ax7.set_xlabel("time (ms)")
    ax7.set_ylabel("Ca in the dendite (nM)")
    fig7.savefig("Ca_cytosol_dendrite_plot.png")
    fig8, ax8 = plt.subplots(1, 1)
    ax8.plot(t.as_numpy(), rec_v_soma.as_numpy(), label="apic[23]")
    ax8.legend()
    ax8.set_xlabel("time (ms)")
    ax8.set_ylabel("Voltage (mM)")
    fig8.savefig("Voltage_in_the_soma.png")

    fig1, ax1 = plt.subplots(1, 1)
    show1 = ax1.imshow(pmcaca_new, extent=[0, t_stop, 0, 19],
                       origin="lower", interpolation="none", aspect="auto")
    fig1.colorbar(show1)
    ax1.set_title("pmcaCa membrane")
    ax1.set_xlabel("time (ms)")
    ax1.set_ylabel("dendrite (mm)")
    fig1.savefig("pmcaca_membrane_dend.png")
    fig3, ax3 = plt.subplots(1, 1)
    show3 = ax3.imshow(ncxca_new, extent=[0, t_stop, 0, 19], origin="lower",
                       interpolation="none", aspect="auto")
    fig3.colorbar(show3)
    ax3.set_title("ncxCa membrane")
    ax3.set_xlabel("time (ms)")
    ax3.set_ylabel("dendrite (mm)")
    fig3.savefig("ncxca_membrane.png")

    fig5, ax5 = plt.subplots(1, 1)
    show5 = ax5.imshow(fixed_new, extent=[0, t_stop, 0, 19], origin="lower",
                       interpolation="none", aspect="auto")
    fig5.colorbar(show5)
    ax5.set_title("Ca bound fixed")
    ax5.set_xlabel("time (ms)")
    ax5.set_ylabel("dendrite (mm)")
    fig5.savefig("fixed_bound_dendrite.png")

    fig7, ax7 = plt.subplots(3, 1)
    mini = min(ca_new_0.min(), ca_new_1.min(), ca_new_2.min())
    maxi = max(ca_new_0.max(), ca_new_1.max(), ca_new_2.max())
        
    show11 = ax7[0].imshow(ca_new_0, extent=[0, t_stop, 0, 19], origin="lower",
                           interpolation="none", aspect="auto", vmin=mini,
                           vmax=maxi)
    ax7[0].set_xticks([])
    show7 = ax7[1].imshow(ca_new_1, extent=[0, t_stop, 0, 19], origin="lower",
                          interpolation="none", aspect="auto", vmin=mini,
                          vmax=maxi)
    ax7[1].set_xticks([])
    show7 = ax7[2].imshow(ca_new_2, extent=[0, t_stop, 0, 19], origin="lower",
                          interpolation="none", aspect="auto", vmin=mini,
                          vmax=maxi)
    fig7.colorbar(show11, ax=ax7, shrink=0.6)
    fig7.suptitle("Ca cytosol")
    ax7[2].set_xlabel("time (ms)")
    ax7[2].set_ylabel("dendrite (mm)")

    fig7.savefig("Ca_cytosol.png")
    
    fig8, ax8 = plt.subplots(1, 1)
    show8 = ax8.imshow(ca_ecs_new, extent=[0, t_stop, 0, 19], origin="lower",
                       interpolation="none", aspect="auto")
    fig8.colorbar(show8)
    ax8.set_title("extracellular Ca")
    ax8.set_xlabel("time (ms)")
    ax8.set_ylabel("dendrite (mm)")
    fig8.savefig("extracellular_Ca_dendrite.png")
    
    fig9, ax9 = plt.subplots(1, 1)
    show9 = ax9.imshow(camn_new, extent=[0, t_stop, 0, 19], origin="lower",
                       interpolation="none", aspect="auto")
    fig9.colorbar(show9)
    ax9.set_title("CalmodulinNCa cytosol")
    ax9.set_xlabel("time (ms)")
    ax9.set_ylabel("dendrite (mm)")
    fig9.savefig("CaMN.png")


    fig10, ax10 = plt.subplots(1, 1)
    show10 = ax10.imshow(camc_new, extent=[0, t_stop, 0, 19], origin="lower",
                       interpolation="none", aspect="auto")
    fig10.colorbar(show10)
    ax10.set_title("CalmodulinCCa cytosol")
    ax10.set_xlabel("time (ms)")
    ax10.set_ylabel("dendrite (mm)")
    fig10.savefig("CaMC.png")
    
    fig11, ax11 = plt.subplots(1, 1)
    show11 = ax11.imshow(calbca_new, extent=[0, t_stop, 0, 19],
                         origin="lower", interpolation="none", aspect="auto")
    fig11.colorbar(show11)
    ax11.set_title("Ca bound Calbindin cytosol")
    ax11.set_xlabel("time (ms)")
    ax11.set_ylabel("dendrite (mm)")
    fig11.savefig("CalbindingCa.png")
    
    plt.show()
