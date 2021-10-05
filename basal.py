import time
import numpy as np
import matplotlib.pyplot as plt

from CA1 import CA1_PC
import utility_functions as uf
from neuron import h
dt = 10
if __name__ == "__main__":
    add_ER = False
    where_spines = []
    where_ca = ["soma", "apical"]
    t_stop = 10000
    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines,
                  buffer_list=["Calmodulin", "Calbindin", "Fixed"])
    section_order = [sec.name() for sec in cell.sections_rxd
                     if "head" not in sec.name()]
    
    sec_list = ["soma", "radTprox1", "radTprox2", "radTmed1",
                "radTmed2", "radTdist1", "radTdist2", "lm_thick1",
                "lm_thick2", "lm_medium1",
                "lm_medium2", "lm_thin1", "lm_thin2", "rad_t1",
                "rad_t2", "rad_t3"]
    t = h.Vector().record(h._ref_t, dt)
    # calcium_indicator = {}
    calcium = {}
    headers = {}
    for sec_name in sec_list:
        #calcium_indicator[sec_name] = []
        calcium[sec_name] = []
        headers[sec_name] = []
        for shell in cell.shells[sec_name]:
            #ind, header = uf.record_specie_vec(cell.indicator_ca, shell, 1)
            #calcium_indicator[sec_name].append(ind)
            ca, header = uf.record_specie_vec(cell.ca, shell, 1)
            calcium[sec_name].append(ca)
            headers[sec_name].append(header)

    t_vec = h.Vector()
    t_vec.record(h._ref_t, 1)
    start = time.time()
    cell.make_a_run(t_stop)
    print(time.time() - start)
    dends = {}
    # for sec_name in sec_list:
    #     fig, ax = plt.subplots(1, 1)
    #     dend = cell.find_sec(sec_name)
    #     indicator = []
    #     for c in calcium_indicator[sec_name]:
    #         shell_data = []
    #         for seg in c:
    #             shell_data.append(seg.as_numpy())
    #         indicator.append(np.array(shell_data))
    #     out = uf.dend_average(indicator, dend, cell.factors[sec_name])
    #     dends[sec_name] = dend
    #     header = headers[sec_name][0].split()
    #     for i, signal in enumerate(out):
    #         ax.plot(t_vec.as_numpy(), signal, label=header[i].split("_")[-1])
    #     ax.set_title("Calcium indicator in %s" % sec_name)
    #     ax.set_xlabel("time [ms]")

    for sec_name in sec_list:
        fig, ax = plt.subplots(1, 1)
        dend = cell.find_sec(sec_name)
        indicator = []
        for c in calcium[sec_name]:
            shell_data = []
            for seg in c:
                shell_data.append(seg.as_numpy())
            indicator.append(np.array(shell_data))
        out = uf.dend_average(indicator, dend, cell.factors[sec_name])
        dends[sec_name] = dend
        header = headers[sec_name][0].split()
        for i, signal in enumerate(out):
            ax.plot(t_vec.as_numpy(), signal, label=header[i].split("_")[-1])
        ax.set_title("Calcium in %s" % sec_name)
        ax.set_xlabel("time [ms]")
    plt.show()
