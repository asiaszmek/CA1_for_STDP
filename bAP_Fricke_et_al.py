# A modified script from Tomko et al 2021import osfrom collections import OrderedDictfrom quantities import mVimport efelimport matplotlib.pyplot as pltimport numpy as npfrom neuron import hfrom CA1 import CA1_PCimport utility_functions as ufsavepath = './figs/'if not os.path.exists(savepath):    os.mkdir(savepath)if __name__ == "__main__":    where_spines = []    dt = 0.25    add_ER = False    where_ca = ["soma", "apical"]    t_stop = 6000    t_basal = 5000    t_stim = 5300    buffer_list = ["Calmodulin", "Calbindin", "Fixed", "OGB1"]    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines,                  buffer_list=buffer_list)        stim = h.IClamp(cell.soma(0.5))    stim.delay = t_stim    stim.amp = 2    stim.dur = 1    h.distance(sec=cell.soma)    sec_for_Vm = {'50': cell.find_sec("radTprox")(0.409091),                  '150': cell.find_sec("radTmed")(0.409091),                  '245.45': cell.find_sec("radTdist")(0.166667),                  '263.63': cell.find_sec("radTdist")(0.261905),                  '336.36': cell.find_sec("radTdist")(0.642857),                  '354.54': cell.find_sec("radTdist")(0.738095)}    sec_Ca_trunk = {"10": cell.find_sec("radTprox")(0.0454545),                    "20": cell.find_sec("radTprox")(0.136364),                    "30": cell.find_sec("radTprox")(0.227273),                    "40": cell.find_sec("radTprox")(0.318182),                    "50": cell.find_sec("radTprox")(0.409091),                    "60": cell.find_sec("radTprox")(0.5),                    "70": cell.find_sec("radTprox")(0.590909),                    "100": cell.find_sec("radTprox")(0.954545),                    "110": cell.find_sec("radTmed")(0.0454545),                    "120": cell.find_sec("radTmed")(0.136364),                    "130": cell.find_sec("radTmed")(0.227273),                    "140": cell.find_sec("radTmed")(0.318182),                    "150": cell.find_sec("radTmed")(0.409091),                    "160": cell.find_sec("radTmed")(0.5),                    "170": cell.find_sec("radTmed")(0.590909),    }        v_vec_soma = h.Vector().record(cell.soma(0.5)._ref_v, dt)    v_vecs_apical_trunk = OrderedDict()    v_vecs_apical_trunk['50'] = h.Vector().record(sec_for_Vm['50']._ref_v, dt)    v_vecs_apical_trunk['150'] = h.Vector().record(sec_for_Vm['150']._ref_v, dt)    v_vecs_apical_trunk['245.45'] = h.Vector().record(sec_for_Vm['245.45']._ref_v, dt)    v_vecs_apical_trunk['263.63'] = h.Vector().record(sec_for_Vm['263.63']._ref_v, dt)    v_vecs_apical_trunk['336.36'] = h.Vector().record(sec_for_Vm['336.36']._ref_v, dt)    v_vecs_apical_trunk['354.54'] = h.Vector().record(sec_for_Vm['354.54']._ref_v, dt)    sec_list = ["soma", "radTprox",                "radTmed", "radTdist", "lm_thick1",                "lm_thick2", "lm_medium1",                "lm_medium2", "lm_thin1", "lm_thin2", "rad_t1",                "rad_t2", "rad_t3"]    soma = cell.soma    ica_soma = h.Vector().record(soma(0.5)._ref_ica, dt)        calcium_indicator = {}    headers = {}    for sec_name in sec_list:        calcium_indicator[sec_name] = []        headers[sec_name] = []        for shell in cell.shells[sec_name]:            ind, header = uf.record_specie_vec(cell.indicator_ca, shell, dt)            headers[sec_name].append(header)            calcium_indicator[sec_name].append(ind)        t_vec = h.Vector()    t_vec.record(h._ref_t, dt)        cell.make_a_run(t_stop)    fig = plt.figure(figsize=(5, 5))    plt.plot(t_vec, v_vec_soma, label='soma')    for vec in v_vecs_apical_trunk:        plt.plot(t_vec, v_vecs_apical_trunk[vec], label=vec + ' um')    plt.xlabel('Time (ms)')    plt.ylabel('Membrane potential (mV)')    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')    plt.savefig(savepath + 'traces_Fricke.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')    fig = plt.figure(figsize=(5, 5))    plt.plot(t_vec, ica_soma, label='soma')    plt.xlabel('Time (ms)')    plt.ylabel('Ca current (nA)')    plt.title('First AP')    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')    plt.savefig(savepath + 'ica_traces_Fricke.png', format='png', bbox_extra_artists=(lgd,),                bbox_inches='tight')    dends = {}    peaks = {}    for sec_name in sec_list:        fig, ax = plt.subplots(1, 1)        dend = cell.find_sec(sec_name)        indicator = []        for c in calcium_indicator[sec_name]:            shell_data = []            for seg in c:                shell_data.append(seg.as_numpy())            indicator.append(np.array(shell_data))        out = uf.dend_average(indicator, dend, cell.factors[sec_name])        dends[sec_name] = dend        header = headers[sec_name][0].split()        for i, signal in enumerate(out):            signal_basal = signal[int(t_basal/dt):int(t_stim/dt)].mean()            new_signal = (signal[int((t_stim-300)/dt):int((t_stim+300)/dt)]                          -signal_basal)/signal_basal            new_t = t_vec.as_numpy()[int((t_stim-300)/dt):int((t_stim+300)/dt)]            ax.plot(new_t, new_signal, label=header[i].split("_")[-1])        ax.set_title(sec_name)        ax.set_xlabel("time [ms]")    plt.show()    