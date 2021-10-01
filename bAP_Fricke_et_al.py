# A modified script from Tomko et al 2021"""produces a simulation of the backpropagation of the action potential into the apical trunk"""# _title_   : bAP.py# _author_  : Matus Tomko# _mail_    : matus.tomko __at__ fmph.uniba.skimport osfrom collections import OrderedDictfrom quantities import mVimport efelimport matplotlib.pyplot as pltimport numpy as npfrom neuron import hfrom CA1 import CA1_PCimport utility_functions as ufsavepath = './figs/'if not os.path.exists(savepath):    os.mkdir(savepath)def extract_AP1_amps(v_vec_soma, v_vecs_apical_trunk, t_vec):    traces = []    trace = {}    trace['T'] = np.array(t_vec)    trace['V'] = np.array(v_vec_soma)    trace['stim_start'] = [500]    trace['stim_end'] = [1500]    traces.append(trace)    efel.reset()    efel.setDoubleSetting('interp_step', 0.025)    efel.setDoubleSetting('DerivativeThreshold', 40.0)    efel_results = efel.getFeatureValues(traces, ['inv_first_ISI', 'AP_begin_time', 'doublet_ISI'])    efel.reset()    soma_AP_begin_time = efel_results[0]['AP_begin_time']    soma_first_ISI = efel_results[0]['doublet_ISI'][0]    s_indices_AP1 = np.where(np.array(t_vec) >= (soma_AP_begin_time[0] - 1.0))    if 10 < soma_first_ISI:        plus = 10    else:        plus = soma_first_ISI - 3    e_indices_AP1 = np.where(np.array(t_vec) >= (soma_AP_begin_time[0] + plus))    start_index_AP1 = s_indices_AP1[0][0]    end_index_AP1 = e_indices_AP1[0][0]    AP1_amps = OrderedDict()    for vec in v_vecs_apical_trunk:        trace = np.array(v_vecs_apical_trunk[vec])        amp1 = float(            np.amax(np.array(trace)[start_index_AP1:end_index_AP1]) - np.array(trace)[start_index_AP1]) * mV        AP1_amps[vec] = amp1    return AP1_ampsif __name__ == "__main__":    where_spines = []    dt = 1    add_ER = False    where_ca = ["soma", "apical"]    t_stop = 26000    t_basal = 23000    t_stim = 25000    buffer_list = ["Calmodulin", "Calbindin", "Fixed", "OGB1"]    cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines,                  buffer_list=buffer_list)        stim = h.IClamp(cell.soma(0.5))    stim.delay = t_stim    stim.amp = 2    stim.dur = 4    sec_for_Vm = {'50': cell.find_sec("radTprox2")(0.166667),                  '150': cell.find_sec("radTmed2")(0.166667),                  '245.45': cell.find_sec("radTdist1")(0.5),                  '263.63': cell.find_sec("radTdist1")(0.7),                  '336.36': cell.find_sec("radTdist2")(0.3),                  '354.54': cell.find_sec("radTdist2")(0.5)}    sec_Ca_trunk = {"10": cell.find_sec("radTprox1")(0.1),                    "20": cell.find_sec("radTprox1")(0.3),                    "30": cell.find_sec("radTprox1")(0.5),                    "40": cell.find_sec("radTprox1")(0.7),                    "50": cell.find_sec("radTprox1")(0.9),                    "60": cell.find_sec("radTprox2")(0.1),                    "70": cell.find_sec("radTprox2")(0.3),                    "100": cell.find_sec("radTprox2")(0.9),                    "110": cell.find_sec("radTmed1")(0.1),                    "120": cell.find_sec("radTmed1")(0.3),                    "130": cell.find_sec("radTmed1")(0.5),                    "140": cell.find_sec("radTmed1")(0.7),                    "150": cell.find_sec("radTmed1")(0.9),                    "160": cell.find_sec("radTmed2")(0.1),                    "170": cell.find_sec("radTmed2")(0.3),    }        v_vec_soma = h.Vector().record(cell.soma(0.5)._ref_v, dt)    v_vecs_apical_trunk = OrderedDict()    v_vecs_apical_trunk['50'] = h.Vector().record(sec_for_Vm['50']._ref_v, dt)    v_vecs_apical_trunk['150'] = h.Vector().record(sec_for_Vm['150']._ref_v, dt)    v_vecs_apical_trunk['245.45'] = h.Vector().record(sec_for_Vm['245.45']._ref_v, dt)    v_vecs_apical_trunk['263.63'] = h.Vector().record(sec_for_Vm['263.63']._ref_v, dt)    v_vecs_apical_trunk['336.36'] = h.Vector().record(sec_for_Vm['336.36']._ref_v, dt)    v_vecs_apical_trunk['354.54'] = h.Vector().record(sec_for_Vm['354.54']._ref_v, dt)    sec_list = ["soma", "radTprox1", "radTprox2", "radTmed1",                "radTmed2", "radTdist1", "radTdist2", "lm_thick1",                "lm_thick2", "lm_medium1",                "lm_medium2", "lm_thin1", "lm_thin2", "rad_t1",                "rad_t2", "rad_t3"]    soma = cell.soma    ica_soma = h.Vector().record(soma(0.5)._ref_ica, 1)        calcium_indicator = {}    headers = {}    for sec_name in sec_list:        calcium_indicator[sec_name] = []        headers[sec_name] = []        for shell in cell.shells[sec_name]:            ind, header = uf.record_specie_vec(cell.indicator_ca, shell, dt)            headers[sec_name].append(header)            calcium_indicator[sec_name].append(ind)        t_vec = h.Vector()    t_vec.record(h._ref_t, 1)        cell.make_a_run(t_stop)    fig = plt.figure(figsize=(5, 5))    plt.plot(t_vec, v_vec_soma, label='soma')    for vec in v_vecs_apical_trunk:        plt.plot(t_vec, v_vecs_apical_trunk[vec], label=vec + ' um')    plt.xlabel('Time (ms)')    plt.ylabel('Membrane potential (mV)')    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')    plt.savefig(savepath + 'traces_Fricke.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')    fig = plt.figure(figsize=(5, 5))    plt.plot(t_vec, ica_soma, label='soma')    plt.xlabel('Time (ms)')    plt.ylabel('Ca current (nA)')    plt.title('First AP')    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')    plt.savefig(savepath + 'ica_traces_Fricke.png', format='png', bbox_extra_artists=(lgd,),                bbox_inches='tight')    dends = {}    peaks = {}    for sec_name in sec_list:        fig, ax = plt.subplots(1, 1)        dend = cell.find_sec(sec_name)        indicator = []        for c in calcium_indicator[sec_name]:            shell_data = []            for seg in c:                shell_data.append(seg.as_numpy())            indicator.append(np.array(shell_data))        out = uf.dend_average(indicator, dend, cell.factors[sec_name])        dends[sec_name] = dend        header = headers[sec_name][0].split()        for i, signal in enumerate(out):            signal_basal = signal[int(t_basal/dt):int(t_stim/dt)].mean()            new_signal = (signal[int((t_stim-300)/dt):int((t_stim+300)/dt)]                          -signal_basal)/signal_basal            new_t = t_vec.as_numpy()[int((t_stim-300)/dt):int((t_stim+300)/dt)]            ax.plot(new_t, new_signal, label=header[i].split("_")[-1])        ax.set_title(sec_name)        ax.set_xlabel("time [ms]")    plt.show()        # amps = extract_AP1_amps(v_vec_soma, v_vecs_apical_trunk, t_vec)    # fig = plt.figure(figsize=(5, 5))    # for d in amps:    #     plt.plot(d, amps[d], marker='o', linestyle='none', label=d)    # plt.xlabel('Distance from the soma (um)')    # plt.ylabel('AP1_amp (mV)')    # lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')    # plt.savefig(savepath + 'AP1_amps.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')   