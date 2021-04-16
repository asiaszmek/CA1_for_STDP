"""produces a simulation of the backpropagation of the action potential into the apical trunk"""

# _title_   : bAP.py
# _author_  : Matus Tomko
# _mail_    : matus.tomko __at__ fmph.uniba.sk

import os
from collections import OrderedDict
from quantities import mV

import efel
import matplotlib.pyplot as plt
import numpy
from neuron import h, gui

savepath = './figs/'
if not os.path.exists(savepath):
    os.mkdir(savepath)


def main():
    h.nrn_load_dll('./Mods/nrnmech.dll')
    h.xopen('pyramidal_cell_weak_bAP_original.hoc')
    cell = h.CA1_PC_Tomko()

    stim = h.IClamp(cell.soma[0](0.5))
    stim.delay = 500
    stim.amp = 0.8
    stim.dur = 1000

    v_vec_soma = h.Vector().record(cell.soma[0](0.5)._ref_v)
    v_vecs_apical_trunk = OrderedDict()
    v_vecs_apical_trunk['50'] = h.Vector().record(cell.radTprox(0.5)._ref_v)
    v_vecs_apical_trunk['150'] = h.Vector().record(cell.radTmed(0.5)._ref_v)
    v_vecs_apical_trunk['245.45'] = h.Vector().record(cell.radTdist(0.22727272727272727)._ref_v)
    v_vecs_apical_trunk['263.63'] = h.Vector().record(cell.radTdist(0.3181818181818182)._ref_v)
    v_vecs_apical_trunk['336.36'] = h.Vector().record(cell.radTdist(0.6818181818181819)._ref_v)
    v_vecs_apical_trunk['354.54'] = h.Vector().record(cell.radTdist(0.7727272727272728)._ref_v)

    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    h.dt = 0.025
    h.tstop = 1700
    h.v_init = -65
    h.celsius = 35
    h.init()
    h.finitialize(-65)
    h.cvode_active(1)
    h.run()


    fig = plt.figure(figsize=(5, 5))
    plt.plot(t_vec, v_vec_soma, label='soma')
    for vec in v_vecs_apical_trunk:
        plt.plot(t_vec, v_vecs_apical_trunk[vec], label=vec + ' um')
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane potential (mV)')
    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    plt.savefig(savepath + 'traces.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(5, 5))
    plt.plot(t_vec, v_vec_soma, label='soma')
    for vec in v_vecs_apical_trunk:
        plt.plot(t_vec, v_vecs_apical_trunk[vec], label=vec + ' um')
    plt.xlim((503, 513))
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane potential (mV)')
    plt.title('First AP')
    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    plt.savefig(savepath + 'AP1_traces.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    amps = extract_AP1_amps(v_vec_soma, v_vecs_apical_trunk, t_vec)
    fig = plt.figure(figsize=(5, 5))
    for d in amps:
        plt.plot(d, amps[d], marker='o', linestyle='none', label=d)
    plt.xlabel('Distance from the soma (um)')
    plt.ylabel('AP1_amp (mV)')
    plt.ylim(0, 70)
    lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    plt.savefig(savepath + 'AP1_amps.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def extract_AP1_amps(v_vec_soma, v_vecs_apical_trunk, t_vec):
    traces = []
    trace = {}
    trace['T'] = numpy.array(t_vec)
    trace['V'] = numpy.array(v_vec_soma)
    trace['stim_start'] = [500]
    trace['stim_end'] = [1500]
    traces.append(trace)

    efel.reset()
    efel.setDoubleSetting('interp_step', 0.025)
    efel.setDoubleSetting('DerivativeThreshold', 40.0)
    efel_results = efel.getFeatureValues(traces, ['inv_first_ISI', 'AP_begin_time', 'doublet_ISI'])
    efel.reset()

    soma_AP_begin_time = efel_results[0]['AP_begin_time']
    soma_first_ISI = efel_results[0]['doublet_ISI'][0]
    s_indices_AP1 = numpy.where(numpy.array(t_vec) >= (soma_AP_begin_time[0] - 1.0))
    if 10 < soma_first_ISI:
        plus = 10
    else:
        plus = soma_first_ISI - 3
    e_indices_AP1 = numpy.where(numpy.array(t_vec) >= (soma_AP_begin_time[0] + plus))
    start_index_AP1 = s_indices_AP1[0][0]
    end_index_AP1 = e_indices_AP1[0][0]

    AP1_amps = OrderedDict()
    for vec in v_vecs_apical_trunk:
        trace = numpy.array(v_vecs_apical_trunk[vec])
        amp1 = float(
            numpy.amax(numpy.array(trace)[start_index_AP1:end_index_AP1]) - numpy.array(trace)[start_index_AP1]) * mV
        AP1_amps[vec] = amp1

    return AP1_amps

if __name__ == '__main__':
    main()
