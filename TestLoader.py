from builtins import range
import os
import numpy
import sciunit
import hippounit.capabilities as cap
from quantities import ms,mV,Hz
from neuron import h
from subprocess import run
import multiprocessing
import zipfile
import collections

import collections

import json

import pkg_resources
import sys



class ModelLoader(sciunit.Model,
                 cap.ProvidesGoodObliques,
                 cap.ReceivesSquareCurrent_ProvidesResponse,
                 cap.ReceivesSynapse,
                 cap.ReceivesMultipleSynapses,
                 cap.ReceivesSquareCurrent_ProvidesResponse_MultipleLocations,
                 cap.ProvidesRecordingLocationsOnTrunk,
                 cap.ProvidesRandomDendriticLocations,
                 cap.ReceivesEPSCstim):

    def __init__(self, model_class, mods_dir, model_params):
        """ Constructor. """

        """ This class should be used with Jupyter notebooks"""
        self.class_name = model_class
        self.model_args = model_params
        self.modelpath = mods_dir
        self.model_args["recompile"] = False
        self.name = "CA1"
        self.max_dist_from_soma = 150
        self.v_init = -70
        self.celsius = 34
        self.c_step_start = 0.00004
        self.c_step_stop = 0.000004
        self.c_minmax = numpy.array([0.00004, 0.04])
        self.threshold = -20
        self.stim = None
        self.soma = None
        self.NMDA_name = "NR2A_Ca"
        self.AMPA_name = "AMPADJ"
        sciunit.Model.__init__(self, name=self.name)
        self.dend_loc = []  
        self.dend_locations = collections.OrderedDict()
        self.base_directory = './validation_results/'   
        self.compile_mod_files()

    def compile_mod_files(self):
        if self.modelpath is None:
            raise Exception("""Please give the path to the mod files (eg. mod_files_path = \'/home/models/CA1_pyr/mechanisms/\') 
            as an argument to the ModelLoader class""")

        #if os.path.isfile(self.modelpath + self.libpath) is False:
        working_dir = os.getcwd()
        os.chdir(self.modelpath)
        p = run('nrnivmodl')
        os.chdir(working_dir)

    def translate(self, sectiontype, distance=0):
        if "soma" in sectiontype:
            return "soma"
        else:
            return False


    def initialize(self, args):
        save_stdout = sys.stdout
        sys.stdout = open('/dev/stdout', 'w')     
        h.load_file("stdrun.hoc")
        cell = self.class_name(**args)
        try:
            self.soma = cell.soma[0]
        except TypeError:
            self.soma = cell.soma
        self.cell = cell
        sys.stdout = save_stdout    #setting output back to normal
        return cell

    def inject_current(self, amp, delay, dur, section_stim,
                       loc_stim, section_rec, loc_rec):

        self.initialize(self.model_args)
        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        stim_s_name = self.translate(section_stim, distance=0)
        rec_sec_name = self.translate(section_rec, distance=0)
        new_sec = self.cell.find_sec(stim_s_name)
        self.sect_loc_stim = new_sec(float(loc_stim))
        print("- running amplitude: %f on model: %s at: %s(%s)" % (amp,
                                                                   self.name,
                                                                   stim_s_name,
                                                                   loc_stim))

        self.stim = h.IClamp(self.sect_loc_stim)
        self.stim.amp = amp
        self.stim.delay = delay
        self.stim.dur = dur
        new_sec = self.cell.find_sec(rec_sec_name)
        self.sect_loc_rec = new_sec(float(loc_rec))
        rec_t = h.Vector()
        rec_t.record(h._ref_t)
        rec_v = h.Vector()
        rec_v.record(self.sect_loc_rec._ref_v)
        h.stdinit()
        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init
        h.celsius = self.celsius
        h.init()
        h.tstop = delay + dur + 200
        h.run()
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        return t, v

    def inject_current_record_respons_multiple_loc(self, amp, delay,
                                                   dur, section_stim,
                                                   loc_stim,
                                                   dend_locations):
        self.initialize(self.model_args)
        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        stim_s_name = self.translate(section_stim, distance=0)
        new_sec = self.cell.find_sec(stim_s_name)
        self.sect_loc_stim = new_sec(float(loc_stim))
        self.sect_loc_rec = new_sec(float(loc_stim))
        print("- running amplitude: %f on model: %s at: %s(%s)" % (amp,
                                                                   self.name,
                                                                   stim_s_name,
                                                                   loc_stim))

        self.stim = h.IClamp(self.sect_loc_stim)

        self.stim.amp = amp
        self.stim.delay = delay
        self.stim.dur = dur

        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v_stim = h.Vector()
        rec_v_stim.record(self.sect_loc_rec._ref_v)

        rec_v = []
        v = collections.OrderedDict()
        self.dend_loc_rec =[]

        #print dend_locations
        for key, value in dend_locations.items():
            print(key, value)
            for x in value:
                print(x)
                new_sec = self.cell.find_sec(x[0])
                self.dend_loc_rec.append(new_sec(x[1]))
                rec_v.append(h.Vector())

        for i, sec in enumerate(self.dend_loc_rec):
            rec_v[i].record(sec._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init#-65

        h.celsius = self.celsius
        h.init()
        h.tstop = delay + dur + 200
        h.run()

        t = numpy.array(rec_t)
        v_stim = numpy.array(rec_v_stim)

        '''
        for i in range(0,len(dend_loc)):
            v.append(numpy.array(rec_v[i]))
        '''

        i = 0
        for key, value in dend_locations.items():
            v[key] = collections.OrderedDict()
            for j in range(len(dend_locations[key])):
                loc_key = (dend_locations[key][j][0], dend_locations[key][j][1])
                # list can not be a key, but tuple can
                v[key][loc_key] = numpy.array(rec_v[i])
                # the list that specifies dendritic location will be a key too.
                i += 1
        return t, v_stim, v


    def find_trunk_locations(self, distances, tolerance, trunk_origin):
        self.initialize(self.model_args)
        locations = collections.OrderedDict()
        actual_distances = {}

        for sec in self.cell.trunk:
            #for seg in sec:
            if not trunk_origin:
                h.distance(sec=self.soma)
            elif len(trunk_origin) == 1:
                h.distance(sec=trunk_origin[0])
            else:
                h.distance(sec=trunk_origin)
            for seg in sec:
                for i in range(0, len(distances)):
                    # if this key doesn't exist it is added with the value: [],
                    # if it exists, value not altered
                    locations.setdefault(distances[i], [])
                    # if the seq is between distance +- 20
                    dist = h.distance(seg.x, sec=sec)
                    if (dist < (distances[i] + tolerance)
                        and dist > (distances[i]- tolerance)): 
                        locations[distances[i]].append([sec.name(),
                                                         seg.x])
                        actual_distances[sec.name(), seg.x] = dist
        return locations, actual_distances

    def get_random_locations(self, num, seed, dist_range, trunk_origin):

        locations=[]
        locations_distances = {}

        self.initialize(self.model_args)
        kumm_length_list = []
        kumm_length = 0
        num_of_secs = 0

        for sec in self.cell.trunk:
            #print sec.L
            num_of_secs += sec.nseg
            kumm_length += sec.L
            kumm_length_list.append(kumm_length)
        #print 'kumm' ,kumm_length_list
        #print num_of_secs

        if num > num_of_secs:
            for sec in self.cell.trunk:
                if not trunk_origin:
                    h.distance(sec=self.soma)
                elif len(trunk_origin) == 1:
                    h.distance(sec=trunk_origin[0])
                else:
                    h.distance(sec=trunk_origin)
                for seg in sec:
                    dist = h.distance(seg.x, sec=sec)
                    if dist > dist_range[0] and dist < dist_range[1]:  
                        locations.append([sec.name(), seg.x])
                        locations_distances[sec.name(), seg.x] = dist
        else:
            norm_kumm_length_list = [i/kumm_length_list[-1] for i in kumm_length_list]
            import random

            _num_ = num  # _num_ will be changed
            num_iterations = 0
            random.seed(seed)

            while len(locations) < num and num_iterations < 50 :
                #print 'seed ', seed
                rand_list = [random.random() for j in range(_num_)]
                #print rand_list

                for rand in rand_list:
                    #print 'RAND', rand
                    for i in range(len(norm_kumm_length_list)):
                        if (rand <= norm_kumm_length_list[i]
                            and (rand > norm_kumm_length_list[i-1]
                                 or i==0)):
                            #print norm_kumm_length_list[i-1]
                            #print norm_kumm_length_list[i]
                            seg_loc = ((rand - norm_kumm_length_list[i-1]) /
                                       (norm_kumm_length_list[i] -
                                        norm_kumm_length_list[i-1]))
                            #print 'seg_loc', seg_loc
                            segs = [seg.x for seg in self.cell.trunk[i]]
                            d_seg = [abs(seg.x - seg_loc) for seg in
                                     self.cell.trunk[i]]
                            min_d_seg = numpy.argmin(d_seg)
                            segment = segs[min_d_seg]
                            #print 'segment', segment
                            if not trunk_origin:
                                h.distance(sec=self.soma)
                            elif len(trunk_origin) == 1:
                                h.distance(sec=trunk_origin[0]) 
                            else:
                                h.distance(sec=trunk_origin)
                            dist = h.distance(segment, sec=self.cell.trunk[i])
                            if ([self.cell.trunk[i].name(), segment] not in locations
                                and dist >= dist_range[0]
                                and dist < dist_range[1]):
                                locations.append([self.cell.trunk[i].name(),
                                                  segment])
                                locations_distances[self.cell.trunk[i].name(),
                                                    segment] = dist
                _num_ = num - len(locations)

                seed += 10
                num_iterations += 1

        return locations, locations_distances

    def find_good_obliques(self, trunk_origin):
        """Used in ObliqueIntegrationTest"""
        self.initialize(self.model_args)
        dend_loc = []
        for sec in self.cell.oblique:
            dend_loc_prox=[]
            dend_loc_dist=[]
            seg_list_prox=[]
            seg_list_dist=[]
 
            h.distance(sec=sec)
            #set the 0 point of the section as the origin
            print(sec.name())


            for seg in sec:
                # print(seg.x, h.distance(seg.x))
                if h.distance(seg.x, sec=sec) > 5 and h.distance(seg.x, sec=sec) < 50:
                    seg_list_prox.append(seg.x)
                if h.distance(seg.x, sec=sec) > 60 and h.distance(seg.x, sec=sec) < 126:
                    seg_list_dist.append(seg.x)

            #print seg_list_prox
            #print seg_list_dist

            if len(seg_list_prox) > 1:
                s = int(numpy.ceil(len(seg_list_prox)/2.0))
                dend_loc_prox.append(sec.name())
                dend_loc_prox.append(seg_list_prox[s])
                dend_loc_prox.append('prox')
            elif len(seg_list_prox) == 1:
                dend_loc_prox.append(sec.name())
                dend_loc_prox.append(seg_list_prox[0])
                dend_loc_prox.append('prox')

            if len(seg_list_dist) > 1:
                s = int(numpy.ceil(len(seg_list_dist)/2.0)-1)
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(seg_list_dist[s])
                dend_loc_dist.append('dist')
            elif len(seg_list_dist) == 1:
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(seg_list_dist[0])
                dend_loc_dist.append('dist')
            elif len(seg_list_dist) == 0:                # if the dendrite is not long enough to meet the criteria, we stimulate its end
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(0.9)
                dend_loc_dist.append('dist')

            if dend_loc_prox:
                dend_loc.append(dend_loc_prox)
            if dend_loc_dist:
                dend_loc.append(dend_loc_dist)

        #print 'Dendrites and locations to be tested: ', dend_loc

        return dend_loc

    def set_netstim_netcon(self, interval, number):
        """Used in ObliqueIntegrationTest"""
        print(interval, number)
        self.presynaptic = []
        self.release = []
        self.nc_list = []
        self.ns_list = []
        print(len(self.cell.ampas))
        for i in range(number):
            self.presynaptic.append(h.Section("PRE_%d" % i))
            self.release.append(h.depletion(self.presynaptic[i](0.5)))

        for i in range(number):
            self.ns_list.append(h.NetStim())
            self.ns_list[i].number = 1
            self.ns_list[i].start = self.start + (i*interval)
            self.nc_list.append(h.NetCon(self.ns_list[i], self.release[i], 0, 0, 1))
            h.setpointer(self.release[i]._ref_T, 'T', self.cell.ampas[i])
            if len(self.cell.nmdas):
                h.setpointer(self.release[i]._ref_T, 'T', self.cell.nmdas[i]) 

    def run_syn(self, dend_loc, interval, number, AMPA_weight):
        """Currently not used - Used to be used in ObliqueIntegrationTest"""
        args = self.model_args
        args["spine_pos"] = {}
        args["spine_pos"][dend_loc[0]] = [dend_loc[1]]
        args["where_spines"] = [dend_loc[0]]
        self.initialise(args)
        self.dendrite = self.cell.find_sec(dend_loc[0])
        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        self.set_netstim_netcon(interval, 1)
        self.set_num_weight(0, 1, 1)

        self.sect_loc=self.soma(0.5)

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/ dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop = 500
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend


    def run_multiple_syn(self, dend_loc, interval, number, weight):
        """Used in ObliqueIntegrationTest"""
        self.start = 300
        args = self.model_args
        args["spine_pos"] = {}
        args["spine_pos"][dend_loc[0]] = []
        dx = 1/150
        for i in range(number):
            args["spine_pos"][dend_loc[0]].append(dend_loc[1]+i*dx)
        args["where_spines"] = [dend_loc[0]]
        self.initialize(args)
        
        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)
        self.dendrite = self.cell.find_sec(dend_loc[0])
        self.xloc = dend_loc[1]

        self.set_netstim_netcon(interval, number)
        self.sect_loc = self.soma(0.5)

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop =500
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend

    def run_EPSCstim(self, dend_loc, weight, tau1, tau2):
        """Used in PSPAttenuationTest"""
        args = self.model_args
        args["spine_pos"] = {}
        args["spine_pos"][dend_loc[0]] = [dend_loc[1]]
        args["where_spines"] = [dend_loc[0]]
        args["receptor_list"] = ["AMPA"]
        self.initialize(args)
        self.start = 300
        
        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        self.set_netstim_netcon(0, 1)
 
        self.sect_loc = self.soma(0.5)
        self.dendrite = self.cell.find_sec(dend_loc[0])
        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop = 450
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend

