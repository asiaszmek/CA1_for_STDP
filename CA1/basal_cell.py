import os
from subprocess import run

from numpy import exp
import neuron


def e_pas_dist(x):
    return -65.73-5*x/150


def ghdbar_dist(x):
    return (1. + 3./100. * x)*1.90e-05


def gkabar_dist(x):
    return (15./(1. + exp((300-x)/50)))* 0.013


my_loc = os.path.dirname(os.path.abspath(__file__))
mechanisms_path = os.path.join(my_loc, "Mods")

      
class CA1_PC_basal:
    sections = []
    def __init__(self, recompile=True, add_to_h=True):
        self.recompile = recompile
        self.add_to_h = add_to_h
        self.load_neuron()

    def make_sections(self):
        self.soma = neuron.h.Section(name="soma")
        self.radTprox1 = neuron.h.Section(name="radTprox1")
        self.radTprox2 = neuron.h.Section(name="radTprox2")
        self.radTmed1 = neuron.h.Section(name="radTmed1")
        self.radTmed2 = neuron.h.Section(name="radTmed2")
        self.radTdist1 = neuron.h.Section(name="radTdist1")
        self.radTdist2 = neuron.h.Section(name="radTdist2")
        self.lm_thick1 = neuron.h.Section(name="lm_thick1")
        self.lm_thick2 = neuron.h.Section(name="lm_thick2")
        self.lm_medium1 = neuron.h.Section(name="lm_medium1")
        self.lm_medium2 = neuron.h.Section(name="lm_medium2")
        self.lm_thin1 = neuron.h.Section(name="lm_thin1")
        self.lm_thin2 = neuron.h.Section(name="lm_thin2")
        self.axon = neuron.h.Section(name="axon")
        self.oriprox1 = neuron.h.Section(name="oriprox1")
        self.oriprox2 = neuron.h.Section(name="oriprox2")
        self.oridist1_1 = neuron.h.Section(name="oridist1_1")
        self.oridist1_2 = neuron.h.Section(name="oridist1_2")
        self.oridist2_1 = neuron.h.Section(name="oridist2_1")
        self.oridist2_2 = neuron.h.Section(name="oridist2_2")
        self.rad_t1 = neuron.h.Section(name="rad_t1")
        self.rad_t2 = neuron.h.Section(name="rad_t2")
        self.rad_t3 = neuron.h.Section(name="rad_t3")
        
    def connections(self):
        self.radTprox1.connect(self.soma(1))
        self.radTprox2.connect(self.radTprox1(1))
        self.radTmed1.connect(self.radTprox2(1))
        self.radTmed2.connect(self.radTmed1(1))
        self.radTdist1.connect(self.radTmed2(1))
        self.radTdist2.connect(self.radTdist1(1))
        self.rad_t1.connect(self.radTprox1(0.5))
        self.rad_t2.connect(self.radTmed1(1))
        self.rad_t3.connect(self.radTdist1(1))
        
        self.lm_thick2.connect(self.radTdist2(1))
        self.lm_medium2.connect(self.lm_thick2(1))
        self.lm_thin2.connect(self.lm_medium2(1))
        self.lm_thick1.connect(self.radTdist2(1))
        self.lm_medium1.connect(self.lm_thick1(1))
        self.lm_thin1.connect(self.lm_medium1(1))

        self.oriprox1.connect(self.soma(0))
        self.oridist1_1.connect(self.oriprox1(1))
        self.oridist1_2.connect(self.oriprox1(1))
        self.oriprox2.connect(self.soma(1))
        self.oridist2_1.connect(self.oriprox2(1))
        self.oridist2_2.connect(self.oriprox2(1))
        self.axon.connect(self.soma(1))

    def make_lists(self):
        self.sections = [self.soma, self.axon, self.radTprox1,
                         self.radTprox2, self.radTmed1, self.radTmed2,
                         self.radTdist1,
                         self.radTdist2, self.lm_thick2, self.lm_medium2,
                         self.lm_thin2,
                         self.lm_thick1,self.lm_medium1, self.lm_thin1,
                         self.oriprox1, self.oridist1_1, self.oridist1_2,
                         self.oriprox2, self.oridist2_1, self.oridist2_2,
                         self.rad_t1, self.rad_t2, self.rad_t3]

        self.somatic = [self.soma]
        self.axonal = [self.axon]
        self.apical = [self.radTprox1, self.radTprox2, self.radTmed1,
                       self.radTmed2, self.radTdist1,
                       self.radTdist2, self.rad_t1, self.rad_t2, self.rad_t3,
                       self.lm_thick2, self.lm_medium2,
                       self.lm_thin2, self.lm_thick1, self.lm_medium1,
                       self.lm_thin1]
        self.basal = [self.oriprox1, self.oridist1_1, self.oridist1_2,
                      self.oriprox2, self.oridist2_1, self.oridist2_2]
        self.trunk = [
            self.radTprox1,
            self.radTprox2,
            self.radTmed1,
            self.radTmed2,
            self.radTdist1,
            self.radTdist2,
        ]
        self.oblique = [self.rad_t2]

    def make_basic_shape(self):
        for sec in self.sections:
            neuron.h.pt3dclear(sec=sec)
        neuron.h.pt3dadd(0, 0, 0, 1, sec=self.soma)
        neuron.h.pt3dadd(15, 0, 0, 1, sec=self.soma)
        neuron.h.pt3dadd(15, 0, 0, 1, sec=self.radTprox1)
        neuron.h.pt3dadd(15, 15, 0, 1, sec=self.radTprox1)
        neuron.h.pt3dadd(15, 15, 0, 1, sec=self.radTprox2)
        neuron.h.pt3dadd(15, 30, 0, 1, sec=self.radTprox2)
        neuron.h.pt3dadd(15, 30, 0, 1, sec=self.radTmed1)
        neuron.h.pt3dadd(15, 45, 0, 1, sec=self.radTmed1)
        neuron.h.pt3dadd(15, 45, 0, 1, sec=self.radTmed2)
        neuron.h.pt3dadd(15, 60, 0, 1, sec=self.radTmed2)
        neuron.h.pt3dadd(15, 60, 0, 1, sec=self.radTdist1)
        neuron.h.pt3dadd(15, 75, 0, 1, sec=self.radTdist1)
        neuron.h.pt3dadd(15, 75, 0, 1, sec=self.radTdist2)
        neuron.h.pt3dadd(15, 90, 0, 1, sec=self.radTdist2)
        neuron.h.pt3dadd(15, 15, 0, 1, sec=self.rad_t1)
        neuron.h.pt3dadd(75, 45, 0, 1, sec=self.rad_t1)
        neuron.h.pt3dadd(15, 45, 0, 1, sec=self.rad_t2)
        neuron.h.pt3dadd(-45, 75, 0, 1, sec=self.rad_t2)
        neuron.h.pt3dadd(15, 75, 0, 1, sec=self.rad_t3)
        neuron.h.pt3dadd(75, 105, 0, 1, sec=self.rad_t3)
        neuron.h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick2)
        neuron.h.pt3dadd(45, 105, 0, 1, sec=self.lm_thick2)
        neuron.h.pt3dadd(45, 105, 0, 1, sec=self.lm_medium2)
        neuron.h.pt3dadd(75, 120, 0, 1, sec=self.lm_medium2)
        neuron.h.pt3dadd(75, 120, 0, 1, sec=self.lm_thin2)
        neuron.h.pt3dadd(105, 135, 0, 1, sec=self.lm_thin2)
       
        neuron.h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick1)
        neuron.h.pt3dadd(-15, 105, 0, 1, sec=self.lm_thick1)
        neuron.h.pt3dadd(-15, 105, 0, 1, sec=self.lm_medium1)
        neuron.h.pt3dadd(-45, 120, 0, 1, sec=self.lm_medium1)
        neuron.h.pt3dadd(-45, 120, 0, 1, sec=self.lm_thin1)
        neuron.h.pt3dadd(-70, 135, 0, 1, sec=self.lm_thin1)

        neuron.h.pt3dadd(0, 0, 0, 1, sec=self.oriprox1)
        neuron.h.pt3dadd(-45, -30, 0, 1, sec=self.oriprox1)
        neuron.h.pt3dadd(-45, -30, 0, 1, sec=self.oridist1_1)
        neuron.h.pt3dadd(-75, -60, 0, 1, sec=self.oridist1_1)
        neuron.h.pt3dadd(-45, -30, 0, 1, sec=self.oridist1_2)
        neuron.h.pt3dadd(-85, -30, 0, 1, sec=self.oridist1_2)
        neuron.h.pt3dadd(15, 0, 0, 1, sec=self.oriprox2)
        neuron.h.pt3dadd(60, -30, 0, 1, sec=self.oriprox2)
        neuron.h.pt3dadd(60, -30, 0, 1, sec=self.oridist2_1)
        neuron.h.pt3dadd(105, -60, 0, 1, sec=self.oridist2_1)
        neuron.h.pt3dadd(60, -30, 0, 1, sec=self.oridist2_2)
        neuron.h.pt3dadd(100, -30, 0, 1, sec=self.oridist2_2)
        neuron.h.pt3dadd(15, 0, 0, 1, sec=self.axon)
        neuron.h.pt3dadd(15, -150, 0, 1, sec=self.axon)
        
    def geometry(self):
        self.soma.L = 10
        self.soma.diam = 10
        self.radTprox1.L = 50
        self.radTprox1.diam = 4
        self.radTprox2.L = 50
        self.radTprox2.diam = 4
        self.radTmed1.L = 50
        self.radTmed1.diam = 3
        self.radTmed2.L = 50
        self.radTmed2.diam = 3
        self.radTdist1.L = 100
        self.radTdist1.diam = 2
        self.radTdist2.L = 100
        self.radTdist2.diam = 2
        self.rad_t1.L = 150
        self.rad_t1.diam = 1
        self.rad_t2.L = 150
        self.rad_t2.diam = 1
        self.rad_t3.L = 150
        self.rad_t3.diam = 1
        self.lm_thick2.L = 100
        self.lm_thick2.diam = 2
        self.lm_medium2.L = 100
        self.lm_medium2.diam = 1.5
        self.lm_thin2.L = 50
        self.lm_thin2.diam = 1
        self.lm_thick1.L = 100
        self.lm_thick1.diam = 2
        self.lm_medium1.L = 100
        self.lm_medium1.diam = 1.5
        self.lm_thin1.L = 50
        self.lm_thin1.diam = 1
        self.oriprox1.L = 100
        self.oriprox1.diam = 2
        self.oridist1_1.L = 200
        self.oridist1_1.diam = 1.5
        self.oridist1_2.L = 200
        self.oridist1_2.diam = 1.5
        self.oriprox2.L = 100
        self.oriprox2.diam = 2
        self.oridist2_1.L = 200
        self.oridist2_1.diam = 1.5
        self.oridist2_2.L = 200
        self.oridist2_2.diam = 1.5
        self.axon.L = 150
        self.axon.diam = 1
        
    def geom_nseg(self):
        for sec in self.sections:
            sec.nseg = 1 + 2*int(sec.L/40)
        
    def load_neuron(self):
        working_dir = os.getcwd()
        os.chdir(mechanisms_path)
        if self.recompile:
            p = run('nrnivmodl')
        if self.add_to_h:
            neuron.load_mechanisms(mechanisms_path)
        os.chdir(working_dir)

        self.make_sections()
        self.connections()
        self.make_lists()
        self.make_basic_shape()
        self.geometry()
        self.geom_nseg()
        self.add_channels()
        self.biophysics()
        
    def add_channels(self):
        for sec in self.sections:
            sec.insert("pas")
            sec.insert("kdr")
            sec.insert("nax")

        for sec in self.somatic:
            sec.insert("kmb")
            sec.insert("kap")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("cal12")
            sec.insert("cal13")
            sec.insert("cav32")

        for sec in self.apical:
            sec.insert("kad")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("cal12")
            sec.insert("cal13")
            sec.insert("cav32")

        for sec in self.basal:
            sec.insert("kad")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("car")
            sec.insert("cal12")
            sec.insert("cal13")
            sec.insert("cav32")

        for sec in self.axonal:
             sec.insert("kmb")
             sec.insert("kap")

    def biophysics(self):
        for sec in self.sections:
            sec.cm = 1
            sec.ena = 50
            sec.ek = -90
        for sec in self.somatic:
            sec.gbar_kap = 0.0074 
            sec.gbar_kmb = 0.001 
            sec.gbar_kdr = 0.0015
            sec.gbar_nax = 0.035 
            sec.gbar_cal12 =  1e-05
            sec.gbar_cal13 =  3e-05
            sec.gbar_can = 2.26e-06
            sec.gbar_cav32 =  4*1.18e-08
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05

        for sec in self.axonal:
            sec.gbar_nax = 0.035 #0.21113423945477339
            sec.gbar_kdr = 0.012
            sec.gbar_kmb = 0.0265
            sec.gbar_kap = 0.164
            sec.Ra = 85.20
            sec.g_pas = 0.00013
            sec.e_pas = -79.92
        for sec in self.apical: 
            sec.gbar_kdr = 0.0043
            sec.gbar_nax = 0.0383
            sec.gbar_cal12 = 2.5e-06
            sec.gbar_cal13 = 2.5e-06
            sec.gbar_can = 1.13e-06
            sec.gbar_cav32 = 1.18e-08
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05
        for sec in self.trunk:
            sec.gbar_kdr = 0.02
            sec.gbar_nax = 0.025
            sec.gbar_cal12 = 1.25e-06
            sec.gbar_cal13 = 1.25e-06
            sec.gbar_can = 2.26e-05
            sec.gbar_cav32 = 2*1.18e-08
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05
        for sec in self.basal:
            sec.gbar_kdr = 0.0043
            sec.gbar_nax = 0.0383
            sec.gbar_cal12 = 1.25e-06
            sec.gbar_cal13 = 2.5e-06
            sec.gbar_can = 1.13e-06
            sec.gbar_cav32 = 4*1.18e-08
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05

        for sec in self.apical+self.basal+self.somatic:
            for seg in sec:
                x = neuron.h.distance(self.soma(0.5), seg)
                value_e_pas = e_pas_dist(x)
                value_ghdbar = ghdbar_dist(x)
                value_gkabar = gkabar_dist(x)
                to_mech = getattr(seg, "pas")
                setattr(to_mech, "e", value_e_pas)
                to_mech = getattr(seg, "hd")
                setattr(to_mech, "gbar", value_ghdbar)
                if "soma" not in sec.name():
                    to_mech = getattr(seg, "kad")
                    setattr(to_mech, "gbar", value_gkabar)
        self.radTprox1.gbar_kad = 0.1
        self.radTprox2.gbar_kad = 0.1
        self.radTmed1.gbar_kad = 0.15
        self.radTmed2.gbar_kad = 0.15
        self.rad_t2.gbar_kad = 0.1
        self.rad_t2.gbar_nax = 0.038
        self.rad_t2.gbar_kdr = 0.002
        self.radTdist1.gbar_kad = 0.2
        self.radTdist2.gbar_kad = 0.2
