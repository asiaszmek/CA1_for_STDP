import os
from subprocess import run

from numpy import exp
import neuron


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
        self.radTprox = neuron.h.Section(name="radTprox")
        self.radTmed = neuron.h.Section(name="radTmed")
        self.radTdist = neuron.h.Section(name="radTdist")
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
        self.radTprox.connect(self.soma(1))
        self.radTmed.connect(self.radTprox(1))
        self.radTdist.connect(self.radTmed(1))
        self.rad_t1.connect(self.radTprox(0.5))
        self.rad_t2.connect(self.radTmed(0.5))
        self.rad_t3.connect(self.radTdist(0.5))
        
        self.lm_thick2.connect(self.radTdist(1))
        self.lm_medium2.connect(self.lm_thick2(1))
        self.lm_thin2.connect(self.lm_medium2(1))
        self.lm_thick1.connect(self.radTdist(1))
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
        self.sections = [self.soma,  self.radTprox,
                         self.radTmed, self.radTdist,
                         self.rad_t1, self.rad_t2, self.rad_t3,
                         self.lm_thick2, self.lm_medium2, self.lm_thin2,
                         self.lm_thick1,self.lm_medium1, self.lm_thin1,
                         self.oriprox1, self.oridist1_1, self.oridist1_2,
                         self.oriprox2, self.oridist2_1, self.oridist2_2,
                         self.axon,]

        self.somatic = [self.soma]
        self.axonal = [self.axon]
        self.apical = [self.radTprox, 
                       self.radTmed,
                       self.radTdist, self.rad_t1, self.rad_t2, self.rad_t3,
                       self.lm_thick2, self.lm_medium2,
                       self.lm_thin2, self.lm_thick1, self.lm_medium1,
                       self.lm_thin1]
        self.basal = [self.oriprox1, self.oridist1_1, self.oridist1_2,
                      self.oriprox2, self.oridist2_1, self.oridist2_2]
        self.trunk = [
            self.radTprox,
            self.radTmed,
            self.radTdist,
        ]
        self.oblique = [self.rad_t2]

    def make_basic_shape(self):
        for sec in self.sections:
            neuron.h.pt3dclear(sec=sec)
        neuron.h.pt3dadd(0, 0, 0, 1, sec=self.soma)
        neuron.h.pt3dadd(15, 0, 0, 1, sec=self.soma)
        neuron.h.pt3dadd(15, 0, 0, 1, sec=self.radTprox)
        neuron.h.pt3dadd(15, 30, 0, 1, sec=self.radTprox)
        neuron.h.pt3dadd(15, 30, 0, 1, sec=self.radTmed)
        neuron.h.pt3dadd(15, 60, 0, 1, sec=self.radTmed)
        neuron.h.pt3dadd(15, 60, 0, 1, sec=self.radTdist)
        neuron.h.pt3dadd(15, 90, 0, 1, sec=self.radTdist)
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
        self.radTprox.L = 100
        self.radTprox.diam = 4
        self.radTmed.L = 100
        self.radTmed.diam = 3
        self.radTdist.L = 200
        self.radTdist.diam = 2
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
            sec.insert("cal")
            sec.insert("cat")

        for sec in self.apical:
            sec.insert("kad")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("cal")
            sec.insert("cat")

        for sec in self.basal:
            sec.insert("kad")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("cal")
            sec.insert("cat")

        for sec in self.axonal:
             sec.insert("kmb")
             sec.insert("kap")

    def biophysics(self):
        for sec in self.sections:
            sec.cm = 1
            sec.ena = 50
            sec.ek = -90
            sec.e_pas = -65
        for sec in self.somatic:
            sec.gbar_kap = 0.0075
            sec.gbar_kmb = 0.001 
            sec.gbar_kdr = 0.0015
            sec.gbar_nax = 0.035 
            sec.gbar_cal =  0.0005
            sec.gbar_can = 2.26e-06
            sec.gbar_cat =  0.00005
            sec.gbar_hd = 1.9e-5
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05
        for sec in self.axonal:
            sec.gbar_nax = 0.035 #0.21113423945477339
            sec.gbar_kdr = 0.012
            sec.gbar_kmb = 0.0265
            sec.gbar_kap = 0.164
            sec.Ra = 85.20
            sec.g_pas = 0.00013
            sec.e_pas = -79.9
        for sec in self.apical: 
            sec.gbar_kdr = 0.0043
            sec.gbar_nax = 0.0383
            sec.gbar_cal = 8.03e-06
            sec.gbar_can = 2.26e-06
            sec.gbar_cat = 1.185e-06
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05
            sec.gbar_hd = 1.9e-5*10
        for sec in self.trunk:
            sec.gbar_kdr = 0.02
            sec.gbar_nax = 0.025
            sec.gbar_cal = 8.03e-06
            sec.gbar_can = 2.26e-06
            sec.gbar_cat = 1.185e-06
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05
        for sec in self.basal:
            sec.gbar_kdr = 0.0043
            sec.gbar_nax = 0.0383
            sec.gbar_cal = 8.03e-06
            sec.gbar_can = 2.26e-06
            sec.gbar_cat = 1.185e-06
            sec.Ra = 115.4
            sec.g_pas = 9.03e-05
            sec.gbar_hd = 1.9e-5*5
        for sec in self.apical:
            if sec.name() in ["radTprox", "rad_t2"]:
                sec.gbar_kad = 0.1
            elif sec.name() in ["redTmed", "rad_t1"]:
                sec.gbar_kad = 0.15
            else:
                sec.gbar_kad = 0.2
        for sec in self.basal:
            if sec.name() in ["oriprox1", "oriprox2"]:
                sec.gbar_kad = 0.002
            else:
                sec.gbar_kad = 0.01

        self.rad_t2.gbar_nax = 0.038
        self.rad_t2.gbar_kdr = 0.002
