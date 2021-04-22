import os
from numpy import pi, exp
from neuron import rxd, h
import neuron
from .ca_params import *
from neuron.units import nM, uM
from subprocess import run
from collections import OrderedDict
# Nomenclature and values adapted from Harris KM, Jensen FE, Tsao BE.
# J Neurosci 1992
#absolute quantities in mM
#reaction rates in uM/ms or 1/ms
my_loc = os.path.dirname(os.path.abspath(__file__))
c_unit = 6.0221409e5
my_loc = os.path.dirname(os.path.abspath(__file__))
mechanisms_path = os.path.join(my_loc, "Mods")

mech_dict = {
    "hd": "ghdbar",
    'kad': "gkabar",
    'kap': "gkabar",
    'kdr': "gkdrbar",
    "kmb": "gbar",
    'na3notrunk': "gbar",
    'nap': "gnabar",
    "cal": "gcalbar",
    "can": "gcanbar",
    "car": "gcabar",
    "cat": "gcatbar",
    "kca": "gbar", # sk - channel
    "cagk": "gbar",
    "nax": "gbar",
    "pas": ["g", "e"], }



SPINE_DIMENSIONS = {
    "mushroom": {
        "head_diam": 1.1,
        "head_len": 0.8,
        "neck_diam": 0.20,
        "neck_len": 0.43,
    },
    "thin": {
        "head_diam": 0.2,
        "head_len": 0.5,
        "neck_diam": 0.1,
        "neck_len": 0.5,
    },
    "stubby": {
        "head_diam": 0.32,
        "head_len": 0.2,
        "neck_diam": 0.32,
        "neck_len": 0.2,
    },
     "generic": {
        "head_diam": 1.,
        "head_len": 1,
        "neck_diam": 0.5,
        "neck_len": 0.5,
    },
    "ball and stick":{
        "head_diam": .5,
        "head_len": .5,
        "neck_diam": 0.2,
        "neck_len": 0.5,
    }
}

def e_pas_dist(x):
    return -65.726902768520958-5*x/150

def ghdbar_dist(x):
    return (1. + 3./100. * x)*1.9042409723832741e-05

def gkabar_dist(x):
    return (15./(1. + exp((300-x)/50)))* 0.012921529390557651

        
def get_ca_init(reg):
    if "ECS" in reg:
        return ca_ECS
    return 1.02e-4


class CA1_PC:
    
    synlist = []
    nclist = []
    spine_dict = OrderedDict()
    sections = []
    heads = []
     
    def cell_filter(self, name, tolist=False):
        out = []
        for sec in self.sections:
            if name in sec.name():
                out.append(sec)
        if len(out) == 1 and tolist is False:
            return out[0]
        return out

    def __init__(self, where_ca=["lm_medium2]"],
                 spine_number=1,
                 where_spines=["lm_medium2"],
                 add_ER=True, buffer_list=["CaM", "Calbindin"]):
        self.load_neuron()
        self.add_ER = add_ER
        if self.add_ER:
            self.ER = OrderedDict()
            self.cyt_er_membrane = OrderedDict()
            self.fc = 0.8 # fraction of cytoplasmic volume
        else:
            self.fc = 1
        
        self.reactions = []
        self.where_spines = []
        self.sections_rxd = []
        if where_spines:
            if not isinstance(where_spines, list):
                where_spines = [where_spines]
            for loc in where_spines:
                self.where_spines.extend(self.cell_filter(loc, tolist=True))
            for spine_loc in self.where_spines:
                for a_spine in range(spine_number):
                    pos = 1/(spine_number+1)*(a_spine+1)
                    head = self.add_head(a_spine, where=spine_loc, position=pos)
                    self.heads.append(head)
            self.compensate_for_spines()
        self.add_calcium(where_rxd=where_ca, buffer_list=buffer_list)
        for head in self.heads:
            self.add_synapse_ampa(head(0.9))
            self.add_synapse_nmda(head(0.9))

    def make_sections(self):
        self.soma = h.Section(name="soma")
        self.radTprox1 = h.Section(name="radTprox1")
        self.radTprox2 = h.Section(name="radTprox2")
        self.radTmed1 = h.Section(name="radTmed1")
        self.radTmed2 = h.Section(name="radTmed2")
        self.radTdist1 = h.Section(name="radTdist1")
        self.radTdist2 = h.Section(name="radTdist2")
        self.lm_thick1 = h.Section(name="lm_thick1")
        self.lm_thick2 = h.Section(name="lm_thick2")
        self.lm_medium1 = h.Section(name="lm_medium1")
        self.lm_medium2 = h.Section(name="lm_medium2")
        self.lm_thin1 = h.Section(name="lm_thin1")
        self.lm_thin2 = h.Section(name="lm_thin2")
        self.axon = h.Section(name="axon")
        self.oriprox1 = h.Section(name="oriprox1")
        self.oriprox2 = h.Section(name="oriprox2")
        self.oridist1_1 = h.Section(name="oridist1_1")
        self.oridist1_2 = h.Section(name="oridist1_2")
        self.oridist2_1 = h.Section(name="oridist2_1")
        self.oridist2_2 = h.Section(name="oridist2_2")
        self.rad_t1 = h.Section(name="rad_t1")
        self.rad_t2 = h.Section(name="rad_t2")
        self.rad_t3 = h.Section(name="rad_t3")
        
    def connections(self):
        self.radTprox1.connect(self.soma(1))
        self.radTprox2.connect(self.radTprox1(1))
        self.radTmed1.connect(self.radTprox2(1))
        self.radTmed2.connect(self.radTmed1(1))
        self.radTdist1.connect(self.radTmed2(1))
        self.radTdist2.connect(self.radTdist1(1))
        self.rad_t1.connect(self.radTprox1(1))
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
        self.trunk_sec_list = [
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
            h.pt3dclear(sec=sec)
        h.pt3dadd(0, 0, 0, 1, sec=self.soma)
        h.pt3dadd(15, 0, 0, 1, sec=self.soma)
        h.pt3dadd(15, 0, 0, 1, sec=self.radTprox1)
        h.pt3dadd(15, 15, 0, 1, sec=self.radTprox1)
        h.pt3dadd(15, 15, 0, 1, sec=self.radTprox2)
        h.pt3dadd(15, 30, 0, 1, sec=self.radTprox2)
        h.pt3dadd(15, 30, 0, 1, sec=self.radTmed1)
        h.pt3dadd(15, 45, 0, 1, sec=self.radTmed1)
        h.pt3dadd(15, 45, 0, 1, sec=self.radTmed2)
        h.pt3dadd(15, 60, 0, 1, sec=self.radTmed2)
        h.pt3dadd(15, 60, 0, 1, sec=self.radTdist1)
        h.pt3dadd(15, 75, 0, 1, sec=self.radTdist1)
        h.pt3dadd(15, 75, 0, 1, sec=self.radTdist2)
        h.pt3dadd(15, 90, 0, 1, sec=self.radTdist2)
        h.pt3dadd(15, 15, 0, 1, sec=self.rad_t1)
        h.pt3dadd(75, 45, 0, 1, sec=self.rad_t1)
        h.pt3dadd(15, 45, 0, 1, sec=self.rad_t2)
        h.pt3dadd(-45, 75, 0, 1, sec=self.rad_t2)
        h.pt3dadd(15, 75, 0, 1, sec=self.rad_t3)
        h.pt3dadd(75, 105, 0, 1, sec=self.rad_t3)
        h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick2)
        h.pt3dadd(45, 105, 0, 1, sec=self.lm_thick2)
        h.pt3dadd(45, 105, 0, 1, sec=self.lm_medium2)
        h.pt3dadd(75, 120, 0, 1, sec=self.lm_medium2)
        h.pt3dadd(75, 120, 0, 1, sec=self.lm_thin2)
        h.pt3dadd(105, 135, 0, 1, sec=self.lm_thin2)
       
        h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick1)
        h.pt3dadd(-15, 105, 0, 1, sec=self.lm_thick1)
        h.pt3dadd(-15, 105, 0, 1, sec=self.lm_medium1)
        h.pt3dadd(-45, 120, 0, 1, sec=self.lm_medium1)
        h.pt3dadd(-45, 120, 0, 1, sec=self.lm_thin1)
        h.pt3dadd(-70, 135, 0, 1, sec=self.lm_thin1)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriprox1)
        h.pt3dadd(-45, -30, 0, 1, sec=self.oriprox1)
        h.pt3dadd(-45, -30, 0, 1, sec=self.oridist1_1)
        h.pt3dadd(-75, -60, 0, 1, sec=self.oridist1_1)
        h.pt3dadd(-45, -30, 0, 1, sec=self.oridist1_2)
        h.pt3dadd(-85, -30, 0, 1, sec=self.oridist1_2)
        h.pt3dadd(15, 0, 0, 1, sec=self.oriprox2)
        h.pt3dadd(60, -30, 0, 1, sec=self.oriprox2)
        h.pt3dadd(60, -30, 0, 1, sec=self.oridist2_1)
        h.pt3dadd(105, -60, 0, 1, sec=self.oridist2_1)
        h.pt3dadd(60, -30, 0, 1, sec=self.oridist2_2)
        h.pt3dadd(100, -30, 0, 1, sec=self.oridist2_2)
        h.pt3dadd(15, 0, 0, 1, sec=self.axon)
        h.pt3dadd(15, -150, 0, 1, sec=self.axon)
        
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
            sec.nseg = 1+2*int(sec.L/40)
        
        
    def load_neuron(self):
        working_dir = os.getcwd()
        os.chdir(mechanisms_path)
        p = run('nrnivmodl')
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
            sec.insert("kca")
            sec.insert("cagk")

        for sec in self.apical:
            sec.insert("kad")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("cal")
            sec.insert("cat")
            sec.insert("kca")
            sec.insert("cagk")

        for sec in self.basal:
            sec.insert("kad")
            sec.insert("hd")
            sec.insert("can")
            sec.insert("cal")
            sec.insert("cat")
            sec.insert("kca")
            sec.insert("cagk")

        for sec in self.axonal:
             sec.insert("kmb")
             sec.insert("kap")

    def biophysics(self):
        h.celsius = 35
        h.v_init = -65
        for sec in self.sections:
            sec.cm = 1
            sec.ena = 50
            sec.ek = -90
        for sec in self.somatic:
            sec.gkabar_kap = 0.0075 
            sec.gbar_kmb = 0.001 
            sec.gkdrbar_kdr = 0.0015 
            sec.gbar_nax = 0.035 
            sec.gcalbar_cal =  0.0005 
            sec.gcanbar_can = 2.2619e-06
            sec.gcatbar_cat =  0.00005 
            sec.gbar_kca = 0.0015 
            sec.gbar_cagk = 4.482e-05
            sec.Ra = 115.3958
            sec.g_pas = 9.0314e-05

        for sec in self.axonal:
            sec.gbar_nax = 0.035 #0.21113423945477339
            sec.gkdrbar_kdr = 0.011664045469379856
            sec.gbar_kmb = 0.026473888790212396
            sec.gkabar_kap = 0.1636942175250268
            sec.Ra = 85.202399381150826
            sec.g_pas = 0.00012898002027660884
            sec.e_pas = -79.917091935442244
        for sec in self.apical: 
            sec.gkdrbar_kdr = 0.0043036502438625682 #0.02
            sec.gbar_nax = 0.038280628170345957 #0.025
            sec.gcalbar_cal = 8.0324964335287e-06
            sec.gcanbar_can = 2.2618914062501833e-06
            sec.gcatbar_cat = 1.184948741542104e-06
            sec.gbar_kca = 9.0311387916396796e-05
            sec.gbar_cagk = 4.4820097108998517e-05
            sec.Ra = 115.3957607556371
            sec.g_pas = 9.031387191839301e-05
        for sec in self.trunk_sec_list:
            sec.gkdrbar_kdr = 0.02
            sec.gbar_nax = 0.025
            sec.gcalbar_cal = 8.0324964335287e-06
            sec.gcanbar_can = 2.2618914062501833e-06
            sec.gcatbar_cat = 1.184948741542104e-06
            sec.gbar_kca = 9.0311387916396796e-05
            sec.gbar_cagk = 4.4820097108998517e-05
            sec.Ra = 115.3957607556371
            sec.g_pas = 9.031387191839301e-05
        for sec in self.basal:
            sec.gkdrbar_kdr = 0.0043036502438625682
            sec.gbar_nax = 0.038280628170345957
            sec.gcalbar_cal = 8.0324964335287e-06
            sec.gcanbar_can = 2.2618914062501833e-06
            sec.gcatbar_cat = 1.184948741542104e-06
            sec.gbar_kca = 9.0311387916396796e-05
            sec.gbar_cagk = 4.4820097108998517e-05
            sec.Ra = 115.3957607556371
            sec.g_pas = 9.031387191839301e-05

        for sec in self.apical+self.basal+self.somatic:
            for seg in sec:
                x = h.distance(self.soma(0.5), seg)
                value_e_pas = e_pas_dist(x)
                value_ghdbar = ghdbar_dist(x)
                value_gkabar = gkabar_dist(x)
                to_mech = getattr(seg, "pas")
                setattr(to_mech, "e", value_e_pas)
                to_mech = getattr(seg, "hd")
                setattr(to_mech, "ghdbar", value_ghdbar)
                if "soma" not in sec.name():
                    to_mech = getattr(seg, "kad")
                    setattr(to_mech, "gkabar", value_gkabar)
        self.radTprox1.gkabar_kad = 0.1
        self.radTprox2.gkabar_kad = 0.1
        self.radTmed1.gkabar_kad = 0.15
        self.radTmed2.gkabar_kad = 0.15
        self.rad_t2.gkabar_kad = 0.1
        self.rad_t2.gbar_nax = 0.038
        self.rad_t2.gkdrbar_kdr = 0.002
        self.radTdist1.gkabar_kad = 0.2
        self.radTdist2.gkabar_kad = 0.2

    def add_synapse_ampa(self, dend):
        syn = h.MyExp2Syn(dend)
        syn.tau1 = 0.5
        syn.tau2 = 3
        syn.e = 0
        syn.gAMPAmax = gAMPA
        self.synlist.append(syn)
        return syn

    def add_synapse_nmda(self, dend):
        syn = h.NMDAca(dend)
        syn.fCa = Ca_per
        syn.tcon = 3
        syn.tcoff = 40
        syn.mgconc = 2
        syn.gamma = 0.08
        syn.gNMDAmax = gNMDA
        self.synlist.append(syn)
        return syn

    def add_spine_mechs(self, spine, segment, spine_mechanisms,
                        spine_params_dict, head=True):
        for mech_name in spine_mechanisms:
            spine.insert(mech_name)
        for seg in spine:
            if "cm" not in spine_params_dict:
                seg.cm = segment.cm
            else:
                seg.cm = spine_params_dict["cm"]
            
            if "ek" not in spine_params_dict:
                seg.ek = segment.ek
            else:
                seg.ek = spine_params_dict["ek"]
            if "Ra" not in spine_params_dict:
                seg.sec.Ra = segment.sec.Ra
            else:
                seg.sec.Ra = spine_params_dict["Ra"]
            for mech_name in spine_mechanisms:
                if mech_name  in spine_params_dict:
                    to_mech = getattr(seg, mech_name)
                    for name in spine_params_dict[mech_name]:
                        if head:
                            if mech_name in head_factor:
                                value = spine_params_dict[mech_name][name]*head_factor
                            else:
                                value = spine_params_dict[mech_name][name]
                        setattr(to_mech, name, value)
                else:
                    gbar_name = mech_dict[mech_name]
                    self.set_val_from_segment(seg, segment, mech_name, gbar_name)

        return spine

    def set_val_from_segment(self, seg, segment, mech_name, gbar_names):
        if not isinstance(gbar_names, list):
            gbar_names = [gbar_names]
        for gbar_name in gbar_names:
            try:
                from_mech = getattr(segment, mech_name)
            except AttributeError:
                continue    
            value = getattr(from_mech, gbar_name)
            to_mech = getattr(seg, mech_name)
            setattr(to_mech, gbar_name, value)
    
    def add_head(self, number, spine_type="ball and stick",
                 where="apical_dendrite[9]",
                 position=0.5, head_mechanisms=["pas", "cal", "can", 
                                                "cat", "kca", "cagk", 'hd',
                                                'kad', 'kap', 'kdr', "kmb",
                                                'nax'],
                 head_params_dict={"Ra":1200, "cad":{"taur": 14, "buffer":20}}):
        if isinstance(where, str): 
            dend = self.cell_filter(where)
        else:
            dend = where
        segment = dend(position)
        head = h.Section(name="%s_%s_%d" % (dend.name(), 'head', number))
        head.L = SPINE_DIMENSIONS[spine_type]["head_len"]
        head.diam = SPINE_DIMENSIONS[spine_type]["head_diam"]
        head.connect(segment)

        self.sections.append(head)
        if where not in self.spine_dict:
            self.spine_dict[where] = []
        self.spine_dict[where] = (head)
        self.add_spine_mechs(head, segment, head_mechanisms,
                             head_params_dict)
        return head

    #Add a number of spines to a dendritic branch.
    #lower conductances, scale resistance and capacitance

    def compensate_for_spines(self):
        pass

    def add_calcium(self, where_rxd=[], buffer_list=["CaM", "Calbindin"]):
        for name in where_rxd:
            if name == "apical":
                self.sections_rxd.extend(self.apical)
            elif name == "basal":
                self.sections_rxd.extend(self.basal)
            elif name == "soma":
                self.sections_rxd.extend(self.soma)
            elif name == "trunk":
                self.sections_rxd.extend(self.trunk_sec_list)
            elif name == "oblique":
                self.sections_rxd.extend(self.oblique)
            
            else:
                self.sections_rxd += self.cell_filter(name, tolist=True)
        self.sections_rxd = list(set(self.sections_rxd))
        self.add_rxd_calcium(buffer_list)
        for sec in self.sections:
            if sec not in self.sections_rxd:
                if h.ismembrane("ca_ion", sec=sec):
                    sec.insert("cad")
                    sec.cao = Ca_Ext
        return

    def _add_rxd_regions(self):
        self.ECS = rxd.Region(self.sections_rxd, name='ECS', nrn_region='o',
                              geometry=rxd.Shell(1, 2))
        self.membrane = OrderedDict()
        self.shells = OrderedDict()
        self.borders = OrderedDict()
        self.dr = OrderedDict() # shell thickness 

        self.factors = OrderedDict()  #  membrane_shell_width/max_radius
        secs_spines = OrderedDict()
        sec_list = self.sections_rxd[:]
        for sec in sec_list:
            if "head" in sec.name():
                sec_list.remove(sec)
        dendrites = sec_list[:]
        all_geom = OrderedDict()
        for sec in dendrites:
            sec_name = sec.name()
            which_dend = sec_name.replace("[", "").replace("]", "")
            if sec not in self.where_spines: #  add first shell with/withour spines
                factor = 2*membrane_shell_width/sec.diam
                self.shells[sec_name] = [rxd.Region(sec, nrn_region='i',
                                                      geometry=rxd.Shell(1-factor, 1),
                                                      name="%s_Shell_0" % which_dend)]
                self.factors[sec_name] = [factor]
                self.membrane[sec_name] = rxd.Region(sec, name='%s_membrane' % which_dend,
                                                 geometry=rxd.membrane())
            
            else:
                #outermost shell is the outermost shell of the dendrite and the spine/spines
                new_secs = self.cell_filter(sec_name, tolist=True)
                secs_spines[sec_name] = new_secs
                all_geom[sec_name] = []
                for new_sec in new_secs:
                    if "head" in new_sec.name():
                        all_geom[sec_name].append(rxd.inside)
                    else:
                        factor = 2*membrane_shell_width/new_sec.diam
                        all_geom[sec_name].append(rxd.Shell(1-factor, 1))
                
        
                self.shells[sec_name] = [rxd.Region(secs_spines[sec_name],
                                                    nrn_region="i",
                                                    name="%s_Shell_0"
                                                    % which_dend,
                                                    geometry=rxd.MultipleGeometry(secs=secs_spines[sec_name],
                                                                                  geos=all_geom[sec_name]))]
                self.factors[sec_name] = [factor]
                self.membrane[sec_name] = rxd.Region(new_secs,
                                                     name='%s_membrane'
                                                     % which_dend,
                                                     geometry=rxd.membrane())
            

        for sec in dendrites:
            sec_name = sec.name()
            which_dend = sec_name.replace("[", "").replace("]", "")
            # add other shells

            factor = self.factors[sec_name][0]
            last_shell = False
            i = 0
            while True:
                new_factor = 2*factor
                i = i + 1
                inner = 1-(sum(self.factors[sec_name])+new_factor)
                outer = 1-sum(self.factors[sec_name])
                self.shells[sec_name].append(rxd.Region(sec, nrn_region='i',
                                                        geometry=rxd.Shell(inner,
                                                                           outer),
                                                        name="%s_Shell_%d" % (which_dend, i)))
                if sec_name not in self.borders:
                    self.borders[sec_name] = []
                self.borders[sec_name].append(rxd.Region(sec, name='%s_Border_%d' % (which_dend, i-1),
                                                         geometry=rxd.ScalableBorder(
                                                             diam_scale=outer)))
                self.factors[sec_name].append(new_factor)
                if last_shell:
                    break
               
                if sum(self.factors[sec_name]) + 2*new_factor >= self.fc:
                    last_shell = True
                    factor = (self.fc - sum(self.factors[sec_name]))/2
                else:
                    factor = new_factor                    

            if self.add_ER:
                last_shell_diam = 1 - sum(self.factors[sec_name])
                self.ER[sec_name] = rxd.Region(sec, nrn_region='i',
                                               geometry=rxd.Shell(0, last_shell_diam),
                                               name="%s_ER" % which_dend)

                self.cyt_er_membrane[sec_name] = rxd.Region(sec,
                                                            name='%s_ER_membrane' % which_dend,
                                                            geometry=rxd.ScalableBorder(
                                                                diam_scale=last_shell_diam))

            self.dr[sec_name] = []
            for j, f in enumerate(self.factors[sec_name][:-1]):
                name = "%s_dr_%d" % (which_dend, j)
                self.dr[sec_name].append(rxd.Parameter(self.borders[sec_name][j], name=name, 
                                               value=lambda nd: nd.segment.diam/2/f))
            
        self.shell_list = []
        for key in self.shells.keys():
            self.shell_list.extend(self.shells[key])

        self.membrane_list = []
        for key in self.membrane.keys():
            self.membrane_list.append(self.membrane[key])

        if self.add_ER:
            self.ER_regions = []
            self.cyt_er_membrane_list = []
            self.shells_to_ER = []
            for sec_name in self.ER.keys():
                self.ER_regions.append(self.ER[sec_name])
                self.cyt_er_membrane_list.append(self.cyt_er_membrane[sec_name])
                self.shells_to_ER.append(self.shells[sec_name][-1])
                


    def _add_species(self, buffer_list):
        regions = self.shell_list + [self.ECS]
        self.ca = rxd.Species(regions, d=caDiff, name='ca', charge=2,
                              initial=lambda nd: get_ca_init(nd.region.name),
                              atolscale=1e-8)
        
        if self.add_ER:
            self.ip3 = rxd.Species(self.shell_list, d=ip3Diff, initial=ip3_init,
                                   atolscale=1e-8)
            self.ca_ER = rxd.Species(self.ER_regions, d=caDiff, name='ca_er', charge=2,
                                     initial=ca_init_ER,
                                     atolscale=1e-8)
            self.calr = rxd.Species(self.ER_regions, initial=calr_tot-calr_bound)
            self.calrca = rxd.Species(self.ER_regions, initial=calr_bound)
        self.add_buffers(buffer_list)
        self.add_surface_pumps()


    def pump_density(self, node, value):

        area = node.sec(node.x).area()
        sec_name = node.sec.name()
        if "head" in sec_name:
            sec_name = node.sec.name().split("_head")[0]
            return value[sec_name]*area*head_factor
        factor = self.factors[sec_name][0]
        shell_vol = pi/4*node.sec.L*node.sec.diam**2*(1-(1-factor)**2)/vol_c
        return value[sec_name]*area*shell_vol

    def add_surface_pumps(self):
        gncx_spine = OrderedDict()
        gpmca_spine = OrderedDict()
        gncx_dend = OrderedDict()
        gpmca_dend = OrderedDict()
        for key in gncx_spine_total.keys():
            gncx_spine[key] = gncx_spine_total[key] - gncx_spine_bound[key]
            gpmca_spine[key] = gpmca_spine_total[key] - gpmca_spine_bound[key]
            gncx_dend[key] =  gncx_dend_total[key] - gncx_dend_bound[key]
            gpmca_dend[key] =  gpmca_dend_total[key] - gpmca_dend_bound[key]
        self.ncx = rxd.Species(self.membrane_list, name='ncx', initial=lambda nd:
                               self.pump_density(nd, gncx_spine)
                               if "head" in nd.sec.name()
                               else self.pump_density(nd, gncx_dend))
        self.ncxca = rxd.Species(self.membrane_list, name='ncxca',
                                 initial=lambda nd:
                                 self.pump_density(nd, gncx_spine_bound)
                                 if "head" in nd.sec.name()
                                 else self.pump_density(nd, gncx_dend_bound))
        self.pmca = rxd.Species(self.membrane_list, name='pmca', initial=lambda nd:
                                self.pump_density(nd, gpmca_spine)
                                if "head" in nd.sec.name()
                                else self.pump_density(nd, gpmca_dend))
        self.pmcaca = rxd.Species(self.membrane_list, name='pmcaca',
                                  initial=lambda nd:
                                  self.pump_density(nd, gpmca_spine_bound)
                                  if "head" in nd.sec.name()
                                  else self.pump_density(nd, gpmca_dend_bound))
        
        
    def add_serca(self):
        # SERCA pump: pumps ca from cyt -> ER
        # Vmax = 500 uM/s = 0.5 mM/s, Kserca = 0.4 uM = 3 e-4 mM (in agreement with the sig. pathways model)
    
        self.gserca = rxd.Parameter(self.cyt_er_membrane_list, initial=gserca)
        #Serca from the sig paths model
        self.serca = []
        for key in self.ER.keys():
            er_region = self.ER[key]
            ca_region = self.shells[key][-1]
            membrane = self.cyt_er_membrane[key]
            self.serca.append(rxd.MultiCompartmentReaction(self.ca[ca_region]> self.ca_ER[er_region],
                                                           self.gserca*kcat_Serca/((Kserca / self.ca[ca_region]) ** 2 + 1),
                                                           membrane=membrane, custom_dynamics=True))

    def add_leak(self):
        # leak channel: bidirectional ca flow btwn cyt <> ER
        self.gleak = rxd.Parameter(self.cyt_er_membrane_list, initial=gleak)
        self.leak = []
        for key in self.ER.keys():
            er_region =	self.ER[key]
            ca_region =	self.shells[key][-1]
            membrane = self.cyt_er_membrane[key]
            self.leak.append(rxd.MultiCompartmentReaction(self.ca_ER[er_region], self.ca[ca_region],
                                                          self.gleak, self.gleak,
                                                          membrane=membrane))

    def add_ip3r(self):
        """
        hopefully corresponding to Neymotin's model however extremely difficult to judge
        Apparently (triple checked) Neymotin used Wagner et al 2004's model and this
        implementation corresponds to Wagner and hopefully is correct
        """
        self.gip3r = rxd.Parameter(self.cyt_er_membrane_list, initial=gip3r)
        self.m = []
        self.n = []
        self.h_inf = []
        self.h_gate = rxd.State(self.cyt_er_membrane_list, initial=0.8)
        self.ip3rg = []
        self.ip3r = []
        for i, key in enumerate(self.ER.keys()):
            er_region =	self.ER[key]
            ca_region =	self.shells[key][-1]
            membrane = self.cyt_er_membrane[key]
            self.m.append(self.ip3[ca_region]/(self.ip3[ca_region] + Kip3))
            self.n.append(self.ca[ca_region]/(self.ca[ca_region]  + Kact))
            self.h_inf.append(k_inh/(self.ca[ca_region] + k_inh))
            minf = self.m[i] * self.n[i]
        

            self.ip3rg.append(rxd.Rate(self.h_gate[membrane], (self.h_inf[i] - self.h_gate[membrane]) / ip3rtau))
            k = self.gip3r[membrane] * (minf * self.h_gate[membrane]) ** 3

            self.ip3r.append(rxd.MultiCompartmentReaction(self.ca_ER[er_region],
                                                         self.ca[ca_region],
                                                         k, k, membrane=membrane))
        
        # IP3 receptor gating



    def add_calmodulin(self):
        self.cam = rxd.Species(self.shell_list, d=camDiff,
                               initial=calmodulin_tot-camn-camc,
                               name='CaM', charge=0, atolscale=1e-8)
        self.camn = rxd.Species(self.shell_list, d=camDiff, initial=camn, name='CaMN',
                                charge=0, atolscale=1e-8)
        self.camc = rxd.Species(self.shell_list, d=camDiff, initial=camc, name='CaMC',
                                charge=0, atolscale=1e-8)
        self.buffers["CaM"] = [self.cam, self.camn, self.camc]

    def add_buffers(self, buffer_names):
        self.buffers = OrderedDict()
        self.indicator = None
        for name in buffer_names:
            if name == "CaM":
                self.add_calmodulin()
            if name == "Calbindin":
                self.calb = rxd.Species(self.shell_list, d=calbDiff,
                                        initial=calbindin_tot-calbca,
                                        name='Calbindin',
                                        charge=0, atolscale=1e-8)
                self.calbca = rxd.Species(self.shell_list, d=calbDiff,
                                          initial=calbca,
                                          name='CalbindinCa',
                                          charge=0, atolscale=1e-8)
                self.buffers["Calb"] = [self.calb, self.calbca]

            if name == "Mg Green":
                self.indicator = rxd.Species(self.shell_list,
                                             initial=tot_magnesium_green_BS -
                                             magnesium_green_bound, d=mggreenDiff,
                                             name='MgGreen',
                                             charge=0, atolscale=1e-8)
                self.indicator_ca = rxd.Species(self.shell_list,
                                             initial=magnesium_green_bound,
                                             name='MgGreenCa', d=mggreenDiff,
                                                charge=0, atolscale=1e-8)
                self.buffers["Mg Green"] = [self.indicator, self.indicator_ca]

        self.fixed = rxd.Species(self.shell_list,
                                 initial=fixed_buffer_tot-fixed_buffer_ca,
                                 name='FixedBuffer',
                                 charge=0, atolscale=1e-8)
        self.fixedca = rxd.Species(self.shell_list,
                                   initial=fixed_buffer_ca,
                                   name='FixedBufferCa',
                                   charge=0, atolscale=1e-8)


        
    def add_buffer_reactions(self, buffer_list):
        fixed_rxn = rxd.Reaction(self.fixed + self.ca, self.fixedca, kf_fixed_b,
                                 kb_fixed_b)
        self.reactions.append(fixed_rxn)
        if "CaM" in buffer_list:
            rn = rxd.Reaction(self.cam + self.ca, self.camn, kf_camn,
                              kb_camn)
            rc = rxd.Reaction(self.cam + self.ca, self.camc, kf_camc, kb_camc)
            self.reactions.extend([rn, rc])
            
        if "Calbindin" in buffer_list:
            calb_rxn = rxd.Reaction(self.calb + self.ca, self.calbca,
                                    kf_calbindin,
                                    kb_calbindin)
            self.reactions.append(calb_rxn)
        if "Mg Green" in buffer_list:
            mg_green_binding = rxd.Reaction(self.indicator + self.ca,
                                            self.indicator_ca,
                                            kf_magnesium_green,
                                            kb_magnesium_green)
            self.reactions.append(mg_green_binding)
        if self.add_ER:
            rxn1 = rxd.Reaction(self.calr + self.ca_ER,
                                self.calrca, kf_calr, kb_calr)
            self.reactions.append(rxn1)


    def add_surface_pump_reactions(self):
        for key in self.shells:
            membrane_shell = self.shells[key][0]
            membrane = self.membrane[key]
            self.pmca_1_r = rxd.MultiCompartmentReaction(self.ca[membrane_shell] +
                                                         self.pmca[membrane],
                                                         self.pmcaca[membrane],
                                                         kf_pmca*c_unit,
                                                         kb_pmca*c_unit,
                                                         membrane=membrane)
            self.pmca_2_r = rxd.MultiCompartmentReaction(self.pmcaca[membrane],
                                                         self.pmca[membrane] +
                                                         self.ca[self.ECS],
                                                         kcat_pmca*c_unit,
                                                         membrane=membrane)
            self.ncx_1_r = rxd.MultiCompartmentReaction(self.ca[membrane_shell] +
                                                        self.ncx[membrane],
                                                        self.ncxca[membrane],
                                                        kf_ncx*c_unit,
                                                        kb_ncx*c_unit,
                                                        membrane=membrane)
            self.ncx_2_r = rxd.MultiCompartmentReaction(self.ncxca[membrane],
                                                        self.ncx[membrane] +
                                                        self.ca[self.ECS],
                                                        kcat_ncx*c_unit,
                                                        membrane=membrane)
            
    def _add_diffusion(self):
        self.diffusions = []
        for sec_name in self.shells.keys():
            for i, shell in enumerate(self.shells[sec_name][:-1]):
                rxn = rxd.MultiCompartmentReaction(self.ca[shell],
                                                   self.ca[self.shells[sec_name][i+1]], 
                                                   c_unit*caDiff/self.dr[sec_name][i], 
                                                   c_unit*caDiff/self.dr[sec_name][i],
                                                   border=self.borders[sec_name][i])
                self.diffusions.append(rxn)
                for buf_name in self.buffers.keys():
                    dif_const = diffusions[buf_name]
                    for buffs in self.buffers[buf_name]:
                        rxn = rxd.MultiCompartmentReaction(buffs[shell],
                                                           buffs[self.shells[sec_name][i+1]], 
                                                           c_unit*dif_const/self.dr[sec_name][i], 
                                                           c_unit*dif_const/self.dr[sec_name][i],
                                                           border=self.borders[sec_name][i])
                        self.diffusions.append(rxn)

    def add_rxd_calcium(self, buffer_names):
        if self.sections_rxd:
            for sec in self.sections_rxd:
                 if "head" in sec.name():
                     sec.nseg = n_seg
            self._add_rxd_regions()
            self._add_species(buffer_names)
            self._add_diffusion()
            if self.add_ER:
                 self.add_serca()
                 self.add_ip3r()
                 self.add_leak()
            self.add_buffer_reactions(buffer_names)
            self.add_surface_pump_reactions()
     

