import os
from numpy import pi
from neuron import rxd, h
from .ca_params import *
from neuron.units import nM, uM
# Nomenclature and values adapted from Harris KM, Jensen FE, Tsao BE.
# J Neurosci 1992
#absolute quantities in mM
#reaction rates in uM/ms or 1/ms
my_loc = os.path.dirname(os.path.abspath(__file__))
hoc_loc = os.path.join(my_loc, 'nrnhoc')

c_unit = 6.0221409e5


mech_dict = {
    "h": "gbar",
    'kad': "gkabar",
    'kap': "gkabar",
    'kdr': "gkdrbar",
    'na3notrunk': "gbar",
    'nap': "gnabar",
    "cal": "gcalbar",
    "calH": "gcalbar",
    "car": "gcabar",
    "cat": "gcatbar",
    "kca": "gbar", # sk - channel
    "mykca": "gkbar",
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

template_path = os.path.join(my_loc, "template",
                             "load_templates.hoc")
morphology_path = os.path.join(my_loc, "morphology",
                             "n123")
lib_path = os.path.join(my_loc, "lib")
generic_path = os.path.join(my_loc, "nrnhoc")

CA1_Bianchi_exec = """
load_file("stdrun.hoc")
load_file("%s")


objref econ  
//////// load needed templates////////////


show_errs=1
debug_lev=1
econ=new ExperimentControl(show_errs, debug_lev)
econ.self_define(econ)

econ.morphology_dir = "%s"       // set location for morphology files
econ.add_lib_dir("Terrence", "%s")   // set location for library files
econ.generic_dir    = "%s"             // set location for cell-setup file 

v_init=-70

econ.xopen_geometry_dependent("cell")
econ.xopen_geometry_dependent("cell-analysis")
cell_analysis(econ)
econ.xopen_generic("cell-setup")
cell_setup(econ)


objref oblique_tip_dendrites, oblique_dendrites
oblique_tip_dendrites = new SectionList()

forsec apical_tip_list {oblique_tip_dendrites.append()}
forsec apical_tip_list_addendum {oblique_tip_dendrites.append()}

oblique_dendrites = new SectionList()

for i = 0, 59 {
	forsec pl[i] { oblique_dendrites.append() }
}
oblique_dendrites.unique()
oblique_dendrites.remove(apical_trunk_list)
"""



        
def get_ca_init(reg):
    if "ECS" in reg:
        return ca_ECS
    return 1.02e-4




class CA1_PC_cAC:
    
    synlist = []
    nclist = []
    spine_dict = {}
    sections = []
    heads = []
    fc = 0.8 # fraction of cytoplasmic volume
    fe = 1-fc # fraction of ER volume    
    def cell_filter(self, name, tolist=False):
        out = []
        for sec in self.sections:
            if name in sec.name():
                out.append(sec)
        if len(out) == 1 and tolist is False:
            return out[0]
        return out

    def __init__(self, where_ca=["apical_dendrite[9]"],
                 spine_number=1,
                 where_spines=["apical_dendrite[9]"],
                 add_ER=True, buffer_list=["CaM", "Calbindin"]):
        self.load_neuron()
        h.distance(sec=self.soma[0])
        self.add_ER = add_ER
        if self.add_ER:
            self.ER = {}
            self.cyt_er_membrane = {}
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
            
    def load_neuron(self):
        h.execute1(CA1_Bianchi_exec % (template_path,
                                       morphology_path,
                                       lib_path,
                                       generic_path))
        for sec in h.allsec():
            self.sections.append(sec)
        self.axon = []
        self.soma = []
        self.apical_non_trunk_list = []
        self.apical_tip_list = []
        self.apical_tip_list_addendum = []
        self.apical_trunk_list = []
        self.basal_tree_list = []
        self.peri_trunk_list = []
        self.oblique_dendrites = []
        self.oblique_tip_dendrites = []
        for sec in h.axon_sec_list:
            self.axon.append(sec)
        for sec in h.soma:
            self.soma.append(sec)
        for sec in h.apical_non_trunk_list:
            self.apical_non_trunk_list.append(sec)
        for sec in h.apical_trunk_list:
            self.apical_trunk_list.append(sec)
        for sec in h.apical_tip_list:
            self.apical_tip_list.append(sec)
        for sec in h.apical_tip_list_addendum:
            self.apical_tip_list_addendum.append(sec)
        for sec in h.basal_tree_list:
            self.basal_tree_list.append(sec)
        for sec in h.peri_trunk_list:
            self.peri_trunk_list.append(sec)
        for sec in h.oblique_dendrites:
            self.oblique_dendrites.append(sec)
        for sec in h.oblique_tip_dendrites:
            self.oblique_tip_dendrites.append(sec)
        
        
        
    def add_synapse_ltpltd(self, dend):
        syn = h.ltpltd(dend)
        self.synlist.append(syn)
        return syn


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
            from_mech = getattr(segment, mech_name)
            value = getattr(from_mech, gbar_name)
            to_mech = getattr(seg, mech_name)
            setattr(to_mech, gbar_name, value)
    
    def add_head(self, number, spine_type="ball and stick", where="apical_dendrite[9]",
                 position=0.5, head_mechanisms=["pas", "calH", "car", 
                                                "cat", "kca", "mykca", 'h',
                                                'kad', 'kap', 'kca', 'kdr',
                                                'na3notrunk', 'nap'],
                 head_params_dict={"Ra":1200, "cacum":{"taur": 14, "buffer":20}}):
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
            self.sections_rxd += self.cell_filter(name, tolist=True)
        
        self.add_rxd_calcium(buffer_list)
        for sec in self.sections:
            if sec in self.sections_rxd:
                for mech in sec.psection()["density_mechs"]:
                    if "ica" in sec.psection()["density_mechs"][mech]:
                        for seg in sec:
                            to_mech = getattr(seg, mech)
                            value = getattr(to_mech, mech_dict[mech])
                            setattr(to_mech, mech_dict[mech], value*ca_i_f)

            elif h.ismembrane("ca_ion", sec=sec):
                sec.insert("cacum")
                sec.cao = Ca_Ext
        return

    def _add_rxd_regions(self):
        self.ECS = rxd.Region(self.sections_rxd, name='ECS', nrn_region='o',
                              geometry=rxd.Shell(1, 2))
        self.membrane = {}
        self.shells = {}
        self.borders = {}
        self.dr = {} # shell thickness 
        #max_radius = max([sec.diam/2 for sec in sec_list])
        self.factors = {}  #  membrane_shell_width/max_radius
        secs_spines = {}
        sec_list = self.sections_rxd[:]
        for sec in sec_list:
            if "head" in sec.name():
                sec_list.remove(sec)
        dendrites = [[sec] for sec in sec_list]
        all_geom = {}
        for sec in dendrites:
            sec_name = sec[0].name()
            which_dend = sec_name.replace("[", "").replace("]", "")
            if sec[0] not in self.where_spines: #  add first shell with/withour spines
                factor = 2*membrane_shell_width/sec[0].diam
                self.shells[sec_name] = [rxd.Region(sec, nrn_region='i',
                                                      geometry=rxd.Shell(1-factor, 1),
                                                      name="Shell_%s_0" % which_dend)]
                self.factors[sec_name] = [factor]
                self.membrane[sec_name] = rxd.Region(sec, name='membrane_%s' % which_dend,
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
                                                    name="Shell_%s_0"
                                                    % which_dend,
                                                    geometry=rxd.MultipleGeometry(secs=secs_spines[sec_name],
                                                                                  geos=all_geom[sec_name]))]
                self.factors[sec_name] = [factor]
                self.membrane[sec_name] = rxd.Region(new_secs, name='membrane_%s' % which_dend,
                                                 geometry=rxd.membrane())
            

        for sec in dendrites:
            sec_name = sec[0].name()           
            which_dend = sec_name.replace("[", "").replace("]", "")
            # add other shells
            i = 1
            factor = self.factors[sec_name][0]
            while True:
                new_factor = 2*factor
                self.factors[sec_name].append(new_factor)
                self.shells[sec_name].append(rxd.Region(sec, nrn_region='i',
                                                        geometry=rxd.Shell(1-(factor+new_factor), 1-factor),
                                                        name="Shell_%s_%d" % (which_dend, i)))
                if sec_name not in self.borders:
                    self.borders[sec_name] = []
                self.borders[sec_name].append(rxd.Region(sec, name='Border_%s_%d' % (which_dend, i-1),
                                                         geometry=rxd.ScalableBorder(
                                                             diam_scale=1-factor)))
                factor = new_factor
                i = i + 1
                if sum(self.factors[sec_name]) + new_factor >= self.fc:
                    break
            
            last_shell_diam = 1 - sum(self.factors[sec_name])
            if self.add_ER:
                self.ER[sec_name] = rxd.Region(sec, nrn_region='i',
                                               geometry=rxd.Shell(0, last_shell_diam),
                                               name="ER_%s" % which_dend)

                self.cyt_er_membrane[sec_name] = rxd.Region(sec,
                                                            name='ER_membrane_%s' % which_dend,
                                                            geometry=rxd.ScalableBorder(
                                                                diam_scale=last_shell_diam))
            else:
                self.shells[sec_name].append(rxd.Region(sec, nrn_region='i',
                                              geometry=rxd.Shell(0, last_shell_diam),
                                              name="Shell_%s_%d" % (which_dend, i)))
                self.borders[sec_name].append(rxd.Region(sec,
                                               name='Border_%s_%d' % (which_dend, i-1),
                                               geometry=rxd.ScalableBorder(
                                                diam_scale=last_shell_diam)))
                self.factors[sec_name].append(last_shell_diam)

            self.dr[sec_name] = []
            for j, f in enumerate(self.factors[sec_name][:-1]):
                name = "dr_%s_%d" % (which_dend, j)
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
                self.ER_regions.append(self.ER[key])
                self.cyt_er_membrane_list.append(self.cyt_er_membrane[key])
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
        if "head" in  node.sec.name():
            return value*area*head_factor
        factor = self.factors[node.sec.name()][0]
        shell_vol = pi/4*node.sec.L*node.sec.diam**2*(1-(1-factor)**2)/vol_c
        return value*area*shell_vol

    def add_surface_pumps(self):
        gncx_spine = ca_i_f*(gncx_spine_total - gncx_spine_bound)
        gpmca_spine = ca_i_f*(gpmca_spine_total - gpmca_spine_bound)
        gncx_dend =  ca_i_f*(gncx_dend_total - gncx_dend_bound)
        gpmca_dend =  ca_i_f*(gpmca_dend_total - gpmca_dend_bound)
        gncxca_spine = ca_i_f*gncx_spine_bound
        gncxca_dend = ca_i_f*gncx_dend_bound
        gpmcaca_spine = ca_i_f*gpmca_spine_bound
        gpmcaca_dend = ca_i_f*gpmca_dend_bound
        self.ncx = rxd.Species(self.membrane_list, name='ncx', initial=lambda nd:
                               self.pump_density(nd, gncx_spine)
                               if "head" in nd.sec.name()
                               else self.pump_density(nd, gncx_dend))
        self.ncxca = rxd.Species(self.membrane_list, name='ncxca',
                                 initial=lambda nd:
                                 self.pump_density(nd, gncxca_spine)
                                 if "head" in nd.sec.name()
                                 else self.pump_density(nd, gncxca_dend))
        self.pmca = rxd.Species(self.membrane_list, name='pmca', initial=lambda nd:
                                self.pump_density(nd, gpmca_spine)
                                if "head" in nd.sec.name()
                                else self.pump_density(nd, gpmca_dend))
        self.pmcaca = rxd.Species(self.membrane_list, name='pmcaca',
                                  initial=lambda nd:
                                  self.pump_density(nd, gpmcaca_spine)
                                  if "head" in nd.sec.name()
                                  else self.pump_density(nd, gpmcaca_dend))
        
        
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
        self.buffers = {}
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
     

