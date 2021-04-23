import os
from collections import OrderedDict

from numpy import pi, exp

from neuron import rxd, h
import neuron
from neuron.units import nM, uM
from .basal_cell import CA1_PC_basal
# Nomenclature and values adapted from Harris KM, Jensen FE, Tsao BE.
# J Neurosci 1992
#absolute quantities in mM
#reaction rates in uM/ms or 1/ms
my_loc = os.path.dirname(os.path.abspath(__file__))
params_file = os.path.join(my_loc, "ca_params.py")
c_unit = 6.0221409e5
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
        
class CA1_PC:
    
    synlist = []
    nclist = []
    spine_dict = OrderedDict()
    heads = []
    sections = []
    def cell_filter(self, name, tolist=False):
        out = []
        for sec in self.sections:
            if name in sec.name():
                out.append(sec)
        if len(out) == 1 and tolist is False:
            return out[0]
        return out

    def __init__(self, model=None, params=params_file, where_ca=["lm_medium2]"],
                 spine_number=1,
                 where_spines=["lm_medium2"],
                 add_ER=True, buffer_list=["CaM", "Calbindin"]):
        if model is None:
            model = CA1_PC_basal()
        for sec in model.sections:
            self.sections.append(sec)
        try:
            self.apical = model.apical
        except AttributeError:
            self.apical = []
        try:
            self.somatic = model.somatic
        except AttributeError:
            self.somatic = self.cell_filter("soma", tolist=True)
        try:
            self.axonal = model.axonal
        except AttributeError:
            self.axonal = self.cell_filter("axon", tolist=True)
        try:
            self.basal = model.basal
        except AttributeError:
            self.basal = []
        try:
            self.trunk = model.trunk
        except AttributeError:
            self.trunk = []
        try:
            self.oblique = model.oblique
        except AttributeError:
            self.oblique = []

        #read in parameters from file
        self.params = {}
        exec(open(params).read(), {}, self.params)

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

    def get_ca_init(self, reg):
        if "ECS" in reg:
            return self.params["ca_ECS"]
        return 1.02e-4

            
    def add_synapse_ampa(self, dend):
        syn = h.MyExp2Syn(dend)
        syn.tau1 = 0.5
        syn.tau2 = 3
        syn.e = 0
        syn.gAMPAmax = self.params["gAMPA"]
        self.synlist.append(syn)
        return syn

    def add_synapse_nmda(self, dend):
        syn = h.NMDAca(dend)
        syn.fCa = self.params["Ca_per"]
        syn.tcon = 3
        syn.tcoff = 40
        syn.mgconc = 2
        syn.gamma = 0.08
        syn.gNMDAmax = self.params["gNMDA"]
        self.synlist.append(syn)
        return syn

    def add_spine_mechs(self, spine, segment, spine_mechanisms,
                        spine_params_dict, head=True):
        for mech_name in spine_mechanisms:
            spine.insert(mech_name)
        head_factor = self.params["head_factor"]
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
                self.sections_rxd.extend(self.trunk)
            elif name == "oblique":
                self.sections_rxd.extend(self.oblique)
            else:
                self.sections_rxd += self.cell_filter(name, tolist=True)
        for sec in self.where_spines:
            if sec in self.sections_rxd:
                sec_list = self.cell_filter(sec.name(), tolist=True)
                self.sections_rxd.extend(sec_list)

        self.sections_rxd = list(set(self.sections_rxd))
        self.add_rxd_calcium(buffer_list)
        for sec in self.sections:
            if sec not in self.sections_rxd:
                if h.ismembrane("ca_ion", sec=sec):
                    sec.insert("cad")
                    sec.cao = self.params["Ca_Ext"]
        return

    def outermost_shell(diam, shell_width):
        return 2*shell_width/diam
    
    def add_shells(self, section, start_from, first_shell):
        name = section.name()
        which_dend = name.replace("[", "").replace("]", "")
        i = start_from
        new_factor = first_shell
        last_shell = False
        while True:
            if i == 0:
                self.membrane[name] = rxd.Region(section,
                                                 name='%s_membrane' %
                                                 which_dend,
                                                 geometry=rxd.membrane())
                self.factors[name] = []
                self.shells[name] = []
                
            inner = 1-(sum(self.factors[name])+new_factor)
            outer = 1-sum(self.factors[name])
            self.shells[name].append(rxd.Region(section, nrn_region='i',
                                                geometry=rxd.Shell(inner,
                                                                   outer),
                                                name="%s_Shell_%d" %
                                                (which_dend, i)))
            if i == 1:
                self.borders[name] = []
            if i > 0:
                self.borders[name].append(rxd.Region(section,
                                                     name='%s_Border_%d'
                                                     % (which_dend, i-1),
                                                     geometry=rxd.ScalableBorder(diam_scale=outer)))

            self.factors[name].append(new_factor)
            if last_shell:
                break
            if sum(self.factors[name]) + 2*new_factor >= self.fc:
                last_shell = True
                new_factor = (self.fc - sum(self.factors[name]))
            else:
                new_factor = 2*new_factor 
            i = i+1               
        
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
        membrane_shell_width = self.params["membrane_shell_width"]
        for sec in dendrites:
            factor = 2*membrane_shell_width/sec.diam
            if sec not in self.where_spines: 
                self.add_shells(sec, 0, factor)
            else:
                #outermost shell is the outermost shell of the dendrite
                # and the spine/spines
                sec_name = sec.name()
                which_dend = sec_name.replace("[", "").replace("]", "")
                new_secs = self.cell_filter(sec_name, tolist=True)
                all_geom[sec_name] = []
                for new_sec in new_secs:
                    if "head" in new_sec.name():
                        all_geom[sec_name].append(rxd.inside)
                    else:
                        factor = 2*membrane_shell_width/new_sec.diam
                        all_geom[sec_name].append(rxd.Shell(1-factor, 1))
                self.shells[sec_name] = [rxd.Region(new_secs,
                                                    nrn_region="i",
                                                    name="%s_Shell_0"
                                                    % which_dend,
                                                    geometry=rxd.MultipleGeometry(secs=new_secs, geos=all_geom[sec_name]))]
                self.factors[sec_name] = [factor]
                self.membrane[sec_name] = rxd.Region(new_secs,
                                                     name='%s_membrane'
                                                     % which_dend,
                                                     geometry=rxd.membrane())
                self.add_shells(sec, 1, 2*factor)
        if self.add_ER:
            self._add_ER_and_membrane()
        self._make_object_lists()

    def _add_ER_and_membrane(self):
        for sec in self.sections_rxd:
            sec_name = sec.name()
            which_dend = sec_name.replace("[", "").replace("]", "")
            l_s_d = 1 - sum(self.factors[sec_name]) #last shell diameter
            self.ER[sec_name] = rxd.Region(sec, nrn_region='i',
                                           geometry=rxd.Shell(0,
                                                              l_s_d),
                                           name="%s_ER" % which_dend)

            self.cyt_er_membrane[sec_name] = rxd.Region(sec,
                                                        name='%s_ER_membrane'
                                                        % which_dend,
                                                        geometry=rxd.ScalableBorder(
                                                            diam_scale=l_s_d))



    def _make_object_lists(self):
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
        caDiff = self.params["caDiff"]
        self.ca = rxd.Species(regions, d=caDiff, name='ca', charge=2,
                              initial=lambda nd: self.get_ca_init(nd.region.name),
                              atolscale=1e-8)
        
        if self.add_ER:
            ip3Diff = self.params["ip3Diff"]
            ip3_init = self.params["ip3_init"]
            calr_tot = self.params["calr_tot"]
            calr_bound = self.params["calr_bound"]
            self.ip3 = rxd.Species(self.shell_list, d=ip3Diff, initial=ip3_init,
                                   atolscale=1e-8)
            self.ca_ER = rxd.Species(self.ER_regions, d=caDiff, name='ca_er', charge=2,
                                     initial=self.params["ca_init_ER"],
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
            head_factor = self.params["head_factor"]
            return value[sec_name]*area*head_factor
        factor = self.factors[sec_name][0]
        vol_c = self.params["vol_c"]
        shell_vol = pi/4*node.sec.L*node.sec.diam**2*(1-(1-factor)**2)/vol_c
        return value[sec_name]*area*shell_vol

    def add_surface_pumps(self):
        gncx_spine = OrderedDict()
        gpmca_spine = OrderedDict()
        gncx_dend = OrderedDict()
        gpmca_dend = OrderedDict()
        gncx_spine_total = self.params["gncx_spine_total"]
        gncx_spine_bound = self.params["gncx_spine_bound"]
        gpmca_spine_total = self.params["gpmca_spine_total"]
        gpmca_spine_bound = self.params["gpmca_spine_bound"]
        gncx_dend_total = self.params["gncx_dend_total"]
        gncx_dend_bound = self.params["gncx_dend_bound"]
        gpmca_dend_total = self.params["gpmca_dend_total"]
        gpmca_dend_bound = self.params["gpmca_dend_bound"]
        for key in gncx_spine_total.keys():
            gncx_spine[key] = gncx_spine_total[key] - gncx_spine_bound[key]
            gpmca_spine[key] = gpmca_spine_total[key] - gpmca_spine_bound[key]
            gncx_dend[key] = gncx_dend_total[key] - gncx_dend_bound[key]
            gpmca_dend[key] = gpmca_dend_total[key] - gpmca_dend_bound[key]
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
    
        self.gserca = rxd.Parameter(self.cyt_er_membrane_list,
                                    initial=self.params["gserca"])
        #Serca from the sig paths model
        kcat_Serca = self.params["kcat_Serca"]
        Kserca = self.params["Kserca"]
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
        self.gleak = rxd.Parameter(self.cyt_er_membrane_list,
                                   initial=self.params["gleak"])
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
        self.gip3r = rxd.Parameter(self.cyt_er_membrane_list,
                                   initial=self.params["gip3r"])
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
            self.m.append(self.ip3[ca_region]/(self.ip3[ca_region] +
                                               self.params["Kip3"]))
            self.n.append(self.ca[ca_region]/(self.ca[ca_region]  +
                                              self.params["Kact"]))
            self.h_inf.append(self.params["k_inh"]/(self.ca[ca_region] +
                                                    self.params["k_inh"]))
            minf = self.m[i] * self.n[i]
        

            self.ip3rg.append(rxd.Rate(self.h_gate[membrane], (self.h_inf[i] - self.h_gate[membrane]) /
                                       self.params["ip3rtau"]))
            k = self.gip3r[membrane] * (minf * self.h_gate[membrane]) ** 3

            self.ip3r.append(rxd.MultiCompartmentReaction(self.ca_ER[er_region],
                                                         self.ca[ca_region],
                                                         k, k, membrane=membrane))
        
        # IP3 receptor gating



    def add_calmodulin(self):
        camDiff = self.params["camDiff"]
        calmodulin_tot = self.params["calmodulin_tot"]
        camn = self.params["camn"]
        camc = self.params["camc"]
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
                calbDiff = self.params["calbDiff"]
                calbindin_tot = self.params["calbindin_tot"]
                calbca = self.params["calbca"]
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
                tot_magnesium_green_BS = self.params["tot_magnesium_green_BS"]
                magnesium_green_bound = self.params["magnesium_green_bound"]
                mggreenDiff = self.params["mggreenDiff"]
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
        fixed_buffer_tot = self.params["fixed_buffer_tot"]
        fixed_buffer_ca = self.params["fixed_buffer_ca"]
        self.fixed = rxd.Species(self.shell_list,
                                 initial=fixed_buffer_tot-fixed_buffer_ca,
                                 name='FixedBuffer',
                                 charge=0, atolscale=1e-8)
        self.fixedca = rxd.Species(self.shell_list,
                                   initial=fixed_buffer_ca,
                                   name='FixedBufferCa',
                                   charge=0, atolscale=1e-8)


        
    def add_buffer_reactions(self, buffer_list):
        kf_fixed_b = self.params["kf_fixed_b"]
        kb_fixed_b = self.params["kb_fixed_b"]
        kf_camn = self.params["kf_camn"]
        kb_camn = self.params["kb_camn"]
        kf_camc = self.params["kf_camc"]
        kb_camc = self.params["kb_camc"]
        kf_calbindin = self.params["kf_calbindin"]
        kb_calbindin = self.params["kb_calbindin"]
        kf_magnesium_green = self.params["kf_magnesium_green"]
        kb_magnesium_green = self.params["kb_magnesium_green"]
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
            kf_calr = self.params["kf_calr"]
            kb_calr = self.params["kb_calr"]
            rxn1 = rxd.Reaction(self.calr + self.ca_ER,
                                self.calrca, kf_calr, kb_calr)
            self.reactions.append(rxn1)


    def add_surface_pump_reactions(self):
        self.pmca_1_r = []
        self.pmca_2_r = []
        self.ncx_1_r = []
        self.ncx_2_r = []
        kf_pmca = self.params["kf_pmca"]
        kb_pmca = self.params["kb_pmca"]
        kcat_pmca = self.params["kcat_pmca"]
        kf_ncx  = self.params["kf_ncx"]
        kb_ncx = self.params["kb_ncx"]
        kcat_ncx = self.params["kcat_ncx"]
       
        for key in self.shells.keys():
            membrane_shell = self.shells[key][0]
            membrane = self.membrane[key]
            self.pmca_1_r.append(rxd.MultiCompartmentReaction(self.ca[membrane_shell] +
                                                              self.pmca[membrane],
                                                              self.pmcaca[membrane],
                                                              kf_pmca*c_unit,
                                                              kb_pmca*c_unit,
                                                              membrane=membrane))
            self.pmca_2_r.append(rxd.MultiCompartmentReaction(self.pmcaca[membrane],
                                                              self.pmca[membrane] +
                                                              self.ca[self.ECS],
                                                              kcat_pmca*c_unit,
                                                              membrane=membrane))
            self.ncx_1_r.append(rxd.MultiCompartmentReaction(self.ca[membrane_shell] +
                                                             self.ncx[membrane],
                                                             self.ncxca[membrane],
                                                             kf_ncx*c_unit,
                                                             kb_ncx*c_unit,
                                                             membrane=membrane))
            self.ncx_2_r.append(rxd.MultiCompartmentReaction(self.ncxca[membrane],
                                                             self.ncx[membrane] +
                                                             self.ca[self.ECS],
                                                             kcat_ncx*c_unit,
                                                             membrane=membrane))
            
    def _add_diffusion(self):
        self.diffusions = []
        caDiff = self.params["caDiff"]
        diffusions = self.params["diffusions"]
        for sec_name in self.shells.keys():
            for i, shell in enumerate(self.shells[sec_name][:-1]):
                dname = sec_name.replace("[", "").replace("]", "")
                f = self.factors[sec_name][i]
                dr = rxd.Parameter(self.borders[sec_name][i],
                                   name=dname, 
                                   value=lambda nd:
                                   nd.segment.diam/2/f)
                rxn = rxd.MultiCompartmentReaction(self.ca[shell],
                                                   self.ca[self.shells[sec_name][i+1]], 
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.borders[sec_name][i])
                self.diffusions.append(rxn)
                for buf_name in self.buffers.keys():
                    dif_const = diffusions[buf_name]
                    for buffs in self.buffers[buf_name]:
                        rxn = rxd.MultiCompartmentReaction(buffs[shell],
                                                           buffs[self.shells[sec_name][i+1]], 
                                                           c_unit*dif_const/dr, 
                                                           c_unit*dif_const/dr,
                                                           border=self.borders[sec_name][i])
                        self.diffusions.append(rxn)

    def add_rxd_calcium(self, buffer_names):
        if self.sections_rxd:
            for sec in self.sections_rxd:
                 if "head" in sec.name():
                     sec.nseg = self.params["n_seg"]
            self._add_rxd_regions()
            self._add_species(buffer_names)
            self._add_diffusion()
            if self.add_ER:
                 self.add_serca()
                 self.add_ip3r()
                 self.add_leak()
            self.add_buffer_reactions(buffer_names)
            self.add_surface_pump_reactions()
     

