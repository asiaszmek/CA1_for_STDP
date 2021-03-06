import os
from collections import OrderedDict
import numpy as np
from neuron import rxd, h
from .ca_params import *
from .basal_cell import CA1_PC_basal
inc_fac = 2


# Nomenclature and values adapted from Harris KM, Jensen FE, Tsao BE.
# J Neurosci 1992
#absolute quantities in mM
#reaction rates in uM/ms or 1/ms
my_loc = os.path.dirname(os.path.abspath(__file__))
params_file = os.path.join(my_loc, "ca_params.py")

c_unit = 6.0221409e5


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
    
    ampas = []
    nmdas = []
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

    def find_sec(self, name):
        for sec in self.sections:
            if sec.name() == name:
                return sec

    def __init__(self, model=None, params=params_file,
                 where_ca=["lm_medium2"],
                 spine_number=1,
                 where_spines=["lm_medium2"],
                 add_ER=False, buffer_list=["Calmodulin", "Calbindin", "Fixed"],
                 pump_list=["ncx", "pmca"], receptor_list=["AMPA", "NMDA"],
                 spine_pos={}, recompile=True, add_to_h=True, v_init=-65, celsius=34):

        h.load_file("stdrun.hoc")
        h.CVode()
        h.CVode().active(True)

        if model is None:
            model = CA1_PC_basal(recompile=recompile,
                                 add_to_h=add_to_h)
        try:
            self.soma = model.soma[0]
        except TypeError:
            self.soma = model.soma
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
        self.pumps = []
        self.leak = []
        self.g = {}
        self.add_ER = add_ER
        if self.add_ER:
            self.ER = OrderedDict()
            self.cyt_er_membrane = OrderedDict()
            self.fc = 0.8 # fraction of cytoplasmic volume
        else:
            self.fc = 1
        self.where_ca = where_ca
        secs = []
        for sec in where_ca:
            secs.extend(self.cell_filter(sec, tolist=True))
        for sec in secs:
            channels = sec.psection()["density_mechs"].keys()
            for key in channels:
                if not key.startswith("ca"):
                    continue
                for seg in sec:
                    from_mech = getattr(seg, key)
                    value = getattr(from_mech, "gbar")
                    setattr(from_mech, "gbar",
                            self.params["ca_factor"]*value)
        self.comp_factors = {}
        self.reactions = []
        self.where_spines = []
        self.sections_rxd = []
        self.add_spines(where_spines, spine_number, spine_pos)
        self.pump_list = pump_list
        self.add_calcium(where_rxd=where_ca, buffer_list=buffer_list,
                         pump_list=pump_list)
        for head in self.heads:
            if "AMPA" in receptor_list:
                self.add_synapse_ampa(head(0.9),
                                      self.params["gAMPA"])
            if "NMDA" in receptor_list:
                self.add_synapse_nmda(head(0.9),
                                      self.params["gNMDA"],
                                      self.params["Ca_per"],
                                      is_ca=True)
        h.celsius = celsius
        self.v_init = v_init

    def add_spines(self, dends, spine_no, spine_pos={}):
        """
        Place either spine_no spines on every dend in dends, or spines provided by spine_no
        (a dictionary of spine position on chosen dend), if dend not in spine_pos, spine_no
        dendrites will be placed.
        """

        if dends:
            new_locs = []
            if not isinstance(dends, list):
                dends = [dends]
            for loc in dends:
                new_locs.extend(self.cell_filter(loc, tolist=True))
                self.where_spines.extend(new_locs)
            for spine_loc in new_locs:
                if spine_loc.name() in spine_pos:
                    positions = spine_pos[spine_loc.name()]
                    for a_spine, pos in enumerate(positions):
                        head = self.add_head(a_spine, where=spine_loc, position=pos)
                        self.heads.append(head)
                else:
                    for a_spine in range(spine_no):
                        pos = 1/(spine_no+1)*(a_spine+1)
                        head = self.add_head(a_spine, where=spine_loc, position=pos)
                        self.heads.append(head)

                self.compensate_for_spines(spine_loc)
        
    def get_ca_init(self, reg):
        if "ECS" in reg:
            return self.params["ca_ECS"]
        return self.params["ca_init"]

    def add_synapse_ampa(self, dend, gmax):
        syn = h.AMPADJ(dend)
        syn.gmax = gmax
        self.ampas.append(syn)
        return syn

    def add_synapse_nmda(self, dend, gmax, ca_per,
                         is_ca=True):
        if not is_ca:
            syn = h.NR2A(dend)
            syn.gmax = gmax
            self.nmdas.append(syn)
            return syn

        syn = h.NR2A_CA(dend)
        syn.fCa = ca_per
        syn.gmax = gmax
        self.nmdas.append(syn)
        return syn

    def add_spine_mechs(self, spine, segment, spine_mechanisms,
                        spine_params_dict, head=True):
        for mech_name in spine_mechanisms:
            spine.insert(mech_name)
        head_factor = self.params["head_factor"]
        for i, seg in enumerate(spine):
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
                    if mech_name == "pas":
                        value_e = spine_params_dict[mech_name]["e"]
                        setattr(to_mech, "e", value_)
                        value_e = spine_params_dict[mech_name]["g"]
                        setattr(to_mech, "g", value_)
                    else:
                        for name in spine_params_dict[mech_name]:
                            if i > spine.nseg/2 - 1 or mech_name == "cad":
                                if mech_name in head_factor:
                                    value = spine_params_dict[mech_name][name]*head_factor
                                else:
                                    value = spine_params_dict[mech_name][name]
                            else:
                                value = 0
                            setattr(to_mech, name, value)
                else:
                    if mech_name == "pas":
                        self.set_val_from_segment(seg, segment,
                                                  mech_name, "e")
                        self.set_val_from_segment(seg, segment,
                                                  mech_name, "g")
                    else:
                        if i > spine.nseg/2 - 1:
                            self.set_val_from_segment(seg, segment,
                                                  mech_name, "gbar")
                        else:
                            to_mech = getattr(seg, mech_name)
                            setattr(to_mech, "gbar", 0)
                
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
                 where="apical_dendrite[10]",
                 position=0.5, head_mechanisms=["pas", "kdr", "nax", "cal12",
                                                "cal13", "can",
                                                "car", "cav33", "cav32", 'hd',
                                                'kad', 'kdr'],
                 head_params_dict={"Ra":120000}):
    
        if isinstance(where, str): 
            dend = self.cell_filter(where)
        else:
            dend = where
            where = dend.name()
        segment = dend(position)
        head = h.Section(name="%s_%s_%d" % (dend.name(), 'head', number))
        if where in self.where_ca:
            head.nseg = 5
        else:
            head.nseg = 3

        head.L = (SPINE_DIMENSIONS[spine_type]["head_len"] +
                  SPINE_DIMENSIONS[spine_type]["neck_len"])
                  
        for i, seg in enumerate(head):
            if i > head.nseg/2 - 1 :
                seg.diam = SPINE_DIMENSIONS[spine_type]["head_diam"]
            else:
                seg.diam = SPINE_DIMENSIONS[spine_type]["neck_diam"]
        head.connect(segment)
        
        self.sections.append(head)
        if segment not in self.spine_dict:
            self.spine_dict[segment] = []
        self.spine_dict[segment].append(head)
        self.add_spine_mechs(head, segment, head_mechanisms,
                             head_params_dict)
        return head

    #Add a number of spines to a dendritic branch.
    #lower conductances, scale resistance and capacitance

    def compensate_for_spines(self, dend):
        for segment in dend:
            if segment in self.spine_dict:
                spines = self.spine_dict[segment]
                section = segment.sec
                seg_surf = segment.area()
                mech_dict = section.psection()["density_mechs"]
                ca_channels = []
                for mech in mech_dict:
                    if 'gbar' in mech_dict[mech].keys():  # is a channel
                        tot_spine_cond = 0
                        for spine in spines:
                            for i, spine_seg in enumerate(spine):
                                if mech in spine.psection()["density_mechs"]:
                                    cond = spine.psection()["density_mechs"][mech]['gbar'][i]
                                    tot_spine_cond += spine_seg.area() * cond
                        seg_mech = getattr(segment, mech)
                        seg_cond = getattr(seg_mech, "gbar")
                        new_cond = (seg_cond*seg_surf - tot_spine_cond)/seg_surf
                        setattr(seg_mech, "gbar", new_cond)
                        if mech.startswith("ca") or mech.startswith("Ca"):
                            if seg_cond:
                                ca_channels.append(new_cond/seg_cond)
                self.comp_factors[segment] = np.array(ca_channels).mean()
                
                spine_g = 0
                spine_cm = 0
                for spine in spines:
                    for spine_seg in spine:
                        spine_g += spine_seg.area()*spine_seg.g_pas
                        spine_cm += spine_seg.area()*spine_seg.cm
                new_g = (segment.g_pas*seg_surf - spine_g)/seg_surf
                segment.g_pas = new_g
                new_cm = (segment.cm*seg_surf - spine_cm)/seg_surf
                segment.cm = new_cm

    def add_calcium(self, where_rxd=[], buffer_list=["Calmodulin", "Calbindin"],
                    pump_list=["pmca", "ncx"]):

        for name in where_rxd:
            if name == "apical":
                self.sections_rxd.extend(self.apical)
            elif name == "basal":
                self.sections_rxd.extend(self.basal)
            elif name == "soma":
                self.sections_rxd.extend(self.somatic)
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
        self.add_rxd_calcium(buffer_list, pump_list)
        for sec in self.sections:
            if h.ismembrane("ca_ion", sec=sec):
                sec.insert("bk")
                sec.gbar_bk = self.params["gbar_bk"][sec.name()]
                sec.insert("sk")
                sec.gbar_sk = self.params["gbar_sk"][sec.name()]
                if sec not in self.sections_rxd:
                    #print("Adding simple Ca dynamics to %s" % sec.name())
                    sec.insert("cad")
                    if len(buffer_list) > 2:
                        sec.Buffer_cad = 25
                        sec.cainf_cad = self.params["ca_init"]
                        #an additional buffer will change Ca dynamics
                    
        h.cao0_ca_ion = self.params["Ca_Ext"]
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
                
            inner = 1 - (sum(self.factors[name])+new_factor)
            if inner < 0:
                inner = 0
            outer = 1 - sum(self.factors[name])
            if i == 0 and inner == 0:
                self.shells[name].append(rxd.Region(section, nrn_region='i',
                                                    name="%s_Shell_%d" %
                                                    (which_dend, i)))
            else:
                if i:
                    self.shells[name].append(rxd.Region(section,
                                                        geometry=rxd.Shell(inner,
                                                                           outer),
                                                        name="%s_Shell_%d" %
                                                        (which_dend, i)))
                else:
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
            if not inner:
                break
            if sum(self.factors[name]) + inc_fac*new_factor >= self.fc:
                last_shell = True
                new_factor = (self.fc - sum(self.factors[name]))
            else:
                new_factor = inc_fac*new_factor 
            i = i+1               

    def _add_rxd_regions(self):
        self.membrane = OrderedDict()
        self.shells = OrderedDict()
        self.borders = OrderedDict()

        self.factors = OrderedDict()  #  membrane_shell_width/max_radius
        secs_spines = OrderedDict()
        sec_list = []
        for sec in self.sections_rxd:
            if "head" in sec.name():
                continue
            sec_list.append(sec)
        dendrites = sec_list[:]
        all_geom = OrderedDict()
        membrane_shell_width = self.params["membrane_shell_width"]

        for sec in dendrites:
            print("Add rxd Ca to %s" % sec.name())
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
                        if factor > 1:
                            factor = 1
                        all_geom[sec_name].append(rxd.Shell(1-factor, 1))
                geom = rxd.MultipleGeometry(secs=new_secs,
                                            geos=all_geom[sec_name])
                self.shells[sec_name] = [rxd.Region(new_secs,
                                                    nrn_region="i",
                                                    name="%s_Shell_0"
                                                    % which_dend,
                                                    geometry=geom)]
                self.factors[sec_name] = [factor]
                self.membrane[sec_name] = rxd.Region(new_secs,
                                                     name='%s_membrane'
                                                     % which_dend,
                                                     geometry=rxd.membrane())
                if factor == 1:
                    break
                self.add_shells(sec, 1, 2*factor)
        if self.add_ER:
            self._add_ER_and_membrane()
        self._make_object_lists()

    def _add_ER_and_membrane(self):
        for sec in self.sections_rxd:
            if "head" in sec.name():
                continue
            sec_name = sec.name()
            which_dend = sec_name.replace("[", "").replace("]", "")
            l_s_d = 1 - sum(self.factors[sec_name]) #last shell diameter
            self.ER[sec_name] = rxd.Region(sec,
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
        for key in self.shells.keys():
            self.membrane_list.append(self.shells[key][0])
        if self.add_ER:
            self.ER_regions = []
            self.cyt_er_membrane_list = []
            self.shells_to_ER = []
            for sec_name in self.ER.keys():
                self.ER_regions.append(self.ER[sec_name])
                self.cyt_er_membrane_list.append(self.cyt_er_membrane[sec_name])
                self.shells_to_ER.append(self.shells[sec_name][-1])

    def _add_species(self, buffer_list, pump_list):
        regions = self.shell_list 
        caDiff = self.params["caDiff"]
        self.ca = rxd.Species(regions, d=caDiff, name='ca', charge=2,
                              initial=lambda nd:
                              self.get_ca_init(nd.region.name),
                              atolscale=1e-9)
        
        if self.add_ER:
            ip3Diff = self.params["ip3Diff"]
            ip3_init = self.params["ip3_init"]
            calr_tot = self.params["calr_tot"]
            calr_bound = self.params["calr_bound"]
            self.ip3 = rxd.Species(self.shell_list, d=ip3Diff, initial=ip3_init,
                                   atolscale=1e-9)
            self.ca_ER = rxd.Species(self.ER_regions, d=caDiff, name='caer',
                                     charge=2,
                                     initial=self.params["ca_init_ER"],
                                     atolscale=1e-9)
            self.calr = rxd.Species(self.ER_regions,
                                    initial=calr_tot-calr_bound)
            self.calrca = rxd.Species(self.ER_regions, initial=calr_bound)
        self.add_buffers(buffer_list)

    def ncx_val(self, node):
        if "head" in node.sec.name():
            return self.params["gncx_spine"]*self.params["kcat_ncx"]
        if "soma" in node.sec.name():
            return self.params["gncx_soma"]*self.params["kcat_ncx"]
        # h.distance(sec=self.soma)
        # dist = h.distance(node.sec(0.5), sec=node.sec)
        # if dist > 100:
        #     return self.params["gncx"]*self.params["kcat_ncx"]/5
        return self.params["gncx"]*self.params["kcat_ncx"]
        

    def pmca_val(self, node):
        if "head" in node.sec.name():
            return self.params["gpmca_spine"]*self.params["kcat_pmca"]
        if "soma" in node.sec.name():
            self.params["gpmca_soma"]*self.params["kcat_pmca"]
        # h.distance(sec=self.soma)
        # dist = h.distance(node.sec(0.5), sec=node.sec)
        # if dist > 100:
        #     return self.params["gpmca"]*self.params["kcat_pmca"]/5
        return self.params["gpmca"]*self.params["kcat_pmca"]

    def add_pump(self, name):
        memb_flux = False
        if name == "Serca":
            membranes = self.cyt_er_membrane
            n = 2
            membrane_list = [membranes[key] for key in membranes.keys()]
            self.g[name] = rxd.Parameter(membrane_list,
                                         value=lambda nd:
                                         self.params["g%s" % name])
        # tu shell a nie membrane
        elif name == "pmca":
            membranes = self.membrane
            memb_flux = True
            n = 1
            membrane_list = self.membrane_list
            #[membranes[key] for key in membranes.keys()]
            self.g[name] = rxd.Parameter(membrane_list,
                                         value=lambda nd: self.pmca_val(nd))
        elif name == "ncx":
            membranes = self.membrane
            memb_flux = True
            n = 1
            membrane_list = self.membrane_list
            # [membranes[key] for key in membranes.keys()]
            self.g[name] = rxd.Parameter(membrane_list,
                                         value=lambda nd: self.ncx_val(nd))
        else:
            membranes = self.membrane
            membrane_list = [membranes[key] for key in membranes.keys()]
            n = 1

        kcat_pump = self.params["kcat_%s" % name]
        Km_pump = self.params["Km_%s" % name]
        for key in self.shells.keys():
            if "head" in key:
                continue
            if name == "Serca":
                outside = self.ca_ER[self.ER[key]]
                inside = self.ca[self.shells[key][-1]]
            else:
                outside = None
                inside = self.ca[self.shells[key][0]]
            membrane = membranes[key]
            if outside is not None:
                rate = self.g[name]*inside**n*kcat_pump/(Km_pump**n+inside**n)
                pump = rxd.MultiCompartmentReaction(inside > outside,
                                                    rate,
                                                    membrane=membrane,
                                                    custom_dynamics=True,
                                                    membrane_flux=memb_flux)
                self.pumps.append(pump)
            else:
                
                # Keener and Sneyd, mathematical physiology, page 32
                basal = self.params["ca_init"]
                extrusion =-self.g[name]/(Km_pump/inside + 1)
                leak = self.g[name]/(Km_pump/basal + 1)
                pump = rxd.Rate(inside, extrusion + leak,
                                regions=[self.shells[key][0]])
                self.pumps.append(pump)
    
    def add_leak(self):
        # leak channel: bidirectional ca flow btwn cyt <> ER
        self.leak = []
        for key in self.ER.keys():
            er_region =	self.ER[key]
            ca_region =	self.shells[key][-1]
            membrane = self.cyt_er_membrane[key]
            gleak = rxd.Parameter(ca_region, value=lambda nd:
                                  self.params["gleak"])
            self.leak.append(rxd.MultiCompartmentReaction(self.ca_ER[er_region],
                                                          self.ca[ca_region],
                                                          gleak,
                                                          gleak,
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
        

            self.ip3rg.append(rxd.Rate(self.h_gate[membrane],
                                       (self.h_inf[i] - self.h_gate[membrane])
                                       /self.params["ip3rtau"]))
            self.k = self.gip3r[membrane] * (minf * self.h_gate[membrane]) ** 3

            self.ip3r.append(rxd.MultiCompartmentReaction(self.ca_ER[er_region],
                                                          self.ca[ca_region],
                                                          self.k, self.k,
                                                          membrane=membrane))
        
            # IP3 degradation
            for shell in self.shells[key]:
                self.reactions.append(rxd.Rate(self.ip3,
                                               (ip3_init-self.ip3[shell])
                                               /ip3degTau,
                                               regions=shell,
                                               membrane_flux=False))


    def add_calmodulin(self):
        camDiff = self.params["camDiff"]
        calmodulin_tot = self.params["calmodulin_tot"]
        camn = self.params["camn"]
        camc = self.params["camc"]
        self.cam = rxd.Species(self.shell_list, d=camDiff,
                               initial=calmodulin_tot-camn-camc,
                               name='Calmodulin', charge=0, atolscale=1e-9)
        self.camn = rxd.Species(self.shell_list, d=camDiff, initial=camn,
                                name='CaMN',
                                charge=0, atolscale=1e-9)
        self.camc = rxd.Species(self.shell_list, d=camDiff, initial=camc,
                                name='CaMC',
                                charge=0, atolscale=1e-9)
        self.buffers["Calmodulin"] = [self.cam, self.camn, self.camc]

    def add_buffers(self, buffer_names):
        self.buffers = OrderedDict()
        self.indicator = None
        self.indicator_ca = None
        for name in buffer_names:
            if name == "Calmodulin":
                self.add_calmodulin()
            elif name == "Calbindin":
                calbDiff = self.params["calbDiff"]
                calbindin_tot = self.params["calbindin_tot"]
                calbca = self.params["calbca"]
                self.calb = rxd.Species(self.shell_list, d=calbDiff,
                                        initial=calbindin_tot-calbca,
                                        name='Calbindin',
                                        charge=0, atolscale=1e-9)
                self.calb_ca = rxd.Species(self.shell_list, d=calbDiff,
                                          initial=calbca,
                                          name='CalbindinCa',
                                          charge=0, atolscale=1e-9)
                self.buffers["Calbindin"] = [self.calb, self.calb_ca]

            elif name == "Mg Green":
                tot_magnesium_green_BS = self.params["tot_magnesium_green_BS"]
                magnesium_green_bound = self.params["magnesium_green_bound"]
                mggreenDiff = self.params["mggreenDiff"]
                self.indicator = rxd.Species(self.shell_list,
                                             initial=tot_magnesium_green_BS -
                                             magnesium_green_bound, d=mggreenDiff,
                                             name='MgGreen',
                                             charge=0, atolscale=1e-9)
                self.indicator_ca = rxd.Species(self.shell_list,
                                             initial=magnesium_green_bound,
                                             name='MgGreenCa', d=mggreenDiff,
                                                charge=0, atolscale=1e-9)
                self.buffers["Mg Green"] = [self.indicator, self.indicator_ca]
            elif name == "Fluo3":
                tot_indicator = self.params["tot_fluo3"]
                indicator_bound = self.params["ca_init"]*tot_indicator*kf_fluo3/kb_fluo3
                indicatorDiff = self.params["fluo3Diff"]
                self.indicator = rxd.Species(self.shell_list,
                                             initial=tot_indicator -
                                             indicator_bound, d=indicatorDiff,
                                             name='Fluo3',
                                             charge=0, atolscale=1e-9)
                self.indicator_ca = rxd.Species(self.shell_list,
                                             initial=indicator_bound,
                                             name='Fluo3Ca', d=indicatorDiff,
                                            charge=0, atolscale=1e-9)
                self.buffers["Fluo3"] = [self.indicator, self.indicator_ca]
            elif name == "BF2":
                tot_indicator = self.params["tot_BF2"]
                indicator_bound = self.params["ca_init"]*tot_indicator*kf_BF2/kb_BF2
                indicatorDiff = self.params["BF2Diff"]
                self.indicator = rxd.Species(self.shell_list,
                                             initial=tot_indicator -
                                             indicator_bound, d=indicatorDiff,
                                             name='BF2',
                                             charge=0, atolscale=1e-9)
                self.indicator_ca = rxd.Species(self.shell_list,
                                             initial=indicator_bound,
                                             name='BF2Ca', d=indicatorDiff,
                                            charge=0, atolscale=1e-9)
                self.buffers["BF2"] = [self.indicator, self.indicator_ca]
            elif name == "OGB1":
                tot_indicator = self.params["tot_OGB1"]
                indicator_bound = 0.0388#self.params["ca_init"]*tot_indicator*kf_OGB1/kb_OGB1
                indicatorDiff = self.params["OGB1Diff"]
                self.indicator = rxd.Species(self.shell_list,
                                             initial=tot_indicator -
                                             indicator_bound, d=indicatorDiff,
                                             name='OGB1',
                                             charge=0, atolscale=1e-9)
                self.indicator_ca = rxd.Species(self.shell_list,
                                             initial=indicator_bound,
                                             name='OGB1Ca', d=indicatorDiff,
                                            charge=0, atolscale=1e-9)
                self.buffers["OGB1"] = [self.indicator, self.indicator_ca]
            
            elif name == "Fixed":
                fixed_buffer_tot = self.params["fixed_buffer_tot"]
                fixed_buffer_ca = self.params["fixed_buffer_ca"]
                self.fixed = rxd.Species(self.shell_list,
                                         initial=fixed_buffer_tot
                                         -fixed_buffer_ca,
                                         name='FixedBuffer',
                                         charge=0, atolscale=1e-9)
                self.fixed_ca = rxd.Species(self.shell_list,
                                           initial=fixed_buffer_ca,
                                           name='FixedBufferCa',
                                           charge=0, atolscale=1e-9)
                self.buffers["Fixed"] = [self.fixed, self.fixed_ca]

        
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
        if "Fixed" in buffer_list:
            fixed_rxn = rxd.Reaction(self.fixed + self.ca, self.fixed_ca,
                                     kf_fixed_b,
                                     kb_fixed_b)
            self.reactions.append(fixed_rxn)
        if "Calmodulin" in buffer_list:
            rn = rxd.Reaction(self.cam + self.ca, self.camn, kf_camn,
                              kb_camn)
            rc = rxd.Reaction(self.cam + self.ca, self.camc, kf_camc, kb_camc)
            self.reactions.extend([rn, rc])
            
        if "Calbindin" in buffer_list:
            calb_rxn = rxd.Reaction(self.calb + self.ca, self.calb_ca,
                                    kf_calbindin,
                                    kb_calbindin)
            self.reactions.append(calb_rxn)
        if "Mg Green" in buffer_list:
            mg_green_binding = rxd.Reaction(self.indicator + self.ca,
                                            self.indicator_ca,
                                            kf_magnesium_green,
                                            kb_magnesium_green)
            self.reactions.append(mg_green_binding)
        if "Fluo3" in buffer_list:
            fluo3_binding = rxd.Reaction(self.indicator + self.ca,
                                         self.indicator_ca,
                                         kf_fluo3,
                                         kb_fluo3)
            self.reactions.append(fluo3_binding)
        if "OGB1" in buffer_list:
            ogb1_binding = rxd.Reaction(self.indicator + self.ca,
                                         self.indicator_ca,
                                         kf_OGB1,
                                         kb_OGB1)
            self.reactions.append(ogb1_binding)
        if "BF2" in buffer_list:
            bf2_binding = rxd.Reaction(self.indicator + self.ca,
                                         self.indicator_ca,
                                         kf_BF2,
                                         kb_BF2)
            self.reactions.append(bf2_binding)

        if self.add_ER:
            kf_calr = self.params["kf_calr"]
            kb_calr = self.params["kb_calr"]
            rxn1 = rxd.Reaction(self.calr + self.ca_ER,
                                self.calrca, kf_calr, kb_calr)
            self.reactions.append(rxn1)


            
    def _add_diffusion(self):
        self.diffusions = []
        caDiff = self.params["caDiff"]
        diffusions = self.params["diffusions"]
        self.drs = []
        for sec_name in self.shells.keys():
            for i, shell in enumerate(self.shells[sec_name][:-1]):
                dname = sec_name.replace("[", "").replace("]", "")
                f = self.factors[sec_name][i]
                dr = rxd.Parameter(self.borders[sec_name][i],
                                   name=dname, 
                                   value=lambda nd:
                                   nd.segment.diam/2/f)
                self.drs.append(dr)
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
                if self.add_ER:
                    rxn = rxd.MultiCompartmentReaction(self.ip3[shell],
                                                       self.ip3[self.shells[sec_name][i+1]], 
                                                       c_unit*ip3Diff/dr, 
                                                       c_unit*ip3Diff/dr,
                                                       border=self.borders[sec_name][i])
                    self.diffusions.append(rxn)


    def add_rxd_calcium(self, buffer_names, pump_list):
        if self.sections_rxd:
            self._add_rxd_regions()
            self._add_species(buffer_names, pump_list)
            self._add_diffusion()
            if self.add_ER:
                 self.add_pump("Serca")
                 self.add_ip3r()
                 self.add_leak()
            self.add_buffer_reactions(buffer_names)
            #self.add_surface_pump_reactions(pump_list)
            for pump in self.pump_list:
                self.add_pump(pump)


    def make_a_run(self, tstop):
      h.CVode().re_init()
      h.finitialize(self.v_init)
      h.fcurrent()
      h.tstop = tstop
      h.run(tstop)
