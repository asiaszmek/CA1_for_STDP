from neuron import h

def record_specie_vec(specie, regions, record_dt, spine=False, special_name=""):
    '''use vec.record to record cyt calcium concentration. Will return list containing                                                                    
    vectors, one vector for each node'''
    Recordl = []
    if not isinstance(regions, list):
        regions = [regions]
    for region in regions:
        for mynode in specie[region].nodes:
            if spine:
                if "head" in mynode.sec.name():
                    cytcavec = h.Vector()
                    cytcavec.record(mynode._ref_concentration, record_dt)
                    Recordl.append(cytcavec)
            else:
                if "head" not in mynode.sec.name():
                    if not special_name:
                        #print(mynode.segment, mynode.x, mynode.sec.trueparentseg())
                        cytcavec = h.Vector()
                        cytcavec.record(mynode._ref_concentration, record_dt)
                        Recordl.append(cytcavec)
                    elif special_name in mynode.sec.name():
                        cytcavec = h.Vector()
                        cytcavec.record(mynode._ref_concentration, record_dt)
                        Recordl.append(cytcavec)

    return Recordl
