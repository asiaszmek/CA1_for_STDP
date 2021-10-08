import os
import pkg_resources
import json
import collections
import numpy

from neuron import h
from hippounit import tests

import CA1
from TestLoader import ModelLoader
base_directory = 'validation_results/'


if __name__ == "__main__":
    add_ER = False
    where_spines = []
    where_ca = ["soma", "apical"]
    mods_path = os.path.join(CA1.path, "Mods")
    my_model = ModelLoader(CA1.CA1_PC, mods_path,
                           {"add_ER": add_ER, "where_ca": where_ca,
                            "where_spines": where_spines}, "CA1_Ca_no_ER")

    my_model.v_init = -70
    my_model.celsius = 34
    with open('target_features/depol_block_target_data.json') as f:
        observation = json.load(f,
                                object_pairs_hook=collections.OrderedDict)
    
    test = tests.DepolarizationBlockTest(observation=observation,
                                         force_run=True,
                                         show_plot=True,
                                         save_all=False,
                                         base_directory=base_directory)
    score = test.judge(my_model)

    print(score.summary)

