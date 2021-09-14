import time
import os
from json2html import *
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
    where_spines = []#["radTprox1"]
    where_ca = ["soma", "apical"]
    #cell = CA1_PC(add_ER=add_ER, where_ca=where_ca, where_spines=where_spines)
    mods_path = os.path.join(CA1.path, "Mods")
    my_model = ModelLoader(CA1.CA1_PC, mods_path, {"add_ER": add_ER,
                                    "where_ca": where_ca,
                                    "where_spines": where_spines})
    #my_model.initialize()
    my_model.v_init = -70
    my_model.celsius = 34
    target_features_file = os.path.join("target_features",
                                        "feat_CA1_pyr_cACpyr_more_features.json")
    with open(target_features_file) as f:
        observation = json.load(f, object_pairs_hook=collections.OrderedDict)


    # Load stimuli file
    ttype = "CA1_pyr_cACpyr"

    stim_file = pkg_resources.resource_filename("hippounit",
                                                "tests/stimuli/somafeat_stim/stim_" + ttype + ".json")
    with open(stim_file, 'r') as f:
        config = json.load(f, object_pairs_hook=collections.OrderedDict)
    
    # Instantiate test class   
    test = tests.SomaticFeaturesTest(observation=observation, config=config,
                                     force_run=True, show_plot=True,
                                     save_all=False,
                                     base_directory=base_directory)
    test.specify_data_set = 'UCL_data'

    # Number of parallel processes
    test.n = 1
    score = test.judge(my_model)
    #score.summarize()
    # try:
    #     #Run the test 
    #     score = test.judge(my_model)
    #     #Summarize and print the score achieved by the model on the test using SciUnit's summarize function
    #     score.summarize()
    # except Exception as e:
    #     print('Model: ' + my_model.name + ' could not be run')
    #     print(e)
    #     pass
    target_features_file = os.path.join("target_features",
                                        "feat_backpropagating_AP_target_data.json")

    with open(target_features_file) as f:
        observation = json.load(f, object_pairs_hook=collections.OrderedDict)
    stim_file = pkg_resources.resource_filename("hippounit",
                                                "tests/stimuli/bAP_stim/stim_bAP_test.json")

    with open(stim_file, 'r') as f:
        config = json.load(f, object_pairs_hook=collections.OrderedDict)

    # Instantiate the test class
    test = tests.BackpropagatingAPTest(config=config, observation=observation,
                                       force_run=True, force_run_FindCurrentStim=False,
                                       show_plot=True, save_all=False,
                                       base_directory=base_directory)

    # Number of parallel processes
    test.npool = 10
    

    score = test.judge(my_model)
    #Summarize and print the score achieved by the model on the test using SciUnit's summarize function
    #score.summarize()
    
