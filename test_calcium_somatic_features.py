import time
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
    test.n = 10
    score = test.judge(my_model)
    print(score.summary)

    # with open('target_features/feat_rat_CA1_JMakara_more_features.json') as f:
    #     observation = json.load(f,
    #                             object_pairs_hook=collections.OrderedDict)
    # stim_file = pkg_resources.resource_filename("hippounit",
    #                                             "tests/stimuli/somafeat_stim/stim_rat_CA1_PC_JMakara.json")
    # with open(stim_file, 'r') as f:
    #     config = json.load(f, object_pairs_hook=collections.OrderedDict)
    
    # # Instantiate test class   
    # test = tests.SomaticFeaturesTest(observation=observation,
    #                                  config=config, force_run=False,
    #                                  show_plot=True,
    #                                  save_all = True,
    #                                  base_directory=base_directory)

    # test.specify_data_set = 'JMakara_data'
    # test.npool = 10

    # score = test.judge(my_model)
    # #Summarize and print the score achieved by the model on the test using SciUnit's summarize function
    # score.summarize()
 
    

    ## Synaptic/spine part is not ready yet
    # with open('target_features/oblique_target_data.json') as f:
    #     observation = json.load(f, object_pairs_hook=collections.OrderedDict)
    # test = tests.ObliqueIntegrationTest(observation=observation,
    #                                     save_all=False, force_run_synapse=True,
    #                                     force_run_bin_search=False,
    #                                     show_plot=True,
    #                                     base_directory=base_directory)

    # test.npool = 10
    # score = test.judge(my_model)
    # print(score.summary)

    # with open("target_features/feat_PSP_attenuation_target_data.json", 'r') as f:
    #     observation = json.load(f, object_pairs_hook=collections.OrderedDict)
    
    # stim_file = pkg_resources.resource_filename("hippounit",
    #                                             "tests/stimuli/PSP_attenuation_stim/stim_PSP_attenuation_test.json")

    # with open(stim_file, 'r') as f:
    #     config = json.load(f, object_pairs_hook=collections.OrderedDict)

    # # Instantiate test class 
    # test = tests.PSPAttenuationTest(config=config, observation=observation,
    #                                 num_of_dend_locations = 15, force_run=True,
    #                                 show_plot=True, save_all=False,
    #                                 base_directory=base_directory)
            
    # # Number of parallel processes
    # test.npool = 10

    # score = test.judge(my_model)
    # #Summarize and print the score achieved by the model on the test using SciUnit's summarize function
    # print(score.summary)

