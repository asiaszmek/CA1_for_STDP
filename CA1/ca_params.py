head_factor = 1
ip3rtau = 2000 #wagner
gip3r = 0.85e-3

gleak = 1.0625e-5
ip3_init = 0 # 50 nM? signalling pathways model

calr_tot = 86 #  calreticulin
calr_bound = 7.11945
kf_calr = 0.1 # 1/mM/ms
kb_calr = 0.2 #  1/ms

Ca_per = .1 # 10% at 2mM, 5% at 1mM, 20% at 4mM extracellular Ca  doi: 10.1073/pnas.90.24.11573

ca_init_ER = .180512
ca_ECS = 2.

caDiff = 0.174 
ip3Diff = 0.1
camDiff = 0.004
calbDiff = 0.009

AtoN_ratio = 2.1 # at 6-8 weeks  doi: 10.1113/jphysiol.2008.160929

gNMDA =2e-4
gAMPA = AtoN_ratio * gNMDA

Kip3 = 0.15e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
Kact = 0.8e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
k_inh = 2e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
Kserca = 0.1e-3 # Michaelis constant for SERCA pump
kcat_Serca = 7.5e-3 #ms
gserca = 0.5e-3 #Avramas hermissenda model
ip3r_gate_state = 0.8
calbindin_tot = 0.150 #(150 uM)
calmodulin_tot = 0.03 #(30 uM)
fixed_buffer_tot = 2 #(2 mM)
fixed_buffer_ca = 4.07965e-3
camn = 2.927e-05
camc = 1.965e-04
calbca = 2.189e-5


kf_calbindin = 2.8e-2
kb_calbindin = 0.0196

kf_camn = 2*7.7e2
kb_camn = 1.6e2
kf_camnn = 3.2e4
kb_camnn = 1.1e1
kf_camc = 2*8.4e1
kb_camc = 2.6
kf_camcc = 2.5e1
kb_camcc = 0.6e-3

kf_fixed_b = 400
kb_fixed_b = 20

kf_pmca = 50
kb_pmca = 7e-3
kcat_pmca = 3.5e-3
gpmca_dend_bound = {
    "CA1_PC_Tomko[0].radTprox1": 1.11e-6,
    "CA1_PC_Tomko[0].radTprox2": 1.11e-6,
    "CA1_PC_Tomko[0].radTmed1": 1.11e-6,
    "CA1_PC_Tomko[0].radTmed2": 1.11e-6,
    "CA1_PC_Tomko[0].radTdist1": 1.11e-6,
    "CA1_PC_Tomko[0].radTdist2": 1.11e-6,
    "CA1_PC_Tomko[0].rad_t1": 1.11e-6,
    "CA1_PC_Tomko[0].rad_t2": 1.11e-6,
    "CA1_PC_Tomko[0].rad_t3": 1.11e-6,
    "CA1_PC_Tomko[0].lm_thick2": 1.11e-6,
    "CA1_PC_Tomko[0].lm_medium2": 1.11e-6,
    "CA1_PC_Tomko[0].lm_thin2": 1.11e-6,
    "CA1_PC_Tomko[0].lm_thick1": 1.11e-6,
    "CA1_PC_Tomko[0].lm_medium1": 1.11e-6,
    "CA1_PC_Tomko[0].lm_thin1": 5.55e-7,
}
gpmca_spine_bound = {
    "CA1_PC_Tomko[0].radTprox1": 1.01975e-8,
    "CA1_PC_Tomko[0].radTprox2": 1.01975e-8,
    "CA1_PC_Tomko[0].radTmed1": 1.01975e-8,
    "CA1_PC_Tomko[0].radTmed2": 1.01975e-8,
    "CA1_PC_Tomko[0].radTdist1": 1.01975e-8,
    "CA1_PC_Tomko[0].radTdist2": 1.01975e-8,
    "CA1_PC_Tomko[0].rad_t1": 1.01975e-8,
    "CA1_PC_Tomko[0].rad_t2": 1.01975e-8,
    "CA1_PC_Tomko[0].rad_t3": 1.01975e-8,
    "CA1_PC_Tomko[0].lm_thick2": 1.01975e-8,
    "CA1_PC_Tomko[0].lm_medium2": 1.01975e-8,
    "CA1_PC_Tomko[0].lm_thin2": 1.01975e-8,
    "CA1_PC_Tomko[0].lm_thick1": 1.01975e-8,
    "CA1_PC_Tomko[0].lm_medium1": 1.01975e-8,
    "CA1_PC_Tomko[0].lm_thin1": 1.01975e-7,
}
gpmca_dend_total = {
    "CA1_PC_Tomko[0].radTprox1": 10.61e-6,
    "CA1_PC_Tomko[0].radTprox2": 10.61e-6,
    "CA1_PC_Tomko[0].radTmed1": 10.61e-6,
    "CA1_PC_Tomko[0].radTmed2": 10.61e-6,
    "CA1_PC_Tomko[0].radTdist1": 10.61e-6,
    "CA1_PC_Tomko[0].radTdist2": 10.61e-6,
    "CA1_PC_Tomko[0].rad_t1": 10.61e-6,
    "CA1_PC_Tomko[0].rad_t2": 10.61e-6,
    "CA1_PC_Tomko[0].rad_t3": 10.61e-6,
    "CA1_PC_Tomko[0].lm_thick2": 10.61e-6,
    "CA1_PC_Tomko[0].lm_medium2": 10.61e-6,
    "CA1_PC_Tomko[0].lm_thin2": 10.61e-6,
    "CA1_PC_Tomko[0].lm_thick1": 10.61e-6,
    "CA1_PC_Tomko[0].lm_medium1": 10.61e-6,
    "CA1_PC_Tomko[0].lm_thin1": 50.61e-7,
}
gpmca_spine_total = {
    "CA1_PC_Tomko[0].radTprox1": 9e-07,
    "CA1_PC_Tomko[0].radTprox2": 9e-07,
    "CA1_PC_Tomko[0].radTmed1": 9e-07,
    "CA1_PC_Tomko[0].radTmed2": 9e-07,
    "CA1_PC_Tomko[0].radTdist1": 9e-07,
    "CA1_PC_Tomko[0].radTdist2": 9e-07,
    "CA1_PC_Tomko[0].rad_t1": 9e-07,
    "CA1_PC_Tomko[0].rad_t2": 9e-07,
    "CA1_PC_Tomko[0].rad_t3": 9e-07,
    "CA1_PC_Tomko[0].lm_thick2": 9e-07,
    "CA1_PC_Tomko[0].lm_medium2": 9e-07,
    "CA1_PC_Tomko[0].lm_thin2": 9e-07,
    "CA1_PC_Tomko[0].lm_thick1": 9e-07,
    "CA1_PC_Tomko[0].lm_medium1": 9e-07,
    "CA1_PC_Tomko[0].lm_thin1": 9e-07
}

kf_ncx = 1.68e1  # 1/mM/ms
kb_ncx = 0.0112
kcat_ncx = 5.6e-3
gncx_dend_bound = {
    "CA1_PC_Tomko[0].radTprox1": 4.66e-8,
    "CA1_PC_Tomko[0].radTprox2": 4.66e-8,
    "CA1_PC_Tomko[0].radTmed1": 4.66e-8,
    "CA1_PC_Tomko[0].radTmed2": 4.66e-8,
    "CA1_PC_Tomko[0].radTdist1": 4.66e-8,
    "CA1_PC_Tomko[0].radTdist2": 4.66e-8,
    "CA1_PC_Tomko[0].rad_t1": 4.66e-8,
    "CA1_PC_Tomko[0].rad_t2": 4.66e-8,
    "CA1_PC_Tomko[0].rad_t3": 4.66e-8,
    "CA1_PC_Tomko[0].lm_thick2": 4.66e-8,
    "CA1_PC_Tomko[0].lm_medium2": 4.66e-8,
    "CA1_PC_Tomko[0].lm_thin2": 4.66e-8,
    "CA1_PC_Tomko[0].lm_thick1": 4.66e-8,
    "CA1_PC_Tomko[0].lm_medium1": 4.66e-8,
    "CA1_PC_Tomko[0].lm_thin1": 23e-9,
}
gncx_spine_bound = {
    "CA1_PC_Tomko[0].radTprox1": 0,
    "CA1_PC_Tomko[0].radTprox2": 0,
    "CA1_PC_Tomko[0].radTmed1": 0,
    "CA1_PC_Tomko[0].radTmed2": 0,
    "CA1_PC_Tomko[0].radTdist1": 0,
    "CA1_PC_Tomko[0].radTdist2": 0,
    "CA1_PC_Tomko[0].rad_t1": 0,
    "CA1_PC_Tomko[0].rad_t2": 0,
    "CA1_PC_Tomko[0].rad_t3": 0,
    "CA1_PC_Tomko[0].lm_thick2": 0,
    "CA1_PC_Tomko[0].lm_medium2": 0,
    "CA1_PC_Tomko[0].lm_thin2": 0,
    "CA1_PC_Tomko[0].lm_thick1": 0,
    "CA1_PC_Tomko[0].lm_medium1": 0,
    "CA1_PC_Tomko[0].lm_thin1": 0
}

gncx_dend_total = {
    "CA1_PC_Tomko[0].radTprox1": 5.01225e-7,
    "CA1_PC_Tomko[0].radTprox2": 5.01225e-7,
    "CA1_PC_Tomko[0].radTmed1": 5.01225e-7,
    "CA1_PC_Tomko[0].radTmed2": 5.01225e-7,
    "CA1_PC_Tomko[0].radTdist1": 5.01225e-7,
    "CA1_PC_Tomko[0].radTdist2": 5.01225e-7,
    "CA1_PC_Tomko[0].rad_t1": 5.01225e-7,
    "CA1_PC_Tomko[0].rad_t2": 5.01225e-7,
    "CA1_PC_Tomko[0].rad_t3": 5.01225e-7,
    "CA1_PC_Tomko[0].lm_thick2": 5.01225e-7,
    "CA1_PC_Tomko[0].lm_medium2": 5.01225e-7,
    "CA1_PC_Tomko[0].lm_thin2": 5.01225e-7,
    "CA1_PC_Tomko[0].lm_thick1": 5.01225e-7,
    "CA1_PC_Tomko[0].lm_medium1": 5.01225e-7,
    "CA1_PC_Tomko[0].lm_thin1": 25e-8,
}
gncx_spine_total = {
    "CA1_PC_Tomko[0].radTprox1": 7.16e-10,
    "CA1_PC_Tomko[0].radTprox2": 7.16e-10,
    "CA1_PC_Tomko[0].radTmed1": 7.16e-10,
    "CA1_PC_Tomko[0].radTmed2": 7.16e-10,
    "CA1_PC_Tomko[0].radTdist1": 7.16e-10,
    "CA1_PC_Tomko[0].radTdist2": 7.16e-10,
    "CA1_PC_Tomko[0].rad_t1": 7.16e-10,
    "CA1_PC_Tomko[0].rad_t2": 7.16e-10,
    "CA1_PC_Tomko[0].rad_t3": 7.16e-10,
    "CA1_PC_Tomko[0].lm_thick2": 7.16e-10,
    "CA1_PC_Tomko[0].lm_medium2": 7.16e-10,
    "CA1_PC_Tomko[0].lm_thin2": 7.16e-10,
    "CA1_PC_Tomko[0].lm_thick1": 7.16e-10,
    "CA1_PC_Tomko[0].lm_medium1": 7.16e-10,
    "CA1_PC_Tomko[0].lm_thin1": 7.16e-10
}  # 1/7 of dends  https://doi.org/10.1073/pnas.0605412104 
Ca_Ext = 2
n_seg = 1


#Ca indicators

#For Sabatinis life cycle of a Ca2+ ion Ca in the spine tuning
tot_magnesium_green_BS = 0.100
magnesium_green_bound = 0.0014972
kf_magnesium_green = 9e1#1/mM/ms from doi: 10.1016/j.ymeth.2008.09.025
kb_magnesium_green = 0.6 #1/ms", Kd in vitro 6uM", doi: 10.1016/j.ymeth.2008.09.025 and Sabatini
mggreenDiff = 15e-3   #15 um^2/s doi:10.1016/S0006-3495(96)79633-9 

diffusions = {"CaM": camDiff,
              "Calb": calbDiff,
              "Mg Green": mggreenDiff}
membrane_shell_width = 0.1
vol_c = 34.19592685263736
