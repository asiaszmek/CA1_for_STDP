head_factor = 1
ca_factor = 1
c_unit = 6.0221409e5
ip3rtau = 2000 #wagner
ip3degTau = 1000
gip3r = 120400

gleak = 15.905e-3
ip3_init = 0 # 50 nM? signalling pathways model

calr_tot = 86 #  calreticulin
calr_bound = 7.11945
kf_calr = 0.1 # 1/mM/ms
kb_calr = 0.2 #  1/ms

Ca_per = .1 # 10% at 2mM, 5% at 1mM, 20% at 4mM extracellular Ca  doi: 10.1073/pnas.90.24.11573
ca_i_f = 1  # increase VGCC conductance by a factor for sections with rxd Calcium
ca_init = 50e-6
ca_init_ER = .180512
ca_ECS = 2.

caDiff = 0.2
ip3Diff = 0.1
camDiff = 0.066
calbDiff = 0.066

AtoN_ratio = 1
gAMPA = 25e-3
gNMDA = gAMPA/AtoN_ratio


Kip3 = 0.15e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
Kact = 0.8e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
k_inh = 2e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
Km_Serca = 0.13e-3 # Michaelis constant for SERCA pump
kcat_Serca = 7.5e-3 #ms
gSerca = 1#signaling pathways model
ip3r_gate_state = 0.8
calbindin_tot = 0.150 #(150 uM)
calmodulin_tot = 0.03 #(30 uM)
fixed_buffer_tot = 2 #(2 mM)


kf_calbindin = 2.8e-2
kb_calbindin = 0.0196

kf_camn = 2*7.7e2
kb_camn = 1.6e2
kf_camc = 2*8.4e1
kb_camc = 2.6
kf_fixed_b = 400
kb_fixed_b = 40
fixed_buffer_ca = kf_fixed_b*ca_init*fixed_buffer_tot/kb_fixed_b

camn = kf_camn*ca_init*calmodulin_tot/kb_camn
camc = kf_camc*ca_init*calmodulin_tot/kb_camc
calbca = kf_calbindin*ca_init*calbindin_tot/kb_calbindin


kf_pmca = 50
kb_pmca = 0.007
kcat_pmca = 0.0035
Km_pmca = (kb_pmca+kcat_pmca)/kf_pmca
gpmca = 0.1
gpmca_spine = 50e-3*ca_factor # {"apical_dendrite[10]": 0.1e-5*ca_factor}

ncx_pow = 1
kf_ncx = 16.8
kb_ncx = 0.0112
kcat_ncx = 0.0056
Km_ncx = (kb_ncx+kcat_ncx)/kf_ncx
#  this dynamics is more similar to quasi-steady state approx
gncx = 3
gncx_spine = gncx/7 #{}
# for key in gncx:
#     gncx_spine[key] = gncx[key]/7
# 1/7 of dends  https://doi.org/10.1073/pnas.0605412104 
Ca_Ext = 2
n_seg = 1


#Ca indicators

#For Sabatinis life cycle of a Ca2+ ion Ca in the spine tuning
tot_magnesium_green_BS = 0.100
magnesium_green_bound = 0.0014972
kf_magnesium_green = 9e1#1/mM/ms from doi: 10.1016/j.ymeth.2008.09.025
kb_magnesium_green = 0.6 #1/ms, Kd in vitro 6uM, doi: 10.1016/j.ymeth.2008.09.025 and Sabatini
mggreenDiff = 15e-3   #15 um^2/s doi:10.1016/S0006-3495(96)79633-9 

#For Marsden et al. 10.1073/pnas.1010346107
#Fluo-3, constants Pflügers Arch – Eur J Physiol (1997) 434:615–631
tot_fluo3 = 1e-3 #  1 uM
kf_fluo3 = 6e2 # In agreement with Roger Tsien
kb_fluo3 = 0.6   # 1/ms
fluo3Diff = 15e-3 # um^2/ms take Fluo5's

tot_BF2 = 0.1 # bis-Fura-2, 100 uM, Frick, Migliore, Johnston
kf_BF2 = 6e2 # 0.53uM DOI: 10.1117/1.NPh.2.2.021010 
kb_BF2 = 0.318
BF2Diff = .66

tot_OGB1 = 0.2 # Oregon Green 488 BAPTA-1
kf_OGB1 = 6e2
kb_OGB1 = 0.126
OGB1Diff = .66

diffusions = {"Calmodulin": camDiff,
              "Calbindin": calbDiff,
              "Mg Green": mggreenDiff,
              "Fluo3": fluo3Diff,
              "BF2": BF2Diff,
              "OGB1": OGB1Diff,
              "Fixed": 0,
              }
              
membrane_shell_width = .1


gbar_kca = {
    "soma": 0.0015,
    "radTprox1": 9.03e-05,
    "radTprox2": 9.03e-05,
    "radTmed1": 9.03e-05,
    "radTmed2": 9.03e-05,
    "radTdist1": 9.03e-05,
    "radTdist2": 9.03e-05,
    "lm_thick1": 9.03e-05,
    "lm_thick2": 9.03e-05,
    "lm_medium1": 9.03e-05,
    "lm_medium2": 9.03e-05,
    "lm_thin1": 9.03e-05,
    "lm_thin2": 9.03e-05,
    "rad_t1":  9.03e-05,
    "rad_t2":  9.03e-05,
    "rad_t3":  9.03e-05,
    "oriprox1": 9.03e-05,
    "oridist1_1": 9.03e-05,
    "oridist1_2": 9.03e-05,
    "oriprox2": 9.03e-05,
    "oridist2_1": 9.03e-05,
    "oridist2_2": 9.03e-05,
}

gbar_cagk = {
    "soma": 4.48e-05,
    "radTprox1": 4.48e-05,
    "radTprox2": 4.48e-05,
    "radTmed1": 4.48e-05,
    "radTmed2": 4.48e-05,
    "radTdist1": 4.48e-05,
    "radTdist2": 4.48e-05,
    "lm_thick1": 4.48e-05,
    "lm_thick2": 4.48e-05,
    "lm_medium1": 4.48e-05,
    "lm_medium2": 4.48e-05,
    "lm_thin1": 4.48e-05,
    "lm_thin2": 4.48e-05,
    "rad_t1":  4.48e-05,
    "rad_t2":  4.48e-05,
    "rad_t3":  4.48e-05,
    "oriprox1": 4.48e-05,
    "oridist1_1": 4.48e-05,
    "oridist1_2": 4.48e-05,
    "oriprox2": 4.48e-05,
    "oridist2_1": 4.48e-05,
    "oridist2_2": 4.48e-05,
}
