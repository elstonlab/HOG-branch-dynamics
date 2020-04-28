from scipy.optimize import fsolve
from scipy.integrate import odeint

        # Functions used by all models
def run_ss(m, inits, total_protein, learned_params):
    ss = fsolve(m, inits, args=(0,total_protein, 0, learned_params))
    return ss

def simulate_wt_experiment(m, inits, total_protein, sig, learned_params, time, run_type=None):
    odes = odeint(m, inits, time, args=(total_protein, sig, learned_params, run_type))
    return odes

# Signal function
def signal_ramp_special(t_step):
    sig = 0
    if t_step >= .001:
        sig = 300000
    if t_step >= 40:
        sig = 450000
    if t_step >= 45:
        sig = 550000
    return sig

# Specific model functions

class Model():
    def __init__(self, m, sho1=None, sln1=None):
        self.m = m
        self.sho1 = sho1
        self.sln1 = sln1


# M16
def M16(initials, t, total_protein, salt, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12 = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = k1 * (base_osmo + salt - Glycerol) * (Sln1_inactive) - k3 * Sln1
    dSho1 = k2 * (base_osmo + salt - Glycerol) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

def simulate_sho1D_experiment_M16(m, inits, total_protein, sig, learned_params, time,  run_type=None):
# ADJUST
    # MAP3K_t, MAP2K_t, MAPK_t, fb = total_protein
    # total_protein = MAP3K_t, MAP2K_t, MAPK_t, 0
    # odes = odeint(m, inits, time, args=(total_protein, sig, learned_params, run_type))
    return odes

# MK
def MK(initials,t,total_protein,sig,params, run_type=None):
    MAP3K, MAP2K, MAPK, gly = initials
    MAP3K_t, MAP2K_t, MAPK_t, _ = total_protein
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6 = params #16

    MAP3K_I = MAP3K_t-MAP3K
    MAP2K_I = MAP2K_t-MAP2K
    MAPK_I = MAPK_t-MAPK
    # PTP_I = PTP_t-PTP

    dMAP3K = (((sig*k1 + kb)/(1+gly/beta_3))*MAP3K_I)/(K_1+MAP3K_I) - (k2*MAP3K/(K_2+MAP3K))
    dMAP2K = (((k3*MAP3K)*MAP2K_I)/(K_3+MAP2K_I)) - (k4*MAP2K/(K_4+MAP2K))
    dMAPK = ((k5*MAP2K + MAPK*alpha)*MAPK_I)/(K_5+MAPK_I) - (k6*MAPK)/(K_6+MAPK)
    dgly = s7*MAPK - d8*gly

    return dMAP3K, dMAP2K, dMAPK, dgly

# M15
def M15(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 = params #15

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = k1 * (base_osmo + sig - Glycerol) * (Sln1_inactive) - k3 * Sln1
    dSho1 = k2 * (base_osmo + sig - Glycerol) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M15 Hill
def M15H(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, beta = params #16

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = (base_osmo + k1 * sig/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = (base_osmo + k2 * sig/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M15 Hill2
def M15H2(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, beta = params #16

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M16 Hill
def M16H(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = (base_osmo + k1 * sig/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = (base_osmo + k2 * sig/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M16 Hill2
def M16H2(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M17
def M17(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k14, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M18
def M18(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k15, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + k15 * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - k15 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M19
def M19(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k14, k15, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + k15 * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - k15 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M20
def M20(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M21
def M21(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k15, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + k15 * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - k15 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M22
def M22(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + k15 * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - k15 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M22_lin
def M22_lin(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15 = params #13

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = k1 * (base_osmo + sig - Glycerol) * (Sln1_inactive) - k3 * Sln1
    dSho1 = k2 * (base_osmo + sig - Glycerol) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + k15 * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - k15 * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M23
def M23(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta, alpha = params #15

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + (k15 + alpha * Hog1_AN) * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - (k15 + alpha * Hog1_AN) * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M24
def M24(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta, alpha = params #15

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + (k15 + alpha * Hog1_AC) * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - (k15 + alpha * Hog1_AC) * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M25
def M25(initials, t, total_protein, sig, params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            salt = signal_ramp_special(t)

    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta, alpha = params #15

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k13 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + (k15 + alpha * (Hog1_AC + Hog1_AN)) * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - (k15 + alpha * (Hog1_AC + Hog1_AN)) * Hog1_AN
    dHog1_IN = k9 * Hog1_AN - k10 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k11 * Hog1_AC - k12 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model 26
def M26(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, beta = params #16

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - k8 * Hog1_AC + k9 * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9 * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9 * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model 27
def M27(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9A, k9B, k10, k11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * Hog1_AC) * Hog1_AC + k9B * Hog1_AN
    dHog1_AN = (k8 + k15 * Hog1_AC) * Hog1_AC - k9B * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9A * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model 28
def M28(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9A, k9B, k10, k11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k9B * Hog1_AN
    dHog1_AN = (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k9B * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9A * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model 29
def M29(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9A, k9B, k10, k11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k9B * Hog1_AN
    dHog1_AN = (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k9B * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9A * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model 30
def M30(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * Hog1_AC) * Hog1_AC + k9 * Hog1_AN
    dHog1_AN = (k8 + k15 * Hog1_AC) * Hog1_AC - k9 * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9 * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model31
def M31(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k9 * Hog1_AN
    dHog1_AN = (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k9 * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9 * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model32
def M32(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k9 * Hog1_AN
    dHog1_AN = (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k9 * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9 * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model33
def M33(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9A, k9B, k10, k11, k12, k13, k14, k15, k16, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC + (k9B + k16 * Hog1_AN) * Hog1_AN
    dHog1_AN = (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC - (k9B + k16 * Hog1_AN) * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - k9A * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model34
def M34(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC - k7 * Hog1_AC - (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC + (k9 + k16 * Hog1_AN) * Hog1_AN
    dHog1_AN = (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC - (k9 + k16 * Hog1_AN) * Hog1_AN - k10 * Hog1_AN
    dHog1_IN = k10 * Hog1_AN - (k9 + k16 * Hog1_AN) * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model29_MM
def M29MM(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model29_Pbs2
def M29Pbs2(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Pbs2, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Pbs2_inactive = Pbs2_tot - Pbs2
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dPbs2 = (k5 * Sln1 + k6 * Sho1) * Pbs2_inactive - k7 * Pbs2
    dHog1_AC = (k8 * Pbs2 + k17 * Hog1_AC) * Hog1_IC - k9 * Hog1_AC - (k10 + k18 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k11 * Hog1_AN
    dHog1_AN = (k10 + k18 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k11 * Hog1_AN - k12 * Hog1_AN
    dHog1_IN = k12 * Hog1_AN - k13 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k15 * Hog1_AC - k16 * Glycerol

    return dSln1, dSho1, dPbs2, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model29, Pbs2, MM
def M29Pbs2MM(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Pbs2, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9, K9, k10, K10, k11, K11, k12, K12, k13, K13, k14, K14, k15, k16, k17, k18, beta = params #33

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Pbs2_inactive = Pbs2_tot - Pbs2
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dPbs2 = (k5 * Sln1 + k6 * Sho1) * Pbs2_inactive / (K56 + Pbs2_inactive) - k7 * Pbs2 / (K7 + Pbs2)
    dHog1_AC = (k8 * Pbs2 + k17 * Hog1_AC) * Hog1_IC / (K8 + Hog1_IC) - k9 * Hog1_AC / (K9 + Hog1_AC) - (k10 + k18 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K10 + Hog1_AC) + k11 * Hog1_AN / (K11 + Hog1_AN)
    dHog1_AN = (k10 + k18 * Hog1_AC * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K10 + Hog1_AC) - k11 * Hog1_AN / (K11 + Hog1_AN) - k12 * Hog1_AN / (K12 + Hog1_AN)
    dHog1_IN = k12 * Hog1_AN / (K12 + Hog1_AN) - k13 * Hog1_IN / (K13 + Hog1_IN) + k14 * Hog1_IC / (K14 + Hog1_IC)
    dGlycerol = k15 * Hog1_AC - k16 * Glycerol

    return dSln1, dSho1, dPbs2, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model28, Pbs2
def M28Pbs2(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Pbs2, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Pbs2_inactive = Pbs2_tot - Pbs2
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) - k3 * Sln1
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) - k4 * Sho1
    dPbs2 = (k5 * Sln1 + k6 * Sho1) * Pbs2_inactive - k7 * Pbs2
    dHog1_AC = (k8 * Pbs2 + k17 * Hog1_AC) * Hog1_IC - k9 * Hog1_AC - (k10 + k18 * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k11 * Hog1_AN
    dHog1_AN = (k10 + k18 * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k11 * Hog1_AN - k12 * Hog1_AN
    dHog1_IN = k12 * Hog1_AN - k13 * Hog1_IN + k14 * Hog1_IC
    dGlycerol = k15 * Hog1_AC - k16 * Glycerol

    return dSln1, dSho1, dPbs2, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model28, Pbs2, MM
def M28Pbs2MM(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Pbs2, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9, K9, k10, K10, k11, K11, k12, K12, k13, K13, k14, K14, k15, k16, k17, k18, beta = params #33

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Pbs2_inactive = Pbs2_tot - Pbs2
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dPbs2 = (k5 * Sln1 + k6 * Sho1) * Pbs2_inactive / (K56 + Pbs2_inactive) - k7 * Pbs2 / (K7 + Pbs2)
    dHog1_AC = (k8 * Pbs2 + k17 * Hog1_AC) * Hog1_IC / (K8 + Hog1_IC) - k9 * Hog1_AC / (K9 + Hog1_AC) - (k10 + k18 * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K10 + Hog1_AC) + k11 * Hog1_AN / (K11 + Hog1_AN)
    dHog1_AN = (k10 + k18 * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K10 + Hog1_AC) - k11 * Hog1_AN / (K11 + Hog1_AN) - k12 * Hog1_AN / (K12 + Hog1_AN)
    dHog1_IN = k12 * Hog1_AN / (K12 + Hog1_AN) - k13 * Hog1_IN / (K13 + Hog1_IN) + k14 * Hog1_IC / (K14 + Hog1_IC)
    dGlycerol = k15 * Hog1_AC - k16 * Glycerol

    return dSln1, dSho1, dPbs2, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model28, MM
def M28MM(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model28.2
def M28_2(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K5, K6, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, K14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = k5 * Sln1 * Hog1_IC / (K5 + Hog1_IC) + k6 * Sho1 * Hog1_IC / (K6 + Hog1_IC) + k14 * Hog1_AC * Hog1_IC / (K14 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model28.3
def M28_3(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K5, K6, k7, K7, k8, k9A, k9B, k10, K10, k11, k12, k13, k14, K14, k15, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = k5 * Sln1 * Hog1_IC / (K5 + Hog1_IC) + k6 * Sho1 * Hog1_IC / (K6 + Hog1_IC) + k14 * Hog1_AC * Hog1_IC / (K14 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC + k9B * Hog1_AN
    dHog1_AN = (k8 + k15 * ((sig)/(1+Glycerol/beta))) * Hog1_AC - k9B * Hog1_AN - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Model28, MM minus osmostress dependency
def M28MMmin(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, beta = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig)/(1+Glycerol/beta)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig)/(1+Glycerol/beta)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - k8 * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = k8 * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#OptimizedModel2
def OM2(initials,t,total_protein,sig,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15 = params #18

    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0

    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * sig- Glycerol) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * sig - Glycerol) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol