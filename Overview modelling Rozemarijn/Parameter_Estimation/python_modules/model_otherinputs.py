from scipy.optimize import fsolve
from scipy.integrate import odeint

        # Functions used by all models
def run_ss_stairs(m, inits, total_protein, learned_params):
    ss = fsolve(m, inits, args=(0, [0,0,0], [0,0,0], total_protein, learned_params))
    return ss

def simulate_wt_experiment_stairs(m, timepoints, levels, inits, total_protein, learned_params, time, run_type=None):
    odes = odeint(m, inits, time, args=(timepoints, levels, total_protein, learned_params, run_type))
    return odes

# Signal functions

def signal_ramp(timepoints, levels, t):
    time1, time2= timepoints
    level1, level2 = levels    
    sig = 0
    if t >= time1:
        leveldiff = level2-level1
        timediff = time2-time1
        ramp = leveldiff/timediff
        sig = ramp*(t-time1)
    if t >= time2:
        sig = level2
    return sig


def signal_stairs(timepoints, levels, t):
    time1, time2, time3 = timepoints
    level1, level2, level3 = levels
    sig = 0
    if t >= time1:
        sig = level1
    if t >= time2:
        sig = level2
    if t >= time3:
        sig = level3
    return sig

# Specific model functions

class Model():
    def __init__(self, m, sho1=None, sln1=None):
        self.m = m
        self.sho1 = sho1
        self.sln1 = sln1

#OptimizedModel1, import & export non-MM, stairs
def OM1_impexp_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9A, k9B, k10, K10, k11, k12, k13, k14 = params #18

    sig = signal_stairs(timepoints, levels, t)   
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
        
    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - k8 * Hog1_AC + k9B * Hog1_AN
    dHog1_AN = k8 * Hog1_AC - k9B * Hog1_AN - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN + k11 * Hog1_IC
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#OptimizedModel2
def OM2_stairs(initials,t,timepoints,levels,total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15 = params #18
    
    sig = signal_stairs(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * sig- Glycerol) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * sig - Glycerol) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#OptimizedModel1, Hill function
def OM1_Hill_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, beta = params #18

    sig = signal_stairs(timepoints, levels, t)   
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


#OptimizedModel, Sho1exp
def OM_Sho1exp_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15 = params #18

    sig = signal_ramp(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0
        
    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - k8 * Hog1_AC / (K8 + Hog1_AC) + (k9B - k15 * Sho1) * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = k8 * Hog1_AC / (K8 + Hog1_AC) - (k9B - k15 * Sho1) * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#OptimizedModel, Sho1imp
def OM_Sho1imp_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15 = params #18

    sig = signal_ramp(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0
        
    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * Sho1) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * Sho1) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#OptimizedModel1
def OM1_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14 = params #18

    sig = signal_stairs(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0
        
    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - k8 * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = k8 * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol


#OptimizedModel2
def OM2_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15 = params #18

    sig = signal_stairs(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * sig- Glycerol) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * sig - Glycerol) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#OptimizedModel, osmo-export
def OM_osmoexp_stairs(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, K8, k9A, K9a, k9B, K9b, k10, K10, k11, K11, k12, k13, k14, k15 = params #18

    sig = signal_ramp(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0
        
    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (base_osmo + k2 * sig - Glycerol) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - k8 * Hog1_AC / (K8 + Hog1_AC) + (k9B - k15 * sig - Glycerol) * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = k8 * Hog1_AC / (K8 + Hog1_AC) - (k9B - k15 * sig - Glycerol) * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

#Translocation123456
def T23456(initials,t, timepoints, levels, total_protein,params, run_type=None):
    Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k10, K10, k11, K11, kA, KA, k6, K6, kB, KB, k4, K4, k2, K2, k5, K5, k12, k13, k3, K3 = params #18
    
    sig = signal_stairs(timepoints, levels, t)
    Hog1_IC = Hog1_tot - Hog1_AC - Hog1_AN - Hog1_IN
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0
        
    dSho1 = (base_osmo + k10 * sig - Glycerol) * (Sho1_inactive) / (K10 + Sho1_inactive) - k11 * Sho1 / (K11 + Sho1)
    dHog1_AC = kA * Sho1 * Hog1_IC / (KA + Hog1_IC) - k6 * Hog1_AC / (K6 + Hog1_AC) - kB * Hog1_AC / (KB + Hog1_AC)
    dHog1_AN = kB * Hog1_AC / (KB + Hog1_AC) - k2 * Hog1_AN / (K2 + Hog1_AN) + k3 * Sho1 * Hog1_IN / (K3 + Hog1_IN)
    dHog1_IN = k2 * Hog1_AN / (K2 + Hog1_AN) - k4 * Hog1_IN / (K4 + Hog1_IN) + k5 * Hog1_IC / (K5 + Hog1_IC) - k3 * Sho1 * Hog1_IN / (K3 + Hog1_IN)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol