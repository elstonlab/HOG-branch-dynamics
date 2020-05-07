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