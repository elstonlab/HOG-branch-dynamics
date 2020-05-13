from scipy.optimize import fsolve
from scipy.integrate import odeint

class Model():
    def __init__(self, m):
        self.m = m

def run_wt_ss(model, inits, total_protein, learned_params):
    ss = fsolve(model, inits, args=(0,total_protein, 0, learned_params))
    return ss

def run_sln1_ss(model, inits, total_protein, learned_params):
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, Sho1_on = total_protein
    total_protein = Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, 0
    ss = fsolve(model, inits, args=(0,total_protein, 0, learned_params))
    return ss

def run_sho1_ss(model, inits, total_protein, learned_params):
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, Sho1_on = total_protein
    total_protein = Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, 0, Sho1_on
    ss = fsolve(model, inits, args=(0,total_protein, 0, learned_params))
    return ss

def simulate_wt_experiment(model, inits, total_protein, sig, learned_params, time):
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def simulate_sln1_experiment(model, inits, total_protein, sig, learned_params, time):
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, Sho1_on = total_protein
    total_protein = Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, 0
    #solve odes:
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def simulate_sho1_experiment(model, inits, total_protein, sig, learned_params, time):
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, Sho1_on = total_protein
    total_protein = Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, 0, Sho1_on
    #solve odes:
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def M8(initials,t,total_protein,sig,params):
    Sln1, Sho1, Pbs2A, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, Sho1_on = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, K5, k6, K6, k7, K7, k8, k10, k11, K11, k12, K12 = params #18

    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Pbs2I = Pbs2_tot - Pbs2A
    Hog1I = Hog1_tot - Hog1A

    # checks for negative glycerol
    if Glycerol < 0:
        Glycerol = 0

    dSln1 = (k1*(base_osmo + sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = (k2*(base_osmo + sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dPbs2 =  Sln1_on*(k5 * Sln1 * Pbs2I / (K5 + Pbs2I)) + Sho1_on*(k6 * Sho1 * Pbs2I / (K6 + Pbs2I)) - k7*Pbs2A / (K7 + Pbs2A)
    dHog1A = (k11 * Pbs2A * Hog1I / (K11 + Hog1I)) - k12*Hog1A / (K12 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol
    # dGlycerol = k8 * Hog1A  - k10 * Glycerol

    return dSln1, dSho1, dPbs2, dHog1A, dGlycerol

def M10(initials,t,total_protein,sig,params):
    Sln1, Sho1, Pbs2A, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Pbs2_tot, Hog1_tot, _, Sln1_on, Sho1_on = total_protein
    kb1, kb2, b1, b2, k1, K1, k2, K2, k3, K3, k4, K4, k5, K5, k6, K6, k7, K7, k8, k10, k11, K11, k12, K12 = params #18

    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Pbs2I = Pbs2_tot - Pbs2A
    Hog1I = Hog1_tot - Hog1A


    if Glycerol < 0:
        Glycerol = 0
    # checks for negative glycerol
    # if Glycerol < 0:
    #     Glycerol = 0

    dSln1 = ((kb1 + k1 *sig)/(1+Glycerol/b1))* Sln1_inactive / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((kb2 + k2 *sig)/(1+Glycerol/b2))* Sho1_inactive / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dPbs2 =  Sln1_on*(k5 * Sln1 * Pbs2I / (K5 + Pbs2I)) + Sho1_on*(k6 * Sho1 * Pbs2I / (K6 + Pbs2I)) - k7*Pbs2A / (K7 + Pbs2A)
    dHog1A = (k11 * Pbs2A * Hog1I / (K11 + Hog1I)) - k12*Hog1A / (K12 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol
    # dGlycerol = k8 * Hog1A  - k10 * Glycerol

    return dSln1, dSho1, dPbs2, dHog1A, dGlycerol
