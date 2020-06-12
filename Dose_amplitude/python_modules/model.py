from scipy.optimize import fsolve
from scipy.integrate import odeint
import numpy as np

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
        sig = 250
    if t_step >= 20:
        sig = 550
    return sig

tttt = 18


def get_ramp_signal(t_step, steps):
    dt = 0.1
    s = 601
    time = np.linspace(0,dt*s,s)
    if t_step == 0:
        sig = 0
        return sig
    # print(steps)
    num_steps = len(steps)
    interval = len(time)/num_steps
    i = int(np.floor(t_step/dt/interval))
    if i >= num_steps:
        sig = sum(steps)
    else:
        sig = sum(steps[0:i+1])
    # print(sig)
    if t_step % 10 == 0:
        print(sig)
    return sig

def get_manual_signal(t_step):
    sig = 0
    if t_step >= .00001:
        sig = 200
    if t_step >= 10:
        sig = 300
    if t_step >= 12:
        sig = 375
    if t_step >= 25:
        sig = 500
    if t_step >= 45:
        sig = 550
    return sig



class Model():
    def __init__(self, m, nopos=None, inhib=None):
        self.m = m
        # self.t100a = t100a
        # self.nopos = nopos
        # self.inhib = inhib
# Specific model functions
# M1
def M_ff_cyto(initials,t,total_protein,sig,params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            sig = signal_ramp_special(t)
        elif run_type[0] == 'rand':
            sig = get_ramp_signal(t, run_type[1])
        elif run_type[0] == 'man':
            sig = get_manual_signal(t)
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
    dHog1_AC = (k5 * Sln1 + k6 * Sho1 + k14 * Hog1_AC) * Hog1_IC / (K56 + Hog1_IC) - k7 * Hog1_AC / (K7 + Hog1_AC) - (k8 + k15 * sig - Glycerol) * Hog1_AC / (K8 + Hog1_AC) + k9B * Hog1_AN / (K9b + Hog1_AN)
    dHog1_AN = (k8 + k15 * sig - Glycerol) * Hog1_AC / (K8 + Hog1_AC) - k9B * Hog1_AN / (K9b + Hog1_AN) - k10 * Hog1_AN / (K10 + Hog1_AN)
    dHog1_IN = k10 * Hog1_AN / (K10 + Hog1_AN) - k9A * Hog1_IN / (K9a + Hog1_IN) + k11 * Hog1_IC / (K11 + Hog1_IC)
    dGlycerol = k12 * Hog1_AC - k13 * Glycerol

    return dSln1, dSho1, dHog1_AC, dHog1_AN, dHog1_IN, dGlycerol

# M2
def M_ff_nuc(initials,t,total_protein,sig,params, run_type=None):
    if run_type:
        if run_type[0] == 'ramp':
            sig = signal_ramp_special(t)
        elif run_type[0] == 'rand':
            sig = get_ramp_signal(t, run_type[1])
        elif run_type[0] == 'man':
            sig = get_manual_signal(t)

    Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol = initials
    Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k10, K10, k11, K11, kA, KA, k6, K6, kB, KB, k4, K4, k2, K2, k5, K5, k12, k13, k3, K3 = params #18

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
