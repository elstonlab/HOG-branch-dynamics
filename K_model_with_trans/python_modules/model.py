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
    def __init__(self, m, t100a=None,sho1=None, sln1=None):
        self.m = m
        self.t100a = t100a
        self.sho1 = sho1
        self.sln1 = sln1

def simulate_t100a_experiment_M2c_kb(m, inits, total_protein, sig, learned_params, time,  run_type=None):
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6 = learned_params #17
    learned_params = beta_3, 0, kb, k1, k3, k5, 0, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6
    #solve odes:
    odes = odeint(m, inits, time, args=(total_protein, sig, learned_params, run_type))
    return odes

def M2c_kb(initials,t,total_protein,sig,params,run_type=None):
    # print(initials)
    if run_type:
        if run_type[0] == 'ramp':
            sig = signal_ramp_special(t)
        elif run_type[0] == 'rand':
            sig = get_ramp_signal(t, run_type[1])
        elif run_type[0] == 'man':
            sig = get_manual_signal(t)

    MAP3K, MAP2K, MAPK, gly = initials
    MAP3K_t, MAP2K_t, MAPK_t, _ = total_protein
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6 = params #17

    MAP3K_I = MAP3K_t-MAP3K
    MAP2K_I = MAP2K_t-MAP2K
    MAPK_I = MAPK_t-MAPK

    dMAP3K = (((sig*k1 + kb)/(1+gly/beta_3))*MAP3K_I)/(K_1+MAP3K_I) - (k2*MAP3K/(K_2+MAP3K))
    dMAP2K = ((k3*MAP3K*MAP2K_I)/(K_3+MAP2K_I)) - (k4*MAP2K/(K_4+MAP2K))
    dMAPK = (((k5*MAP2K + MAPK*alpha))*MAPK_I)/(K_5+MAPK_I) - (k6*MAPK)/(K_6+MAPK)  #bug
    dgly = s7*MAPK - d8*gly
    return dMAP3K, dMAP2K, dMAPK, dgly

def simulate_t100a_experiment_M2c_kb_nuc(m, inits, total_protein, sig, learned_params, time,  run_type=None):
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6, n1, n2  = learned_params #17
    learned_params = beta_3, 0, kb, k1, k3, k5, 0, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6, 0, n2
    #solve odes:
    odes = odeint(m, inits, time, args=(total_protein, sig, learned_params, run_type))
    return odes


def M2c_kb_nuc(initials,t,total_protein,sig,params,run_type=None):
    # print(initials)
    if run_type:
        if run_type[0] == 'ramp':
            sig = signal_ramp_special(t)
        elif run_type[0] == 'rand':
            sig = get_ramp_signal(t, run_type[1])
        elif run_type[0] == 'man':
            sig = get_manual_signal(t)

    MAP3K, MAP2K, MAPK, gly, MAPK_n = initials
    MAP3K_t, MAP2K_t, MAPK_t, _ = total_protein
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6, n1, n2 = params #17

    MAP3K_I = MAP3K_t-MAP3K
    MAP2K_I = MAP2K_t-MAP2K
    MAPK_I = MAPK_t-MAPK-MAPK_n

    dMAP3K = (((sig*k1 + kb)/(1+gly/beta_3))*MAP3K_I)/(K_1+MAP3K_I) - (k2*MAP3K/(K_2+MAP3K))
    dMAP2K = ((k3*MAP3K*MAP2K_I)/(K_3+MAP2K_I)) - (k4*MAP2K/(K_4+MAP2K))
    dMAPK = (((k5*MAP2K + MAPK*alpha))*MAPK_I)/(K_5+MAPK_I) - (k6*MAPK)/(K_6+MAPK)  + n2*MAPK_n - n1*MAPK
    dMAPK_n = n1*MAPK - n2*MAPK_n
    dgly = s7*MAPK - d8*gly
    return dMAP3K, dMAP2K, dMAPK, dMAPK_n, dgly

def simulate_t100a_experiment_M2c_kb_nuc_ptrans(m, inits, total_protein, sig, learned_params, time,  run_type=None):
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6, n1, n2  = learned_params #17
    learned_params = beta_3, 0, kb, k1, k3, k5, 0, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6, n1, n2
    #solve odes:
    odes = odeint(m, inits, time, args=(total_protein, sig, learned_params, run_type))
    return odes


def M2c_kb_nuc_ptrans(initials,t,total_protein,sig,params,run_type=None):
    # print(initials)
    if run_type:
        if run_type[0] == 'ramp':
            sig = signal_ramp_special(t)
        elif run_type[0] == 'rand':
            sig = get_ramp_signal(t, run_type[1])
        elif run_type[0] == 'man':
            sig = get_manual_signal(t)

    MAP3K, MAP2K, MAPK, gly, MAPK_n = initials
    MAP3K_t, MAP2K_t, MAPK_t, _ = total_protein
    beta_3, alpha, kb, k1, k3, k5, s7, k2, k4, k6, d8, K_1, K_3, K_5, K_2, K_4, K_6, n1, n2 = params #17

    MAP3K_I = MAP3K_t-MAP3K
    MAP2K_I = MAP2K_t-MAP2K
    MAPK_I = MAPK_t-MAPK-MAPK_n

    dMAP3K = (((sig*k1 + kb)/(1+gly/beta_3))*MAP3K_I)/(K_1+MAP3K_I) - (k2*MAP3K/(K_2+MAP3K))
    dMAP2K = ((k3*MAP3K*MAP2K_I)/(K_3+MAP2K_I)) - (k4*MAP2K/(K_4+MAP2K))
    dMAPK = (((k5*MAP2K + MAPK*alpha))*MAPK_I)/(K_5+MAPK_I) - (k6*MAPK)/(K_6+MAPK)  + n2*MAPK_n - n1*MAPK
    dMAPK_n = n1*MAPK - n2*MAPK_n
    dgly = s7*MAPK - d8*gly
    return dMAP3K, dMAP2K, dMAPK, dMAPK_n, dgly
