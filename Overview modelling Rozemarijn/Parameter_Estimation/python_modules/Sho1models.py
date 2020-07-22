from scipy.optimize import fsolve
from scipy.integrate import odeint

        # Functions used by all models
def run_ss(inits, total_protein, learned_params):
    ss = fsolve(model, inits, args=(0,total_protein, 0, learned_params))
    return ss

def simulate_wt_experiment(model, inits, total_protein, sig, learned_params, time):
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def simulate_sln1_experiment(model, inits, total_protein, sig, learned_params, time):
    #base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = learned_params #16
    #learned_params = base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, 0, K56, k7, K7, k8, k9, k10
    base_osmo, alpha, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k10 = learned_params #16
    learned_params = base_osmo, alpha, k1, K1, 0, 0, k3, K3, 0, 0, k5, 0, K56, k7, K7, k8, k10
    #solve odes:
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def simulate_sho1_experiment(model, inits, total_protein, sig, learned_params, time):
    #base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = learned_params #16
    #learned_params = base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, 0, k6, K56, k7, K7, k8, k9, k10
    base_osmo, alpha, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k10 = learned_params #16
    learned_params = base_osmo, alpha, 0, 0, k2, K2, 0, 0, k4, K4, 0, k6, K56, k7, K7, k8, k10
    #solve odes:
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

# Specific model functions

class Model():
    def __init__(self, m):
        self.m = m

# Model1
def M1(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    dSln1 = ((base_osmo + sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I)
    dGlycerol = k8 * Hog1A + k9 * (base_osmo + sig - Glycerol) - k10 * Glycerol
    # dGlycerol = k8 * Hog1A  - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model1_corr
def M1_corr(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + k9 * (base_osmo + sig - Glycerol) - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model1_corr2
def M1_corr2(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K5, K6, k7, K7, k8, k9, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    dSln1 = ((base_osmo + k1 * sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = k5 * Sln1 * Hog1I / (K5 + Hog1I) + k6 * Sho1 * Hog1I / (K6 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + (base_osmo + k9 * sig - Glycerol) - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Branchesmodel2 (complex formation for the two branches)
def BM2(initials, t, total_protein, sig, params, run_type=None):
    Sln1Pbs2A, Sho1Pbs2A, Hog1IPbs2A, Hog1A, Glycerol = initials
    Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5f, k5r, k6f, k6r, k7, K7, k8, k9, k10, kcat = params #18

    Hog1I = Hog1_tot - Hog1A - Hog1IPbs2A
    Pbs2I = Pbs2_tot - Sln1Pbs2A - Sho1Pbs2A - Hog1IPbs2A
    
    dSln1Pbs2A = ((base_osmo + k1 * sig - Glycerol)) * (Pbs2I) / (K1 + Pbs2I) - k3 * Sln1Pbs2A / (K3 + Sln1Pbs2A)
    dSho1Pbs2A = ((base_osmo + k2 * sig - Glycerol)) * (Pbs2I) / (K2 + Pbs2I) - k4 * Sho1Pbs2A / (K4 + Sho1Pbs2A)
    dHog1IPbs2A = k5f * Hog1I * Sln1Pbs2A + k6f * Hog1I * Sho1Pbs2A - k5r * Hog1IPbs2A - k6r * Hog1IPbs2A - kcat * Hog1IPbs2A
    dHog1A = kcat * Hog1IPbs2A - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + (base_osmo + k9 * sig - Glycerol) - k10 * Glycerol

    return dSln1Pbs2A, dSho1Pbs2A, dHog1IPbs2A, dHog1A, dGlycerol

# Model1_corr, no basal signaling Sho1
def M1_corr_kb(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + k9 * (base_osmo + sig - Glycerol) - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model1_corr, positive feedback
def M1_corr_pfb(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, alpha, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    dSln1 = ((base_osmo + k1 * sig - Glycerol) + alpha*Hog1A) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + k9 * (base_osmo + sig - Glycerol) - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model1_corr, negative feedback Sho1, no basal signaling Sho1
def M1_corr_kb_nfb(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, beta, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0

    dSln1 = ((base_osmo + k1 * sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((k2 * sig - Glycerol) - beta*Hog1A) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model1_corr, negative feedback Sho1
def M1_corr_nfb(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, beta, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1

    if Glycerol < 0:
        Glycerol = 0

    dSln1 = ((base_osmo + k1 * sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol) - beta*Hog1A) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model2_corr2
def M2_corr2(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Fsp1, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, Fsp1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K5, K6, k7, K7, k8, k9, K9, k10, K10, k11, k12 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Fsp1_inactive = Fsp1_tot- Fsp1
    
    if Glycerol < 0:
        Glycerol = 0
    
    dSln1 = ((base_osmo + k1 * sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = k5 * Sln1 * Hog1I / (K5 + Hog1I) + k6 * Sho1 * Hog1I / (K6 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dFsp1 = (base_osmo + k9 * sig - Glycerol) * Fsp1_inactive / (K9 + Fsp1_inactive) - k10 * Fsp1 / (K10 + Fsp1)
    dGlycerol = k8 * Hog1A + k11 * Fsp1 - k12 * Glycerol

    return dSln1, dSho1, dHog1A, dFsp1, dGlycerol\

# Model2_corr
def M2_corr(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Fps1, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, Fps1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, K9, k10, K10, k11, k12 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    Fps1_inactive = Fps1_tot - Fps1

    if Glycerol < 0:
        Glycerol = 0
    
    dSln1 = ((base_osmo + k1 * sig - Glycerol)) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dFps1 = k9 * (base_osmo + sig - Glycerol) * Fps1_inactive / (K9 + Fps1_inactive) - k10 * Fps1 / (K10 + Fps1) 
    dGlycerol = k8 * Hog1A + k11 * Fps1 - k12 * Glycerol

    return dSln1, dSho1, dHog1A, dFps1, dGlycerol

# Model1_corr, no basal signaling Sho1, no HOG-independent glycerol feedback
def M1_corr_kb_v2(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    if Glycerol < 0:
        Glycerol = 0

    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Model1_corr, positive feedback, no HOG-independent glycerol feedback
def M1_corr_pfb_v2(initials, t, total_protein, sig, params, run_type=None):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, alpha, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    if Glycerol < 0:
        Glycerol = 0

    dSln1 = ((base_osmo + k1 * sig - Glycerol) + alpha*Hog1A) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((base_osmo + k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

# Branchesmodel2 (complex formation for the two branches), normal glycerol-independent feedback
def BM2_v2(initials, t, total_protein, sig, params, run_type=None):
    Sln1Pbs2A, Sho1Pbs2A, Hog1IPbs2A, Hog1A, Glycerol = initials
    Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5f, k5r, k6f, k6r, k7, K7, k8, k9, k10, kcat = params #18

    Hog1I = Hog1_tot - Hog1A - Hog1IPbs2A
    Pbs2I = Pbs2_tot - Sln1Pbs2A - Sho1Pbs2A - Hog1IPbs2A
    
    if Glycerol < 0:
        Glycerol = 0
    
    dSln1Pbs2A = ((base_osmo + k1 * sig - Glycerol)) * (Pbs2I) / (K1 + Pbs2I) - k3 * Sln1Pbs2A / (K3 + Sln1Pbs2A)
    dSho1Pbs2A = ((base_osmo + k2 * sig - Glycerol)) * (Pbs2I) / (K2 + Pbs2I) - k4 * Sho1Pbs2A / (K4 + Sho1Pbs2A)
    dHog1IPbs2A = k5f * Hog1I * Sln1Pbs2A + k6f * Hog1I * Sho1Pbs2A - k5r * Hog1IPbs2A - k6r * Hog1IPbs2A - kcat * Hog1IPbs2A
    dHog1A = kcat * Hog1IPbs2A - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + k9 * (base_osmo + sig - Glycerol) - k10 * Glycerol

    return dSln1Pbs2A, dSho1Pbs2A, dHog1IPbs2A, dHog1A, dGlycerol

# Branchesmodel2 (complex formation for the two branches), no glycerol-independent feedback
def BM2_v3(initials, t, total_protein, sig, params, run_type=None):
    Sln1Pbs2A, Sho1Pbs2A, Hog1IPbs2A, Hog1A, Glycerol = initials
    Pbs2_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5f, k5r, k6f, k6r, k7, K7, k8, k10, kcat = params #18

    Hog1I = Hog1_tot - Hog1A - Hog1IPbs2A
    Pbs2I = Pbs2_tot - Sln1Pbs2A - Sho1Pbs2A - Hog1IPbs2A
    
    if Glycerol < 0:
        Glycerol = 0
    
    dSln1Pbs2A = ((base_osmo + k1 * sig - Glycerol)) * (Pbs2I) / (K1 + Pbs2I) - k3 * Sln1Pbs2A / (K3 + Sln1Pbs2A)
    dSho1Pbs2A = ((base_osmo + k2 * sig - Glycerol)) * (Pbs2I) / (K2 + Pbs2I) - k4 * Sho1Pbs2A / (K4 + Sho1Pbs2A)
    dHog1IPbs2A = k5f * Hog1I * Sln1Pbs2A + k6f * Hog1I * Sho1Pbs2A - k5r * Hog1IPbs2A - k6r * Hog1IPbs2A - kcat * Hog1IPbs2A
    dHog1A = kcat * Hog1IPbs2A - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A - k10 * Glycerol

    return dSln1Pbs2A, dSho1Pbs2A, dHog1IPbs2A, dHog1A, dGlycerol