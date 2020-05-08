# 200306 -> cd to joe's longleaf account
#           which affects data dir and ea sims dir

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import pandas as pd
import pathlib
import collections
import h5py
import os
import sys

sys.path.insert(1, '../python_modules/')
import model as model
import model_supp as model_supp

def get_data():

    # cluster
    base_folder = '/pine/scr/j/m/jmcgirr/rozemarijn/Data/Prescaled_input/'

    # local
    # base_folder = 'C:/Users/sksuzuki/Downloads/Rozemarijn/Prescaled_input/'


    wt_phospho_folder = base_folder + 'WT_phospho'
    wt_nuc_folder = base_folder + '90percent/WT_nuc'

    phospho_time, wt_phospho_data = load_csv_data(wt_phospho_folder)
    nuc_time, wt_nuc_data = load_csv_data(wt_nuc_folder)

    data = [wt_phospho_data, wt_nuc_data]
    time = [phospho_time, nuc_time]
    return data, time

def load_csv_data(folder):
    data = []
    doses = []
    for csv in pathlib.Path(folder).glob('*.csv'):
        f_data = pd.read_csv(csv)
        # time = f_data['Time'].tolist()
        time = f_data.iloc[:,0]
        dose = f_data['Dose'][0]
        doses.append(dose)
        # f_data=f_data.set_index('Time')
        f_data = f_data.iloc[:,1:4].mean(axis=1)
        f_data = f_data.tolist()
        data.append(f_data)
    data = np.array(data)
    re_idx = sorted(range(len(doses)), key=lambda k: doses[k])
    data = data[re_idx]
    return time, list(data)


# def run_ss(m, inits, total_protein, learned_params):
#     ss = fsolve(m, inits, args=(0,total_protein, 0, learned_params))
#     return ss
#
# def simulate_wt_experiment(m, inits, total_protein, sig, learned_params, time, run_type=None):
#     odes = odeint(m, inits, time, args=(total_protein, sig, learned_params, run_type))
#     return odes

def calc_mse(model_fxns, theta, exp_data, exp_time, params_constants, initials):
    return model_supp.calc_sim_score(model_fxns, theta, exp_data, exp_time, params_constants, initials)
    # return sum(mse) ##AMY (insert your own error function)

# def draw_theta2():
#     return 10**(-4+(4-(-4))*np.random.random(17))#np.random.uniform(.0001,1000,17)

def draw_thetas(sorted_params):
    idx = np.random.choice(range(len(sorted_params)), 1)
    return sorted_params[idx][0]

def step_theta(theta):
    log_theta = np.log10(theta)
    theta_prime = np.concatenate([10**(np.random.uniform(x-.1,x+.1,1)) for x in log_theta], axis=0)
    return theta_prime

def run_schedule_i(prior_thetas, ei, num_theta_primes, model_fxns):
    thetas_ei = []
    mses_ei = []
    c = collections.Counter({'Pass': 0, 'Fail': 0})
    while len(thetas_ei) < num_theta_primes:
        theta = draw_thetas(prior_thetas)
        theta_prime = step_theta(theta)
        mse = calc_mse(model_fxns, theta_prime, exp_data, exp_time, params_constants, initials) ##AMY error fxn
        if mse < ei:
            c['Pass'] += 1
            thetas_ei.append(theta_prime)
            mses_ei.append(mse)
            if len(mses_ei) % int(num_theta_primes*.1) == 0:
                print(str(int(len(mses_ei)/num_theta_primes*100)) + "% complete.")
        else:
            c['Fail'] += 1
    c = np.array(list(c.values()))
    return np.asarray(mses_ei), np.asarray(thetas_ei), c

def check_dir_exist():
    stripped_name = strip_filename(save_filename)
    dir_to_use = os.getcwd() + '/' + stripped_name
    #check if dir exists and make
    if not os.path.isdir(dir_to_use):
        os.makedirs(dir_to_use)
        fn = dir_to_use + '/' + 'output.txt'
        file = open(fn, 'w')
        script_name = os.path.basename(__file__)#__file__)
        open_script = open(script_name, 'r')
        write_script = open_script.read()
        file.write(write_script)
        open_script.close()

        file.close()
    return dir_to_use, stripped_name

def get_filename(dir_to_use, stripped_name, val):
    filename_base = dir_to_use + '/' + stripped_name + '_'
    if val < 10:
        toret = '000' + str(val)
    elif 10 <= val < 100:
        toret = '00' + str(val)
    elif 100 <= val < 1000:
        toret = '0' + str(val)
    else:
        toret = str(val)
    return filename_base + toret + '.hdf5'

def strip_filename(fn):
    if '/' in fn:
        fn = fn.split('/')[-1]
    fn = fn.split('.')[0]
    return fn

def add_info(fn, gens, inds, mutationrate, crossoverrate):
    #get current date:
    cur_date = timeski.strftime('%y%m%d')
    # setup EA data:
    ea_data = str(gens) + 'g' + str(inds) + 'i' + str(int(mutationrate*100)) + 'm' + str(int(crossoverrate*100)) + 'c'
    #put it all together:
    #new_fn = cur_date + '_' + fn + '_' + ea_data
    new_fn = cur_date + '_' + os.path.basename(fn).split('.')[0].split('_')[-1] + '_' + ea_data
    return new_fn

def data_to_hdf5(dir_to_use, stripped_name, mses, thetas, c=None):
    # arr_to_hdf5 = [arr_best_score, arr_best_ind]
    counter = 0
    filename = get_filename(dir_to_use, stripped_name, counter)
    while os.path.isfile(filename) == True:
        counter += 1
        filename = get_filename(dir_to_use, stripped_name, counter)
    print(filename)
    with h5py.File(filename, 'w') as f:
        f.create_dataset("mses", data = mses)
        f.create_dataset("thetas", data = thetas)
        if not c is None:
            f.create_dataset("count", data=c)

def recalc_mses(model_fxns, EA_theta_set, exp_data, exp_time, params_constants, initials):
    mses = []
    for params in EA_theta_set:
        mses.append(model_supp.calc_sim_score(model_fxns, params, exp_data, exp_time, params_constants, initials))
    re_idx = sorted(range(len(mses)), key=lambda k: mses[k])
    thetas = EA_theta_set[re_idx]
    mses = np.sort(mses)
    return mses, thetas

def def_schedules(sorted_mses):
    best_mse = sorted_mses[0]
    worst_mse = sorted_mses[-1]

    # e1 = (best_mse+worst_mse)/2 #will take longer to run
    e1 = worst_mse
    e2 = (e1+best_mse)/2
    e3 = (e2+best_mse)/2
    e4 = (e3+best_mse)/2
    e5 = (e4+best_mse)/2
    e6 = best_mse
    return e1, e2, e3, e4, e5, e6

def get_ea_data(f, dir_to_use, stripped_name):
    all_params, last_params, all_mses, last_mses = model_supp.get_sim_data(f, num_sims=2000)
    mses_EA, thetas_EA = recalc_mses(model_fxn, last_params, exp_data, exp_time, params_constants, initials)
    data_to_hdf5(dir_to_use, stripped_name, mses_EA, thetas_EA)
    return mses_EA, thetas_EA

def main(f, number_eas, particle_num):
    dir_to_use, stripped_name = check_dir_exist()
    mses_EA, thetas_EA = get_ea_data(f, dir_to_use, stripped_name)
    # f = h5py.File('t190924_kb_M2_ea_abc_smc/t190924_kb_M2_ea_abc_smc_0000.hdf5','r')
    # mses_EA = f["mses"]
    # thetas_EA = f["thetas"]
    mses_EA = mses_EA[:number_eas]
    thetas_EA = thetas_EA[:number_eas]
    eis = def_schedules(mses_EA)
    print(eis)
    prior_thetas = thetas_EA
    for i, ei in enumerate(eis):
        print("Now running schedule #" + str(i))
        ABC_SMC_mses, ABC_SMC_thetas, ABC_SMC_c = run_schedule_i(prior_thetas, ei, particle_num, model_fxn)
        print(ABC_SMC_c)
        data_to_hdf5(dir_to_use, stripped_name, ABC_SMC_mses, ABC_SMC_thetas, ABC_SMC_c)
        prior_thetas = ABC_SMC_thetas


if __name__ == '__main__':
    exp_data, exp_time = get_data()
    Sln1_tot = model_supp.molarity_conversion(1176)
    Sho1_tot = model_supp.molarity_conversion(1534)
    Hog1_tot = model_supp.molarity_conversion(8225)
    params_constants = [Sln1_tot, Sho1_tot, Hog1_tot, 0] #mM

    Sln1 = 0
    Sho1 = 0
    Hog1_AC = 0
    Hog1_AN = 0
    Hog1_IN = 0.23 * Hog1_tot
    Glycerol = 0.0001
    initials = [Sln1, Sho1, Hog1_AC, Hog1_AN, Hog1_IN, Glycerol]

    model_fxn = model.Model(model.M28MMmin)
    save_filename = '2005_combined_90_smc.txt'
    number_eas = 400 #out of 500 based on mses from EA
    particle_num = 1000
    main('./200426_Model28_MM_90perc', number_eas, particle_num)
    # main('C:/Users/sksuzuki/Desktop/killdevil/runs_for_paper/190924_kb_M2/', number_eas, particle_num)
