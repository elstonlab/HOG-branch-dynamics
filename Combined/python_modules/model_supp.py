import pathlib
import pickle
import numpy as np
import pickle
import os
import pandas as pd
import model

def molarity_conversion(molecules):
    Na = 6.02214076*10**23
    cell_volume = 44                                                                 # volume of a yeast cell
    return molecules/(Na*cell_volume*10**-15)*1000                                # returns mM


def get_sim_data(folder, num_sims=None):
    try:
        my_abs_path = os.path.exists(folder)
    except FileNotFoundError:
        print("Folder not found ¯\_(ツ)_/¯")
    else:
        all_mses = []
        all_params = []
        empty_data = 0
        for i, loaded_data in enumerate(pathlib.Path(folder).glob('*.pickled')):
            if num_sims:
                if i > num_sims-1:
                    break
            if os.path.getsize(loaded_data) > 0:
            # with open(str(loaded_data), 'rb') as f: # python35
                with open(loaded_data, 'rb') as f:
                    new_data = pickle.load(f)
                    all_mses.append(np.asarray(new_data[0]))
                    all_params.append(np.asarray(new_data[1]))
            else:
                empty_data += 1

    print("Number of runs collected: " + str(len(all_params)))

    last_mses = [all_mses[i][len(all_mses[0])-1] for i in range(len(all_mses))]
    last_params = [all_params[i][len(all_params[0])-1] for i in range(len(all_params))]

    last_mses = np.array(last_mses)
    last_params = np.array(last_params)
    all_mses = np.array(all_mses)
    all_params = np.array(all_params)

    print('Best last gen MSE: ' + str(np.sort(last_mses)[0]))
    print('Mean last gen MSEs of top 5%: ' + str(np.mean(np.sort(last_mses)[:round(len(last_mses)*0.05)])))
    return all_params, last_params, all_mses, last_mses

def load_csv_data(folder):
    data = []
    doses = []
    for csv in pathlib.Path(folder).glob('*.csv'):
        f_data = pd.read_csv(csv)
        time = f_data['Time'].tolist()
        dose = f_data['Dose'][0]
        doses.append(dose)
        f_data=f_data.set_index('Time')
        f_data = f_data.iloc[:,:4].mean(axis=1)
        f_data = f_data.tolist()
        data.append(f_data)
    data = np.array(data)
    re_idx = sorted(range(len(doses)), key=lambda k: doses[k])
    data = data[re_idx]
    return time, list(data)
#
# def get_data(input):
#     if input == 'long':
#         base_folder = 'C:/Users/Rozemarijn/Documents/Universiteit/Masterstage2/Research/Modelling/Inputs/Input/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + 'WT_nuc'
#     elif input == 'short':
#         base_folder = 'C:/Users/Rozemarijn/Documents/Universiteit/Masterstage2/Research/Modelling/Inputs/Inputshort/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + 'WT_nuc'
#     elif input == 'norm_max':
#         base_folder = 'C:/Users/Rozemarijn/Documents/Universiteit/Masterstage2/Research/Modelling/Inputs/NormalizedInputs/Normalized_Input_max_mean/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + 'WT_nuc'
#     elif input == 'norm_ss':
#         base_folder = 'C:/Users/Rozemarijn/Documents/Universiteit/Masterstage2/Research/Modelling/Inputs/NormalizedInputs/Normalized_Input_ss_mean/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + 'WT_nuc'
#     elif input == '30perc':
#         base_folder = 'C:/Users/Rozemarijn/Documents/Universiteit/Masterstage2/Research/Modelling/Inputs/Prescaled_Input/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + '30percent/WT_nuc'
#     elif input == '60perc':
#         base_folder = 'C:/Users/Rozemarijn/Documents/Universiteit/Masterstage2/Research/Modelling/Inputs/Prescaled_Input/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + '60percent/WT_nuc'
#     elif input == '90perc':
#         base_folder = '../exp_data/Prescaled_Input/'
#         wt_phospho_folder = base_folder + 'WT_phospho'
#         wt_nuc_folder = base_folder + '90percent/WT_nuc'
#     else:
#         print('wrong input')
#
#
#     # load experimental data
#     phospho_time, wt_phospho_data = load_csv_data(wt_phospho_folder)
#     nuc_time, wt_nuc_data = load_csv_data(wt_nuc_folder)
#
#     data = [wt_phospho_data, wt_nuc_data]
#     time = [phospho_time, nuc_time]
#
#     return data, time


def calc_sim_score(model_fxns, params, data, exp_time, total_protein, inits):
    wt_phospho_data, wt_nuc_data = data
    phospho_time, nuc_time = exp_time

    hog1_doses = [150, 350,550]

    wt_ss_inits = model.run_ss(model_fxns.m, inits, total_protein, params)
    if (wt_ss_inits < 0).any():
        return ((63+69+18*2+27)*100)**2

    dt = 0.1
    steps = 601
    time = np.linspace(0,dt*steps,steps)

    closest_idxs_phospho = [np.abs(time - t).argmin() for t in phospho_time]
    closest_idxs_nuc = [np.abs(time-t).argmin() for t in nuc_time]

    # check if wt steady state exists, returns maximal MSE if not
    Hog1_ss = np.sum(wt_ss_inits[2:4])
    wt_ss_check = np.concatenate([wt_ss_inits[:2], [Hog1_ss]], axis=0)
    check = total_protein[:-1] - wt_ss_check
    if (check < 0).any():
        return ((63+69+18*2+27)*100)**2

    mse_total = 0

    # WILDTYPE
    for dose, data_phospho, data_nuc in zip(hog1_doses, wt_phospho_data, wt_nuc_data):
        odes = model.simulate_wt_experiment(model_fxns.m,wt_ss_inits, total_protein, dose, params, time)
        Hog1_phospho = (odes[:,2]+odes[:,3])/total_protein[2]*100 #calc as a percentage
        error_phospho = np.sum((data_phospho - Hog1_phospho[closest_idxs_phospho])**2) #sum of squares to keep the same
        mse_total += error_phospho
        Hog1_nuc = (odes[:,3]+odes[:,4])/total_protein[2]*100 # percentage
        error_nuc = np.sum((data_nuc- Hog1_nuc[closest_idxs_nuc])**2) #calc sum of squares since between 0 and 1
        mse_total += error_nuc
    return mse_total



def get_saved_thetas(f):
    # df = pd.read_csv(str(f))
    # df = df.drop(df.columns[0], axis=1)
    return np.array(pd.read_csv(f).drop(['Unnamed: 0'], axis=1))


def sort_mses_thetas(mses, thetas):
    re_idx = sorted(range(len(mses)), key=lambda k: mses[k])
    thetas = thetas[re_idx]
    return np.sort(mses), thetas
