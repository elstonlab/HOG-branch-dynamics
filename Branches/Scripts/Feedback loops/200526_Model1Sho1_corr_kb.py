# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:27:15 2020

@author: Rozemarijn
"""
###################################################################
#IMPORT PACKAGES
###################################################################
import numpy as np
from scipy.integrate import odeint
# from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from deap import base, creator, tools, algorithms
import os
import sys
import pickle
import time as timeski
import math
# from itertools import product
import pathlib
import pandas as pd
# import multiprocessing
# import sys
# import pathlib


###########################################################################
#LOAD EXPERIMENTAL DATA
###########################################################################
def load_csv_data(folder):
    data = []
    doses = []
    for csv in pathlib.Path(folder).glob('*.csv'):
        f_data = pd.read_csv(csv)
        time = f_data['Time'].tolist()
        dose = f_data['Dose'][0]
        doses.append(dose)
        f_data=f_data.set_index('Time')
        f_data = f_data.iloc[:,:3].mean(axis=1)
        f_data = f_data.tolist()
        data.append(f_data)
    data = np.array(data)
    re_idx = sorted(range(len(doses)), key=lambda k: doses[k])
    data = data[re_idx]
    return time, list(data)



#############################################################################
# Convert molecules to molar concentrations
#############################################################################
def molarity_conversion(molecules):
    Na = 6.02214076*10**23
    cell_volume = 44                                 # volume of a yeast cell
    return molecules/(Na*cell_volume*10**-15)*1000 # returns uM


###################################################################
#MATRIX FOR VARIABLES TO INTERP AND EXPONENTIATE
###################################################################

def make_conversion_matrix(number_of_params, maximums, minimums):
    # want easily savable matrix to hold this info
    # interp boolean, interp range (min,max), power boolean, power number (y)
    arr_IandP = np.zeros((5,number_of_params))
    # Set all interp booleans to 1 - everything is going to be interpreted
    arr_IandP[0,:] = 1
    # Set all power booleans to 1 - everything is in the form of powers
    arr_IandP[3,:] = 1
    # Set all power numbers to 10 - everything has a base of 10
    arr_IandP[4,:] = 10
    # Set minimums and maximums for all parameters. Parameters are in the following order:

    for i in range(len(minimums)):
        arr_IandP[1,i] = minimums[i] #interp_range_min
        arr_IandP[2,i] = maximums[i] #interp_range_max

    return arr_IandP

#############################################################################
#EA FUNCTIONS
#############################################################################

def convert_individual(ea_individual, conversion_matrix, number_of_params):
    # copy and get len of individual
    arr_params_conv = np.zeros(number_of_params)#np.copy(arr_parameters)
    len_ind = len(ea_individual)

    # Interp:
    for idx in np.nonzero(conversion_matrix[0])[0]:
        ea_val = ea_individual[idx]
        r_min = conversion_matrix[1][idx]
        r_max = conversion_matrix[2][idx]
        arr_params_conv[idx] = np.interp(ea_val, (0,1), (r_min, r_max))
        # print(ea_val)

    # Exponentiate:
    for idx in np.nonzero(conversion_matrix[3])[0]:
        ea_val = arr_params_conv[idx]
        base_val = conversion_matrix[4][idx]
        arr_params_conv[idx] = np.power(base_val, ea_val)

    # arr_params_conv[-4:] = np.round(arr_params_conv[-4:],0)

    return arr_params_conv


def scorefxn_helper(individual):
    # just a helper function that pulls all of scorefxn1 dependencies together
    # note the (), <--using single optimization in DEAP for now
    # scorefxn1 is taking care of the multiple optimizations for now
    return scorefxn1(inits, total_protein, individual),

###################################################################
#CHECK FOR / CREATE DIR FOR DATA
###################################################################
def strip_filename(fn):
    #input = full path filename
    #output = filename only
    #EX input = '/home/iammoresentient/phd_lab/data/data_posnegfb_3cellsum.pickled'
    #EX output = 'data_posnegfb_3cellsum'
    if '/' in fn:
        fn = fn.split('/')[-1]
    fn = fn.split('.')[0]
    return fn

def make_dirs(save_filename):
    stripped_name = strip_filename(save_filename)
    dir_to_use = os.getcwd() + '/' + stripped_name

    #check if dir exists and make
    if not os.path.isdir(dir_to_use):
        os.makedirs(dir_to_use)
        # print('Made: ' + dir_to_use)
        # and make README file:
        fn = dir_to_use + '/' + 'output.txt'
        file = open(fn, 'w')

        # write pertinent info at top
        file.write('OUTPUT\n\n')
        file.write('Filename: ' + stripped_name + '\n')
        file.write('Directory: ' + dir_to_use + '\n')
        file.write('Data file: ' + save_filename + '\n\n')
        file.write('\n\n\n\n')

        #write script to file
        #script_name = os.getcwd() + '/' + 'EA_1nf1pf.py'
        script_name = os.path.basename(__file__)#__file__)
        open_script = open(script_name, 'r')
        write_script = open_script.read()
        file.write(write_script)
        open_script.close()

        file.close()
    return stripped_name, dir_to_use

###################################################################
#LOOP: EVOLUTIONARY ALGORITHM + SAVE DATA
###################################################################
def run():


    for i in range(number_of_runs):
        ###################################################################
        #EVOLUTIONARY ALGORITHM
        ###################################################################
        #TYPE
        #Create minimizing fitness class w/ single objective:
        creator.create('FitnessMin', base.Fitness, weights=(-1.0,))
        #Create individual class:
        creator.create('Individual', list, fitness=creator.FitnessMin)

        #TOOLBOX
        toolbox = base.Toolbox()
        #Register function to create a number in the interval [1-100?]:
        #toolbox.register('init_params', )
        #Register function to use initRepeat to fill individual w/ n calls to rand_num:
        toolbox.register('individual', tools.initRepeat, creator.Individual,
                         np.random.random, n=number_of_params)
        #Register function to use initRepeat to fill population with individuals:
        toolbox.register('population', tools.initRepeat, list, toolbox.individual)

        #GENETIC OPERATORS:
        # Register evaluate fxn = evaluation function, individual to evaluate given later
        toolbox.register('evaluate', scorefxn_helper)
        # Register mate fxn = two points crossover function
        toolbox.register('mate', tools.cxTwoPoint)
        # Register mutate by swapping two points of the individual:
        toolbox.register('mutate', tools.mutPolynomialBounded,
                         eta=0.1, low=0.0, up=1.0, indpb=0.2)
        # Register select = size of tournament set to 3
        toolbox.register('select', tools.selTournament, tournsize=3)

        #EVOLUTION!
        pop = toolbox.population(n=number_of_individuals)
        hof = tools.HallOfFame(1)

        stats = tools.Statistics(key = lambda ind: [ind.fitness.values, ind])
        stats.register('all', np.copy)

        # using built in eaSimple algo
        pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=crossover_rate,
                                           mutpb=mutation_rate,
                                           ngen=number_of_generations,
                                           stats=stats, halloffame=hof,
                                           verbose=False)
        # print(f'Run number completed: {i}')

        ###################################################################
        #MAKE LISTS
        ###################################################################
        # Find best scores and individuals in population
        arr_best_score = []
        arr_best_ind = []
        for a in range(len(logbook)):
            scores = []
            for b in range(len(logbook[a]['all'])):
                scores.append(logbook[a]['all'][b][0][0])
            #print(a, np.nanmin(scores), np.nanargmin(scores))
            arr_best_score.append(np.nanmin(scores))
            #logbook is of type 'deap.creator.Individual' and must be loaded later
            #don't want to have to load it to view data everytime, thus numpy
            ind_np = np.asarray(logbook[a]['all'][np.nanargmin(scores)][1])
            ind_np_conv = convert_individual(ind_np, arr_conversion_matrix, number_of_params)
            arr_best_ind.append(ind_np_conv)
            #arr_best_ind.append(np.asarray(logbook[a]['all'][np.nanargmin(scores)][1]))


        # print('Best individual is:\n %s\nwith fitness: %s' %(arr_best_ind[-1],arr_best_score[-1]))

        ###################################################################
        #PICKLE
        ###################################################################
        arr_to_pickle = [arr_best_score, arr_best_ind]

        def get_filename(val):
            filename_base = dir_to_use + '/' + stripped_name + '_'
            if val < 10:
                toret = '000' + str(val)
            elif 10 <= val < 100:
                toret = '00' + str(val)
            elif 100 <= val < 1000:
                toret = '0' + str(val)
            else:
                toret = str(val)
            return filename_base + toret + '.pickled'

        counter = 0
        filename = get_filename(counter)
        while os.path.isfile(filename) == True:
            counter += 1
            filename = get_filename(counter)

        pickle.dump(arr_to_pickle, open(filename,'wb'))


#############################################################################
# MODEL
#############################################################################


def model(initials,t,total_protein,sig,params):
    Sln1, Sho1, Hog1A, Glycerol = initials
    Sln1_tot, Sho1_tot, Hog1_tot, _ = total_protein
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = params #18

    Hog1I = Hog1_tot - Hog1A
    Sln1_inactive = Sln1_tot - Sln1
    Sho1_inactive = Sho1_tot - Sho1
    
    if Glycerol < 0:
        Glycerol = 0

    dSln1 = (base_osmo + k1 * sig - Glycerol) * (Sln1_inactive) / (K1 + Sln1_inactive) - k3 * Sln1 / (K3 + Sln1)
    dSho1 = ((k2 * sig - Glycerol)) * (Sho1_inactive) / (K2 + Sho1_inactive) - k4 * Sho1 / (K4 + Sho1)
    dHog1A = (k5 * Sln1 + k6 * Sho1) * Hog1I / (K56 + Hog1I) - k7 * Hog1A / (K7 + Hog1A)
    dGlycerol = k8 * Hog1A + k9 * (base_osmo + sig - Glycerol) - k10 * Glycerol

    return dSln1, dSho1, dHog1A, dGlycerol

def run_ss(inits, total_protein, learned_params):
    ss = fsolve(model, inits, args=(0,total_protein, 0, learned_params))
    return ss

def simulate_wt_experiment(inits, total_protein, sig, learned_params, time):
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def simulate_sln1_experiment(inits, total_protein, sig, learned_params, time):
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = learned_params #16
    learned_params = base_osmo, k1, K1, 0, 0, k3, K3, 0, 0, k5, 0, K56, k7, K7, k8, k9, k10
    #solve odes:
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

def simulate_sho1_experiment(inits, total_protein, sig, learned_params, time):
    base_osmo, k1, K1, k2, K2, k3, K3, k4, K4, k5, k6, K56, k7, K7, k8, k9, k10 = learned_params #16
    learned_params = base_osmo, 0, 0, k2, K2, 0, 0, k4, K4, 0, k6, K56, k7, K7, k8, k9, k10
    #solve odes:
    odes = odeint(model, inits, time, args=(total_protein, sig, learned_params))
    return odes

#############################################################################
# SCOREFXN
#############################################################################

def scorefxn1(inits, total_protein, learned_params):
    params = convert_individual(learned_params, arr_conversion_matrix, number_of_params)

    dt = 0.1
    steps = 601
    time = np.linspace(0,dt*steps,steps)

    closest_idxs_wt = [np.abs(time - t).argmin() for t in wt_time]
    closest_idxs_sln1 = [np.abs(time - t).argmin() for t in sln1_time]
    closest_idxs_sho1 = [np.abs(time - t).argmin() for t in sho1_time]

    # check if wt steady state exists, returns maximal MSE if not
    wt_ss_inits = run_ss(inits, total_protein, params)
    check = total_protein[:-1] - wt_ss_inits[:-1]
    if (check < 0).any():
        return ((9*5)*100)**2 #if sims were the oppsite of the data (not possible)

    mse_total = 0

    # WILDTYPE
    for dose, data_wt in zip(hog1_doses, wt_data):
        odes = simulate_wt_experiment(wt_ss_inits, total_protein, dose, params, time)
        Hog1_wt = odes[:,2]/total_protein[2]*100
        error_wt = ((data_wt - Hog1_wt[closest_idxs_wt])**2).mean()
        mse_total += error_wt

    # Sln1
    for dose, data_sln1 in zip(hog1_doses, sln1_data):
        odes = simulate_sln1_experiment(wt_ss_inits, total_protein, dose, params, time)
        Hog1_sln1 = odes[:,2]/total_protein[2]*100
        error_sln1 = ((data_sln1 - Hog1_sln1[closest_idxs_sln1])**2).mean()
        mse_total += error_sln1

    # Sho1
    for dose, data_sho1 in zip(hog1_doses, sho1_data):
        odes = simulate_sho1_experiment(wt_ss_inits, total_protein, dose, params, time)
        Hog1_sho1 = odes[:,2]/total_protein[2]*100
        error_sho1 = ((data_sho1 - Hog1_sho1[closest_idxs_sho1])**2).mean()
        mse_total += error_sho1

    return mse_total


if __name__ == '__main__':

    #base filename
    save_filename = '200526_Model1Sho1_kb.txt'

    # Paths to longleaf data
    base_folder = '/nas/longleaf/home/sksuzuki/rvdv/Sho1_branch/'
    # Paths to local data
    # base_folder = 'C:/Users/sksuzuki/Downloads/Rozemarijn/Input/'
    wt_folder = base_folder + 'WT_phospho'
    sln1_folder = base_folder + 'Sln1_phopho'
    sho1_folder = base_folder + 'Sho1_phospho'

    # load experimental data
    wt_time, wt_data = load_csv_data(wt_folder)
    sln1_time, sln1_data = load_csv_data(sln1_folder)
    sho1_time, sho1_data = load_csv_data(sho1_folder)


    # Protein concentrations (mM)
    Sln1_tot = molarity_conversion(1176)
    Sho1_tot = molarity_conversion(1534)
    Hog1_tot = molarity_conversion(8225)
    total_protein = [Sln1_tot, Sho1_tot, Hog1_tot, 0] #mM

    # initial values
    Sln1 = 0
    Sho1 = 0
    Hog1A = 0.01 * Hog1_tot
    Glycerol = 0.0001
    inits = [Sln1, Sho1, Hog1A, Glycerol]

    # doses
    hog1_doses = [150, 350, 550]

    # Parameter ranges
    number_of_params = 17
    minimums = [1, -4, -4, -4, -4, -4,
                -4, -4, -4, -4, -4, -4,
                -4, -4, -4, -4, -4
        ]

    maximums = [2, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4
        ]

    # EA params
    number_of_runs = 5
    number_of_generations = 200
    number_of_individuals = 100
    mutation_rate = 0.2
    crossover_rate = 0.5

    stripped_name, dir_to_use = make_dirs(save_filename)
    arr_conversion_matrix = make_conversion_matrix(number_of_params, maximums, minimums)

    run()
