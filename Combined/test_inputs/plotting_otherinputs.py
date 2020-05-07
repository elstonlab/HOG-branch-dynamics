import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import seaborn as sns
from model_otherinputs import *
import pandas as pd
import warnings

phospho_palette = {150:'#8ace88', 350: '#319a50', 550:'#005723'}
nuc_palette = {150:'#84bcdb', 350: '#196789', 550:'#084082'}

palettes = {'phospho':phospho_palette,
           'nuc':nuc_palette,
           'nuc2':nuc_palette,
           'nuc3':nuc_palette,
           'nuc4':nuc_palette}

x_labels = {'phospho': 'pp Hog1',
          'nuc': 'nuc Hog1'}

param_colors = {'M16': ['#be202e','#606060', '#851747','#33669a','#05582d', '#851747','#33669a','#05582d', '#851747','#33669a','#05582d', '#851747','#33669a'],
          }

def plt_param_behaviors(model_fxns, top_params, plt_top, params_constants, initials,  doses, time, param,
                        mapk_wt_data=None, mapk_time=None, ss=False, plt_bad=0,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 20}
    plt.rc('font', **font)

    fig, (ax1,ax2,ax3,ax4) = plt.subplots(1, 4, figsize=(20,4))


    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    palette = palettes.get(param)

    ax1.set_yticks(np.arange(0, 101, step=25))
    ax1.set_xticks(np.arange(0, 61, step=15))

    if mapk_wt_data:
        params = top_params.mean(axis=0)
       # for sig, wt_data in zip(doses, mapk_wt_data):
       #     ax1.plot(mapk_time, wt_data, 'o', markersize=10, color=palette.get(sig), label = str(int(sig/1000))+'mM KCl')

    if param == 3:
        dt = 0.1
        steps = 2001
        time = np.linspace(0,dt*steps,steps)
    else:
        ax1.set_ylim(0,100)
           

    for sig in doses:
        for params in top_params[0:plt_top]:
            #                 if ss:
                #ss_data = run_ss(model_fxns.m, initials, params_constants, params)
#                     data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, params, time)
#                 else:
                data = simulate_wt_experiment(model_fxns.m, initials, params_constants, sig, params, time)
                phospho = (data[:,2]+data[:,3])/params_constants[2]*100
                nuc = (data[:,3]+data[:,4])/params_constants[2]*100
                Hog1_IC = (params_constants[2] - data[:,2] - data[:,3] - data[:,4])/params_constants[2]*100
                Hog1_AC = data[:,2]/params_constants[2]*100
                Hog1_AN = data[:,3]/params_constants[2]*100
                Hog1_IN = data[:,4]/params_constants[2]*100
                sln1 = data[:,0]/params_constants[0]*100
                sho1 = data[:,1]/params_constants[1]*100
                comp = data[:,3]/params_constants[2]*100
                if param == 'phospho':
                    ax1.plot(time, Hog1_IC, color=palette.get(sig))
                    ax2.plot(time,Hog1_AC, color=palette.get(sig))
                    ax3.plot(time,Hog1_AN, color=palette.get(sig))
                    ax4.plot(time,Hog1_IN, color=palette.get(sig))
                elif param == 'nuc':
                    ax1.plot(time, nuc, color=palette.get(sig))
                    ax1.set_ylim([0,100])
                    ax2.plot(time,data[:,5])
                else:
                    print('wrong param')              

    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=1)
       
    if plt_bad:
        plt_thresh_behavior(model_fxns, top_params, plt_bad, params_constants, initials,  doses, time, param, ax1, ax2)
        
    if save_fig:
        plt.savefig("C:/Users/sksuzuki/Documents/Research/figures/simulations/"+save_fig+".png",
        dpi=300,bbox_inches='tight')
    plt.show()

# def plt_ramp_behaviors(model_fxns, top_params, plt_top, params_constants, initials, time, param,
#                         ss = False, hog1_ramp_data=None, mapk_ramp_time=None,
#                         save_fig=''):
#     fig, (ax1) = plt.subplots(1, 1, figsize=(9,4))

#     ax1.plot(mapk_ramp_time, hog1_ramp_data[0], 'o', markersize=10, color='Black')

#     colors = sns.color_palette("bone")
#     pal2 = sns.set_palette(colors)
#     ax1.set_ylabel('% pp Hog1', fontsize=16)
#     ax1.set_xlabel('Time (min)', fontsize=16)
#     ax1.set_yticks(np.arange(0, 101, step=25))
#     ax1.set_xticks(np.arange(0, 61, step=15))

#     for params in top_params[:plt_top]:
# #         for sig in doses:
#             if ss:
#                 ss_data = run_ss(model_fxns.m, initials, params_constants, params)
#                 data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, 0, params, time, run_type=['ramp'])
#             else:
#                 data = simulate_wt_experiment(model_fxns.m, initials, params_constants, 0, params, time, run_type=['ramp'])
#             active = data[:,param]/params_constants[param]*100
#             ax1.plot(time, active)
#     if save_fig:
#         plt.savefig("C:/Users/sksuzuki/Documents/Research/figures/simulations/"+save_fig+".png",
#         dpi=300,bbox_inches='tight')
#     plt.show()

def fit_data_to_list_stairs(model_fxns, timepoints, levels, top_params, params_constants, initials, time, param,
                        ss = False):
        sims = []
        for idx, params in enumerate(top_params): 
            if ss:
                ss_data = run_ss_stairs(model_fxns.m, initials, params_constants, params)
                data = simulate_wt_experiment_stairs(model_fxns.m, timepoints, levels, ss_data, params_constants, params, time)
            else:
                data = simulate_wt_experiment_stairs(model_fxns.m, timepoints, levels, initials, params_constants, params, time)
            if param == 'phospho':
                species = (data[:,2]+data[:,3])/params_constants[2]*100
            elif param == 'nuc':
                species = (data[:,3]+data[:,4])/params_constants[2]*100
            else:
                print('no data')
            sims.append(species)
        return sims
    
def plt_param_cis_stairs(model_fxns, timepoints, levels, top_params, params_constants, initials,  time, param,
                        ss=False, ci= 95,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 16}
    plt.rc('font', **font)

    if save_fig:
        fig, (ax1) = plt.subplots(1, 1, figsize=(2.25,2))
    else:
        fig, (ax1) = plt.subplots(1, 1, figsize=(5,3))

    plt.rc('ytick', labelsize=16)
    plt.rc('xtick', labelsize=16)

    palette = palettes.get(param)


    # dashes = None
    # if t100a:
    dashes= (2,2)

    sims = fit_data_to_list_stairs(model_fxns, timepoints, levels, top_params, params_constants, initials, time, param, ss=ss)
    if param == 'nuc':
        ax1 = sns.tsplot(sims, time,  ci = ci, dashes = dashes, color='#084082')
    elif param == 'phospho':
        ax1 = sns.tsplot(sims, time,  ci = ci, dashes = dashes, color='#005723')
    else:
        print('define param')

    mstyles = {0: ['full'],
                1: ['^', 'o'],
                2: ['D', 'o', '^', 'v', 's', 'D', 'o'],
                3: ['full'],
                4: ['full']}
    fstyles = {0: ['full'],
                1: ['full', 'bottom'],
                2: ['none', 'none'],
                3: ['full'],
                4: ['full']}
    # if plt_bad:
        # plt_thresh_behavior(model_fxns, top_params, plt_bad, params_constants, initials,  doses, time, param, ax1, ax2)

#     if param == 3:
#         ax1.set_xticks(np.arange(0, 201, step=50))
#         ax1.set_yticks(np.arange(0, 5, step=1))
#         ax1.set_ylim(0,4.25)

#     elif doses == [0]:
#         ax1.set_ylim(-5,105)
#         ax1.set_yticks(np.arange(0, 101, step=25))
#         ax1.set_xticks(np.arange(0, 201, step=50))
#     else:
#         ax1.set_yticks(np.arange(0, 101, step=25))
#         ax1.set_xticks(np.arange(0, 61, step=15))
#         ax1.set_ylim(-5,105)
#         ax1.set_xlim(-2,62)

    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=.5)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.legend(bbox_to_anchor=[1, 0.5], loc='best')

    if save_fig:
        plt.savefig(save_fig+".pdf", dpi=300,bbox_inches='tight')

    plt.show()






