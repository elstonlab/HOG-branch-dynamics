import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import seaborn as sns
from model import *
import pandas as pd
import warnings

phospho_palette = {150:'#8ace88', 350: '#319a50', 550:'#005723'}
nuc_palette = {150:'#84bcdb', 350: '#196789', 550:'#084082'}

palettes = {'phospho':phospho_palette,
           'nuc':phospho_palette,
           'nuc2':nuc_palette,
           'nuc3':nuc_palette,
           'nuc4':nuc_palette,
           'sho1':nuc_palette,
           'sln1':nuc_palette,
           'all Hog1': phospho_palette}

x_labels = {'phospho': 'pp Hog1',
          'nuc': 'nuc Hog1'}

param_colors = {'M16': ['#be202e','#606060', '#851747','#33669a','#05582d', '#851747','#33669a','#05582d', '#851747','#33669a','#05582d', '#851747','#33669a'],
          }

def plt_param_integrals(model_fxns, top_params, plt_top, params_constants, initials,  doses, time, param,
                        mapk_wt_data=None, mapk_time=None, ss=False, plt_bad=0,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)

    fig, (ax1) = plt.subplots(1, 1, figsize=(3,2.5))


    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    palette = palettes.get(param)

    # ax1.set_yticks(np.arange(0, 101, step=25))
    ax1.set_xticks(np.arange(0, 61, step=15))

    if mapk_wt_data:
        params = top_params.mean(axis=0)
        for sig, wt_data in zip(doses, mapk_wt_data):
            ax1.plot(mapk_time, wt_data, 'o', markersize=10, color=palette.get(sig), label = str(int(sig/1000))+'mM KCl')

    # if param == 3:
    #     dt = 0.1
    #     steps = 2001
    #     time = np.linspace(0,dt*steps,steps)
    # else:
    #     ax1.set_ylim(0,100)


    for sig in doses:
        if param == 'all Hog1':
            ss_data = run_ss(model_fxns.m, initials, params_constants, top_params[0:plt_top][-1])
            data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, top_params[0:plt_top][-1], time)
            labels = ['Cytosolic', 'Cytosolic', 'Nuclear', 'Nuclear']
            areas1=np.zeros(601)
            areas2=np.zeros(601)
            for i in range(1,601):
                areas1[i-1]=(data[i-1,2]+data[i,2])*10/params_constants[2]
                areas2[i-1]=(data[i-1,3]+data[i,3])*10/params_constants[2]
            integral1=np.zeros((600))
            integral2=np.zeros((600))
            integral1[0]=areas1[0]
            integral2[0]=areas2[0]
            for j in range(1,600):
                integral1[j]=areas1[j]+integral1[j-1]
                integral2[j]=areas2[j]+integral2[j-1]
            ax1.plot(time[1:], integral1, label = labels[1], color = '#228833', linewidth=2)
            ax1.plot(time[1:], integral2,  label = labels[2], color = '#4477AA', linewidth=2)
            ax1.legend(title='Active',loc='best')
            
        for params in top_params[0:plt_top]:
            #                 if ss:
                ss_data = run_ss(model_fxns.m, initials, params_constants, params)
#                     data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, params, time)
#                 else:
                data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, params, time)
                phospho = (data[:,2]+data[:,3])/params_constants[2]*100
                nuc = (data[:,3]+data[:,4])/params_constants[2]*100
                if param == 'phospho':
                    ax1.plot(time, phospho, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                elif param == 'nuc':
                    ax1.plot(time, nuc, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                    # ax2.plot(time,data[:,5])
                elif param == 'sln1':
                    ax1.plot(time, data[:,0]/params_constants[0]*100, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                elif param == 'sho1':
                    ax1.plot(time, data[:,1]/params_constants[1]*100, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                elif param == 'all Hog1':
                    # labels = ['Inactive Cytosolic', 'Active Cytosolic', 'Active Nuclear', 'Inactive Nuclear']
                    # labels = ['Cytosolic', 'Cytosolic', 'Nuclear', 'Nuclear']
                    ax1.plot(time[1:], integral1, label = labels[1], color = '#228833', linewidth=2)
                    ax1.plot(time[1:], integral2,  label = labels[2], color = '#4477AA', linewidth=2)
                    # fig.legend([ax1,ax2], labels, loc='center left', bbox_to_anchor=(1, 0.5))
                    #ax1.set_ylim([-2,102])
                else:
                    print('wrong param')

    # ax1.legend(title='Active',loc='best')
    # ax2.legend(title='Inactive', loc='best')
    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=1)


    if plt_bad:
        plt_thresh_behavior(model_fxns, top_params, plt_bad, params_constants, initials,  doses, time, param, ax1, ax2)

    if save_fig:
        plt.savefig("C:/Users/sksuzuki/Documents/Research/figures/simulations/"+save_fig+".png",
        dpi=300,bbox_inches='tight')
    plt.show()

def plt_param_behaviors(model_fxns, top_params, plt_top, params_constants, initials,  doses, time, param,
                        mapk_wt_data=None, mapk_time=None, ss=False, plt_bad=0,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)

    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(6,2.5))


    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    palette = palettes.get(param)

    # ax1.set_yticks(np.arange(0, 101, step=25))
    ax1.set_xticks(np.arange(0, 61, step=15))

    if mapk_wt_data:
        params = top_params.mean(axis=0)
        for sig, wt_data in zip(doses, mapk_wt_data):
            ax1.plot(mapk_time, wt_data, 'o', markersize=10, color=palette.get(sig), label = str(int(sig/1000))+'mM KCl')

    # if param == 3:
    #     dt = 0.1
    #     steps = 2001
    #     time = np.linspace(0,dt*steps,steps)
    # else:
    #     ax1.set_ylim(0,100)


    for sig in doses:
        if param == 'all Hog1':
            ss_data = run_ss(model_fxns.m, initials, params_constants, top_params[0:plt_top][-1])
            data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, top_params[0:plt_top][-1], time)
            labels = ['Cytosolic', 'Cytosolic', 'Nuclear', 'Nuclear']
            ax2.plot(time, (params_constants[2]-data[:,2]-data[:,3]-data[:,4])/params_constants[2]*100, label = labels[0], color = '#CCBB44', linewidth=2)
            ax1.plot(time, data[:,2]/params_constants[2]*100, label = labels[1], color = '#228833', linewidth=2)
            ax1.plot(time, data[:,3]/params_constants[2]*100,  label = labels[2], color = '#4477AA', linewidth=2)
            ax2.plot(time, data[:,4]/params_constants[2]*100,  label = labels[3], color = '#66CCEE', linewidth=2)
            ax1.legend(title='Active',loc='best')
            ax2.legend(title='Inactive', loc='center left', bbox_to_anchor=(1, 0.5))
        for params in top_params[0:plt_top]:
            #                 if ss:
                ss_data = run_ss(model_fxns.m, initials, params_constants, params)
#                     data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, params, time)
#                 else:
                data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, sig, params, time)
                phospho = (data[:,2]+data[:,3])/params_constants[2]*100
                nuc = (data[:,3]+data[:,4])/params_constants[2]*100
                if param == 'phospho':
                    ax1.plot(time, phospho, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                elif param == 'nuc':
                    ax1.plot(time, nuc, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                    # ax2.plot(time,data[:,5])
                elif param == 'sln1':
                    ax1.plot(time, data[:,0]/params_constants[0]*100, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                elif param == 'sho1':
                    ax1.plot(time, data[:,1]/params_constants[1]*100, color=palette.get(sig))
                    ax1.set_ylim([-2,102])
                elif param == 'all Hog1':
                    # labels = ['Inactive Cytosolic', 'Active Cytosolic', 'Active Nuclear', 'Inactive Nuclear']
                    # labels = ['Cytosolic', 'Cytosolic', 'Nuclear', 'Nuclear']
                    ax2.plot(time, (params_constants[2]-data[:,2]-data[:,3]-data[:,4])/params_constants[2]*100, label = labels[0], color = '#CCBB44', linewidth=2)
                    ax1.plot(time, data[:,2]/params_constants[2]*100, label = labels[1], color = '#228833', linewidth=2)
                    ax1.plot(time, data[:,3]/params_constants[2]*100,  label = labels[2], color = '#4477AA', linewidth=2)
                    ax2.plot(time, data[:,4]/params_constants[2]*100,  label = labels[3], color = '#66CCEE', linewidth=2)

                    # fig.legend([ax1,ax2], labels, loc='center left', bbox_to_anchor=(1, 0.5))
                    ax1.set_ylim([-2,102])
                    ax2.set_ylim([-2,102])
                else:
                    print('wrong param')

    # ax1.legend(title='Active',loc='best')
    # ax2.legend(title='Inactive', loc='best')
    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=1)
    ax2.grid(color='grey', linestyle='-', axis='y', linewidth=1)


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




def plt_mses_gen(gen, mses_all, idx_top, save_fig=''):
    plt.clf()
    fig, (ax3) = plt.subplots(1, 1, figsize=(9,4))
    colors2 = sns.color_palette("Greys", 20)[10:]
    pal2 = sns.set_palette(colors2)
    ax3.set_xlabel('Generation', fontsize=20)
    mses_all_sort = mses_all[mses_all[:,-1].argsort()]
    for mses in mses_all_sort[:idx_top]:
        # ax3.semilogy([x for x in range(gen)], mses[:gen]) # for plotting on log scale
        ax3.plot([x for x in range(gen)], mses[:gen])
    ax3.yaxis.grid(True)
    ax3.set_ylabel('MSE', fontsize=20)
    ax3.set_xlim([0,gen])
    # ax3.set_ylim(1000,5000)

    if save_fig:
        plt.savefig(save_fig+".png",
        dpi=300,bbox_inches='tight')
    plt.show()

def plt_param_ranges(labelnames, m_name, dims, param_data, single_theta=pd.Series(),num=0, save_fig=''):
    # fig, (ax1) = plt.subplots(1, 1, figsize=(3,3))
    fig, (ax1) = plt.subplots(1, 1, figsize=(6,3))

    pal = sns.set_palette(param_colors.get(m_name))

    # Hide the right and top spiness
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    major_ticks = np.arange(-4, 9, 2)
    ax1.set_yticks(major_ticks)

    with sns.axes_style("whitegrid"):
        plt.bar(range(0,len(labelnames)),height=dims[0],bottom=dims[1],align='center',tick_label=labelnames, color='#dcdcdc',alpha = 0.8)
        ax1 = sns.swarmplot(x='param',y='vals', data = param_data, size=5) #size 3
        ax1.set_xticklabels(labelnames,rotation=90)
        plt.xlabel('Parameters', fontsize=20, fontweight='medium')
        ax1.set_ylabel('')

    plt.grid(color='#606060', which='major', axis='y', linestyle='solid')
    if single_theta.any:
        single_theta = pd.DataFrame(single_theta.loc[num])
        single_theta['param']=single_theta.index
        single_theta.melt(var_name='param', value_name='vals')
        ax1 = sns.swarmplot(x='param', y=num, data = single_theta, color = 'black', size=8)
    # a=ax1.get_yticks().tolist()
    # y_names=['$10^{-6}$','$10^{-4}$','$10^{-2}$','$1$','$10^{2}$','$10^{4}$','$10^{6}$',]
    y_names=['$10^{-4}$','$10^{-2}$','$1$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$']
    ax1.set_yticklabels(y_names)

    if save_fig:
        plt.savefig(save_fig,
        dpi=300,bbox_inches='tight')
    plt.show()

def plt_mses_swarm(mses, save_fig=''): # hog1=True, t100a=False, pbs2=False, pbs2_t100a=False, ramp=False, ptp23D=False,
    mses = mses.apply(np.log10).melt(var_name='Dataset', value_name='MSEs')
    # mses = mses
    fig, (ax1) = plt.subplots(1, 1, figsize=(9,4))
    sns.swarmplot(x='Dataset', y="MSEs",
                            data=mses, palette="muted") #
    # ax1.set_xticklabels(['WT Hog1', 'Hog1-as Hog1', 'WT Pbs2', 'Hog1-as Pbs2', 'Ramp', 'Ptp2/3$\Delta$'],rotation=45)
    plt.show()

def plt_mses(mses, datasets, size, ptpD=True,save_fig=''): # hog1=True, t100a=False, pbs2=False, pbs2_t100a=False, ramp=False, ptp23D=False,
    # mses = mses.apply(np.log10).melt(var_name='Dataset', value_name='$\Sigma$ MSEs')
    hog1_doses = [0, 50000,150000,250000,350000,450000,550000]
    ptp_doses = [0, 150000,550000]
    mses_means = mses.mean()
    r1 = 'Hog1 WT'

    barWidth = .8
    r2 = 'Hog1-as'
    r3 = 'Pbs2 WT'
    r4 = 'Pbs2-as'
    r5 = 'Ramp'
    # r6 = 'Ptp2/3$\Delta$'

    fig, ax1 = plt.subplots(1, 1, figsize=size)


    if datasets[0]:
        ax1.bar(r1, mses_means[0], width=barWidth, edgecolor='white', color=MAPK_palette.get(hog1_doses[0]))#, label='Hog1'
        for i, dose in zip(range(1,7),hog1_doses[1:]) :
            ax1.bar(r1, mses_means[i], bottom=sum(mses_means[0:i]),  width=barWidth, edgecolor='white', color=MAPK_palette.get(dose))

    if datasets[1]:
        ax1.bar(r2, mses_means[7], width=barWidth, edgecolor='white', color=MAPK_palette.get(hog1_doses[0]))
        for i, dose in zip(range(1,7),hog1_doses[1:]) :
            ax1.bar(r2, mses_means[i+7], bottom=sum(mses_means[7:i+7]),  width=barWidth, edgecolor='white', color=MAPK_palette.get(dose))

    if datasets[2]:
        ax1.bar(r3, mses_means[14], width=barWidth, edgecolor='white', color=MAP2K_palette.get(150000))#, label='Hog1'
        ax1.bar(r3, mses_means[15], bottom=mses_means[14],  width=barWidth, edgecolor='white', color=MAP2K_palette.get(550000))

    if datasets[3]:
        ax1.bar(r4, mses_means[16], width=barWidth, edgecolor='white', color=MAP2K_palette.get(150000))#, label='Hog1'
        ax1.bar(r4, mses_means[17], bottom=mses_means[16],  width=barWidth, edgecolor='white', color=MAP2K_palette.get(550000))

    if datasets[4]:
        ax1.bar(r5, mses_means[18], width=barWidth, edgecolor='white', color='#008080')#, label='Hog1'

    if ptpD:
        if datasets[5]:
            ax1.bar(r6, mses_means[19], width=barWidth, edgecolor='white', color=ptp_palette.get(ptp_doses[0]))
            for i, dose in zip(range(20,22), ptp_doses[1:]):
                ax1.bar(r6, mses_means[i], bottom=mses_means[i-1],  width=barWidth, edgecolor='white',color=ptp_palette.get(dose))
        # ax1.set_xticklabels([r6],rotation=45)
        ax1.set_xticklabels([r1,r2,r3,r4,r5],rotation=45)



    # if predicted:
    #     plt.yticks([])
    #     ax1.set_yscale("log", nonposy='clip')
    plt.tight_layout()


    # plt.savefig(figure_dir+'fig3_D'+".pdf",dpi=200)
    # ax1.legend()
    plt.ylabel("MSE")
    # plt.xticklabels(labelnames,rotation=45)
    print(sum(mses_means[:-1]))
    if save_fig:
        plt.savefig("C:/Users/sksuzuki/Documents/Research/figures/simulations/"+save_fig+".png",
        dpi=300,bbox_inches='tight')
    plt.show()



def plt_corr(labelnames, df_top_params, save_fig=''):
    df_top_params.columns = labelnames
    df_top_params.head()
    cmap = sns.cm.rocket

    fig, (ax1) = plt.subplots(1, 1, figsize=(8,6))

    ax1 = sns.heatmap(df_top_params.corr(), cmap = cmap,  vmin=-.1, vmax=1.1) #annot=True,  fmt='d'
    plt.yticks(rotation=0)
    ax1.set_xticklabels(labelnames,rotation=90)
    ax1.tick_params(labelsize=12)
    # ax1.set_xlabel('Parameters', fontweight='medium', fontsize=18)
    # sns.set(font_scale=1.4)
    if save_fig:
        plt.savefig("C:/Users/sksuzuki/Documents/Research/figures/simulations/"+save_fig+".png",
            dpi=300,bbox_inches='tight')
    s = df_top_params.corr().unstack()
    so = s.sort_values(kind="quicksort")
    print(so[-10-len(labelnames):-len(labelnames)])
    plt.show()


def plt_rand_behaviors(model_fxns, top_params, plt_top, params_constants, initials, time, ramp_vals, average, selection='', save_fig=''):
    fig, (ax1) = plt.subplots(1, 1, figsize=(9,4))

    colors = sns.color_palette("Set2", 10)
    pal2 = sns.set_palette(colors)
    ax1.set_ylabel('% pp Hog1', fontsize=16)
    ax1.set_xlabel('Time (min)', fontsize=16)

    for params in top_params[:plt_top]:
        ss_data = run_ss(model_fxns.m, initials, params_constants, params)
        if selection == 't100a':
            data = model_fxns.t100a(model_fxns.m, ss_data, params_constants, 0, params, time, run_type=['rand', ramp_vals[1]])
        elif selection == 'nopos':
            data = model_fxns.nopos(model_fxns.m, ss_data, params_constants, 0, params, time, run_type=['rand', ramp_vals[1]])
        else:
            data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, 0, params, time, run_type=['rand', ramp_vals[1]])
        active = data[:,2]/params_constants[2]*100
        ax1.plot(time, active)
    ax1.plot(time, average, linewidth=4, color='black')
    if save_fig:
        plt.savefig("C:/Users/sksuzuki/Documents/Research/figures/simulations/"+save_fig+".png",
        dpi=300,bbox_inches='tight')
    ax1.set_ylim(0,100)
    plt.show()


def plt_deriv(mses, zoom, elim, deriv=1, window=False):
    a1 = np.sort(mses)[:elim]
    if window:
        a2 = np.concatenate([a1[window:], a1[-window:]])
    else:
        a2 = np.concatenate([a1[1:], [a1[-1]]])
    z = abs(a2-a1)
    if deriv == 2:
       b2 = np.concatenate([z[1:], [z[-1]]])
       z = b2-z
       # print(z)
    thresh = max(z)*0.01
    idx_thresh = np.argwhere(z<thresh)[0]
    fig, (ax1) = plt.subplots(1, 1, figsize=(9,4))
    plt.axhline(y=thresh, color='black', linewidth=3, label='1% of max')
    # plt.axvline(x=idx_thresh, color='black', linewidth=3)
    ax1.plot(range(zoom),z[:zoom], linewidth=3, color='grey')
    # ax1.plot(range(len(z)),z, 'o')
    ax1.set_xlabel('Index')
    ax1.set_ylabel('First derivative, win='+str(window))
    plt.legend()
    return thresh, idx_thresh

def plt_idx_vs_mse(mses, zoom, idx_thresh=False):
    mses = np.sort(mses)[:zoom]
    fig, ax1 = plt.subplots(1, 1, figsize=(9,4))
    plt.plot(range(len(mses)),mses, 'o', markersize=3, color='grey')
    plt.xlabel('Index')
    # plt.ylabel('$\sum MSE$')
    if idx_thresh:
        plt.axvline(x=idx_thresh, color='black', linewidth=3, label='1% max slope of the derivative')
    plt.legend()
    # ax1.semilogy(range(len(mses)),mses)


def fit_data_to_list(model_fxns, top_params, params_constants, initials, time, param, dose, ptpD=False,
                        ss = False):
        if param == 'all Hog1':
            IC = []
            AC = []
            AN = []
            IN = []
            for params in top_params:
                ss_data = run_ss(model_fxns.m, initials, params_constants, params)
                data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, dose, params, time)
                IC.append((params_constants[2]-data[:,2]-data[:,3]-data[:,4])/params_constants[2]*100)
                AC.append(data[:,2]/params_constants[2]*100)
                AN.append(data[:,3]/params_constants[2]*100)
                IN.append(data[:,4]/params_constants[2]*100)
            return [IC, AC, AN, IN]
        else:
            sims = []
            for idx, params in enumerate(top_params):
                if ss:
                    ss_data = run_ss(model_fxns.m, initials, params_constants, params)
                    data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, dose, params, time)
                else:
                    data = simulate_wt_experiment(model_fxns.m, initials, params_constants, dose, params, time)
                if param == 'phospho':
                    species = (data[:,2]+data[:,3])/params_constants[2]*100
                elif param == 'nuc':
                    species = (data[:,3]+data[:,4])/params_constants[2]*100
                else:
                    print('no data')
                sims.append(species)
            print('Dose: ' + str(dose) + ' complete.')
            return sims

def plt_param_cis(model_fxns, top_params, params_constants, initials,  doses, time, param,
                        exp_data=None, exp_time=None, ss=False, ptpD=False, ci= 95,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
    plt.rc('font', **font)

    if save_fig:
        fig, (ax1) = plt.subplots(1, 1, figsize=(4,2))
    else:
        fig, (ax1) = plt.subplots(1, 1, figsize=(5,3))

    plt.rc('ytick', labelsize=12)
    plt.rc('xtick', labelsize=12)

    palette = palettes.get(param)


    # dashes = None
    # if t100a:
    dashes= (2,2)

    for sig in doses:
        sims = fit_data_to_list(model_fxns, top_params, params_constants, initials, time, param, sig, ptpD=False, ss=ss)
        ax1 = sns.tsplot(sims, time,  ci = ci, dashes = dashes, color=palette.get(sig))

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
    if exp_data:
        params = top_params.mean(axis=0)
        for sig, data in zip(doses, exp_data[::-1]):
                if sig == 0:
                    mark = 'o'
                elif sig == 150:
                    mark = '^'
                elif sig == 550:
                    mark = 's'
                else:
                    mark = 'o'
                ax1.plot(exp_time, data, marker=mark, markersize=6,  linestyle="None", color=palette.get(sig), fillstyle='full', mec='black', label = str(int(sig)))#, fillstyle=fill, linestyle="None"), palette.get(sig), mec='black' (outside edge), linestyle="-",
                # if param == 'nuc':
                    # ax1.set_ylim([0,100])
                    # ax1.set_xlim([-2, 61])
                # else:
                    # ax1.set_ylim([0,100])


    ax1.set_ylim(21,65)
    # ax1.set_ylim(-5,105)

    # ax1.set_xlim(-1,31)
    ax1.set_xlim(-2,61)

    ax1.set_yticks(np.arange(20, 61, step=10))
    # ax1.set_yticks(np.arange(0, 101, step=25))

    # ax1.set_xticks(np.arange(0, 31, step=10))
    ax1.set_xticks(np.arange(0, 61, step=15))



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
    ax1.legend(title='KCl (mM)',bbox_to_anchor=[1, 0.5], loc='best')

    if save_fig:
        plt.savefig(save_fig, dpi=300,bbox_inches='tight')

    plt.show()





def plt_all_hog1_cis(model_fxns, top_params, params_constants, initials,  sig, time, param,
                        exp_data=None, exp_time=None, ss=False, ptpD=False, ci= 95,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)

    if save_fig:
        fig, (ax1) = plt.subplots(1, 1, figsize=(2.25,2))
    else:
        fig, (ax1) = plt.subplots(1, 1, figsize=(5,3))

    plt.rc('ytick', labelsize=10)
    plt.rc('xtick', labelsize=10)

    palette = palettes.get(param)

    llabels = ['Inactive Cytosolic', 'Active Cytosolic', 'Active Nuclear', 'Inactive Nuclear']


    # dashes = None
    # if t100a:
    dashes= (2,2)

    sims = fit_data_to_list(model_fxns, top_params, params_constants, initials, time, param, sig, ptpD=False, ss=ss)
    # for sim in sims:
    ax1 = sns.tsplot(sims[0], time,  ci = ci, dashes = dashes, color = '#CCBB44')
    ax1 = sns.tsplot(sims[1], time,  ci = ci, dashes = dashes, color = '#228833')

    ax1 = sns.tsplot(sims[2], time,  ci = ci, dashes = dashes, color = '#4477AA')
    ax1 = sns.tsplot(sims[3], time,  ci = ci, dashes = dashes, color = '#66CCEE')

    # ax1.legend(labels=labels, title='Active',loc='best')
    plt.legend(labels=llabels, title='Hog1', loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.legend(labels=labels, title='Active', loc='center left', bbox_to_anchor=(0.5, 0.5))



    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=.5)
    # ax.grid(color='grey', linestyle='-', axis='y', linewidth=.5)

    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # ax1.legend(bbox_to_anchor=[1, 0.5], loc='best')
    ax1.set_ylim(-2,102)

    if sig == 150:
        ax1.set_xlim(-1,21)
    elif sig == 350:
        ax1.set_xlim(-1.5,45)
    elif sig == 550:
        ax1.set_xlim(-1.5,60)

    if save_fig:
        plt.savefig(save_fig+".png", dpi=300,bbox_inches='tight')

    plt.show()

def plt_nuc_vs_cyto(model_fxns, top_params, params_constants, initials,  sig, time, param,
                        exp_data=None, exp_time=None, ss=False, ptpD=False, ci= 95,
                        save_fig=''):
    plt.clf()
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)

    if save_fig:
        fig, (ax1) = plt.subplots(1, 1, figsize=(2.25,2))
    else:
        fig, (ax1) = plt.subplots(1, 1, figsize=(5,3))

    plt.rc('ytick', labelsize=10)
    plt.rc('xtick', labelsize=10)

    palette = palettes.get(param)

    # llabels = ['Inactive Cytosolic', 'Active Cytosolic', 'Active Nuclear', 'Inactive Nuclear']


    # dashes = None
    # if t100a:
    dashes= (2,2)

    sims = fit_data_to_list(model_fxns, top_params, params_constants, initials, time, param, sig, ptpD=False, ss=ss)
    # for sim in sims:
    # ax1 = sns.tsplot(sims[0], time,  ci = ci, dashes = dashes, color = '#CCBB44')

    # active v active
    # ax1.plot(np.average(sims[1], axis=0), np.average(sims[2], axis=0),  color = '#228833')
    # print(len(np.average(sims[0]+sims[1], axis=0)))
    # nuc vs cytos
    closest_idxs = [np.abs(time - t).argmin() for t in exp_time]

    IC = np.average(sims[0], axis=0)[closest_idxs]
    AC = np.average(sims[1], axis=0)[closest_idxs]
    AN = np.average(sims[2], axis=0)[closest_idxs]
    IN = np.average(sims[3], axis=0)[closest_idxs]
    # print(IC+AC)
    # print(len(IC+AC))
    ax1.plot(AN+AC, AN+IN,  '-', color = '#414042' )
    if sig ==150:
        m = 'd'
    elif sig ==350:
        m = 'o'
    elif sig ==550:
        m = 's'
    else:
        m = 'o'
    ax1.plot(AN+AC, AN+IN,  m, color = palette.get(sig),mec='black' )


    ax1.set_ylim(19,61)
    # ax1.set_ylim(-5,105)

    # ax1.set_xlim(-1,31)
    ax1.set_xlim(-1,102)

    ax1.set_yticks(np.arange(10, 61, step=10))
    # ax1.set_yticks(np.arange(0, 101, step=25))

    # ax1.set_xticks(np.arange(0, 31, step=10))
    ax1.set_xticks(np.arange(0, 105, step=25))
    # ax1 = sns.tsplot(sims[2], time,  ci = ci, dashes = dashes, color = '#4477AA')
    # ax1 = sns.tsplot(sims[3], time,  ci = ci, dashes = dashes, color = '#66CCEE')

    # ax1.legend(labels=labels, title='Active',loc='best')
    # plt.legend(labels=llabels, title='Hog1', loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.legend(labels=labels, title='Active', loc='center left', bbox_to_anchor=(0.5, 0.5))



    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=.5)
    ax1.grid(color='grey', linestyle='-', axis='x', linewidth=.5)

    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # ax1.legend(bbox_to_anchor=[1, 0.5], loc='best')
    # ax1.set_ylim(-2,102)

    # if sig == 150:
    #     ax1.set_xlim(-1,21)
    # elif sig == 350:
    #     ax1.set_xlim(-1.5,45)
    # elif sig == 550:
    #     ax1.set_xlim(-1.5,60)

    if save_fig:
        plt.savefig(save_fig+".png", dpi=300,bbox_inches='tight')

    plt.show()



def simdata_to_list(model_fxns, top_params, params_constants, initials, time, param,
                        ss = False):
    sims = []
    for idx, params in enumerate(top_params):
        try:
            if ss:
                ss_data = fsolve(model_fxns.m, initials, args=(0,params_constants, 0, params))
                data = simulate_wt_experiment(model_fxns.m, ss_data, params_constants, 350000, params, time, run_type=['ramp'])
                active = data[:,2]/params_constants[2]*100
                sims.append(active)
            else:
                data = simulate_wt_experiment(model_fxns.m, initials, params_constants, 0, params, time, run_type=['man'])
        except RuntimeWarning:
            print("Runtimewarning at idx: "+str(idx))

        if idx % int(len(top_params)*.1) == 0:
            print(str(int(idx/len(top_params)*100)) + "% complete.")

    return sims

def plt_ramp_cis(sims, time, hog1_ramp_data=None, mapk_ramp_time=None, ci=95,
                        save_fig=''):

    fig, (ax1) = plt.subplots(1, 1, figsize=(5,4))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax1 = sns.tsplot(sims, time,  ci = ci)
    ax1.set_ylim(-5,105)
    ax1.set_xlim(-2,61)

    ax1.set_yticks(np.arange(0, 101, step=25))
    ax1.set_xticks(np.arange(0, 61, step=15))


    ax1.grid(color='grey', linestyle='-', axis='y', linewidth=1)


    if hog1_ramp_data:
        ax1.plot(mapk_ramp_time, hog1_ramp_data[0], 'o', markersize=10, color='Black')


    if save_fig:
        plt.savefig(save_fig, dpi=300,bbox_inches='tight')
    # plt.show()


def inhibdata_to_list(model_fxns, top_params, params_constants, initials, time, param, sig, run_type=None,
                        ss = False):
    sims = []
    for idx, params in enumerate(top_params):
        if ss:
            ss_data = run_ss(model_fxns.inhib, initials, params_constants, params)
            data = simulate_inhib_experiment(model_fxns.inhib, ss_data, params_constants, sig, params, time, run_type)
        else:
            data = simulate_inhib_experiment(model_fxns.inhib, initials, params_constants, sig, params, time, run_type)
        active = data[:,2]/params_constants[2]*100
        sims.append(active)
        if idx % int(len(top_params)*.1) == 0:
            print(str(int(idx/len(top_params)*100)) + "% complete.")
    return sims
