import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pyemma
import mdtraj as md
import csv
import glob
from matplotlib import rc


def barplot2(H,ST):
    fig, ax1 = plt.subplots(1, 1,figsize=(22,7))
    enthalpy = np.mean(H['CB1'])-np.mean(H['CB2'])
    entropy = np.mean(ST['CB1'])-np.mean(ST['CB2'])
    ax1.barh([0.75,1.25],[-1*enthalpy,entropy],height=0.4,color=['darkblue','magenta'])

    ax1.set_xlim(-2,7)
    ax1.set_ylim(0.5,1.5)

    ax1.set_yticks([0.75,1.25],['$\Delta$H$_{CB1}^{L-P}$-$\Delta$H$_{CB2}^{L-P}$','(-T$\Delta$S$_{CB1}$)-(-T$\Delta$S$_{CB2}$)'],fontsize=24,fontname='Helvetica')

    ax1.set_xticks([-2,0,2,4,6,8],[2,0,-2,-4,-6,-8],fontsize=22,fontname='Helvetica')
    ax1.set_xlabel('Free Energy Contribution (kcal/mol)',fontsize=26,fontname='Helvetica')

    plt.tight_layout()
    plt.savefig('./Main_Figure_9D.png',dpi=300)
    plt.close()



if __name__ == '__main__':
    H = {'CB1': [],'CB2': []}
    ST = {'CB1': [],'CB2': []}

    for protein in ['CB1','CB2']:
        cluster =  1000
        dim = 6
        tic_lag_time = 250
        macrostate_names = {'CB2' :['TM45','TM56','TM17','B1','Bulk','B2'], 'CB1':['TM17','Bulk1','Ex','Bulk2','AB','B1']}
        state = 'B1'

        clusters = pickle.load(open('./' + protein + '_macro_to_micro_indices.pkl','rb'))[macrostate_names[protein].index(state)]
        for bt in range(3):
            temp = []
            filenames = glob.glob('../residue_ligand_interaction/' + protein + '_assym_binding_' + state + '_' + str(bt) + '_*.out')
            for filename in filenames:
                temp.append(np.mean(np.sum(np.loadtxt(filename)[:,1:],axis=1)))
            
            H[protein].append(np.sum(temp))
        
            msm_obj = pickle.load(open('./' + protein + '_assym_MSM.pkl','rb'))
            temp = msm_obj.stationary_distribution
            active_set = msm_obj.active_set
            stat_dis = np.zeros(cluster)
            stat_dis[active_set] = temp
            R = 0.001987
            T = 300
            flag = stat_dis[clusters]/np.sum(stat_dis[clusters])
            S = 0
            for c,_ in enumerate(clusters):
                if flag[c] == 0:
                    continue
                else:
                    S -= flag[c]*np.log(flag[c])
            ST[protein].append(R*T*S)
        
    barplot2(H,ST)
