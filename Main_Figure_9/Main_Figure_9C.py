import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mdtraj as md
import csv
import glob
import seaborn as sns
import pandas as pd
from matplotlib import rc

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


if __name__ == '__main__':
    energy = {}
    proteins = ['CB1','CB2']
    states = ['B1','B1']
    structural_element = ['N-term','TM1','TM2','ECL1','TM3','TM4','ECL2','TM5','TM6','ECL3','TM7']

    for i,protein in enumerate(proteins):
        energy[protein] = {key: 0 for key in structural_element}
        df = pd.read_csv('./' + protein + '_resid.csv',sep=',')
        for bt in range(3):
            filenames = glob.glob('../residue_ligand_interaction/' + protein + '_assym_binding_' + states[i] + '_' + str(bt) + '_*.out')
            for filename in filenames:
                residue = filename.split('../residue_ligand_interaction/' + protein + '_assym_binding_' + states[i] + '_' + str(bt) + '_')[1].split('.out')[0]
                if '.' in df['Position'][int(residue)-1]:
                    energy[protein]['TM' + df['Position'][int(residue)-1][0]] -= np.mean(np.sum(np.loadtxt(filename)[:,1:],axis=1))/3
                else:
                    energy[protein][df['Position'][int(residue)-1]] -= np.mean(np.sum(np.loadtxt(filename)[:,1:],axis=1))/3

        print(protein,energy[protein])
            #energy[protein][se]
    df = pd.DataFrame.from_dict(energy,orient='index')
    fig, axs = plt.subplots(1,1,figsize=(20,7))
    sns.heatmap(df,vmin=0,vmax=20,cbar=True,ax=axs,cmap=mpl.colormaps['Blues'],linewidths=1,square=True,linecolor='Black',cbar_kws={"shrink": 1.0,"fraction":0.15,"aspect":30,"label":'Interaction Energy (kcal/mol)'})
    axs.tick_params(labelsize=30)
    axs.set_xticklabels(structural_element,rotation = 0,**hfont,fontsize=16)
    axs.set_yticklabels(['CB1','CB2'],rotation = 0,**hfont,fontsize=16)
    plt.tight_layout()
    plt.savefig('./Main_Figure_9C.png',transparent=False,dpi =300)

