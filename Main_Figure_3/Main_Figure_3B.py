import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import tol_colors as tc
import textwrap
import mdtraj as md
import pickle
import os

hfont = {'fontname':'Helvetica','fontweight':'bold'}

def heatmap_allpanel(contact):
    lig_moieties = ['Polar     Head','Propyl Linker','Polyene Linker','Tail']
    headgroups = ['PC','PE','SM','PS','GL','PI','CER','CHL']

    fig = plt.figure(figsize=(14,7))
    fig,axs = plt.subplots(1,3,figsize=(14,7),sharey=True,gridspec_kw=dict(left=0.07, right=0.97, bottom=0.10, top=0.99,wspace=0.05, hspace=0.01))

    cbaxes = fig.add_axes([0.05, 0.05, 0.92, 0.04])

    states = ['B1','AB','TM17'] #state B1 represents the B$_{CB1}^{Pro}$; AB represents the I2$_{CB1}^{TM17}$; TM17 represents the I1$_{CB1}^{Pro}$
    for i,state in enumerate(states):
        df = pd.DataFrame.from_dict(contact[state],orient='index')
        sns.heatmap(df,vmin=0,vmax=0.8,cbar=True,ax=axs[i],square=True,cmap=mpl.colormaps['Blues'],linewidths=1,linecolor='Black',cbar_ax=cbaxes,cbar_kws={'orientation':'horizontal','ticks':[0.0,0.2,0.4,0.6,0.8]})
        axs[i].set_xticklabels(headgroups,rotation = 0,**hfont,fontsize=16)

    axs[0].set_yticklabels([textwrap.fill(text,11) for text in lig_moieties],rotation = 0,**hfont,fontsize=16)
    cbaxes.set_xticklabels([0.0,0.2,0.4,0.6,0.8],fontsize=18,fontname='Helvetica')

    plt.tight_layout()
    plt.savefig('./Main_Figure_4B.png',transparent=False,dpi =300)



if __name__=='__main__':
    contact = pickle.load(open('CB1_lipid_contact.pkl','rb'))
    heatmap_allpanel(contact)
