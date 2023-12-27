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

def heatmap(contact):
    state = 'TM45'
    lig_moieties = ['Polar     Head','Propyl Linker','Polyene Linker','Tail']
    headgroups = ['PC','PE','SM','PS','GL','PI','PA','DAG','CER','CHL']
    df = pd.DataFrame.from_dict(contact[state],orient='index')
    fig, axs = plt.subplots(1,1,figsize=(10,5),sharey=True,gridspec_kw=dict(left=0.13, right=0.97, bottom=0.10, top=0.99,wspace=0.05, hspace=0.01))
    cbaxes = fig.add_axes([0.13, 0.07, 0.84, 0.04])
    sns.heatmap(df,vmin=0,vmax=0.8,cbar=True,ax=axs,cmap=mpl.colormaps['Blues'],linewidths=1,square=True,linecolor='Black',cbar_ax=cbaxes,cbar_kws={'orientation':'horizontal','ticks':[0.0,0.2,0.4,0.6,0.8]})
    axs.set_yticklabels([textwrap.fill(text,11) for text in lig_moieties],rotation = 0,**hfont,fontsize=16)
    axs.set_xticklabels(headgroups,rotation = 0,**hfont,fontsize=16)
    cbaxes.set_xticklabels([0.0,0.2,0.4,0.6,0.8],fontsize=18,fontname='Helvetica')
    plt.tight_layout()
    plt.savefig('./Main_Figure_10C.png',transparent=False,dpi =300)


if __name__=='__main__':
    contact = pickle.load(open('CB2_lipid_contact.pkl','rb'))
    heatmap(contact)
