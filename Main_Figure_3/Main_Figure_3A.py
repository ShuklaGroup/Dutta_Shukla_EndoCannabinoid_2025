import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import mdtraj as md
from matplotlib import rc

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


def bar_plot(stable_residues):
    colors = ['mediumpurple','gold','salmon']
    contact_prob = pickle.load(open('CB1_state_wise_contact_probability.pkl','rb'))
    contact_energy = pickle.load(open('CB1_state_wise_contact_energy.pkl','rb'))
    random_number = np.random.normal(0,0.05,3)
    residues_ticks = pickle.load(open('CB1_residues_ticks.pkl','rb'))

    fig,ax1 = plt.subplots(1,1,figsize=(14,5))

    ax2 = ax1.twinx()
    ax2.invert_yaxis()

    x1 = np.array([i-0.2 for i in range(1,len(stable_residues)*states+1)])
    x2 = np.array([i+0.2 for i in range(1,len(stable_residues)*states+1)])
    width = 0.4

    ax1.bar(x1,np.mean(contact_prob,axis=2).reshape(-1),yerr=np.std(contact_prob,axis=2).reshape(-1),width=width,color='blue',edgecolor='black',linewidth=0.3)

    ax2.bar(x2,np.mean(contact_energy,axis=2).reshape(-1),yerr=np.std(contact_energy,axis=2).reshape(-1),width=width,color='orange',edgecolor='black',linewidth=0.3)

    ax1.plot([x1[0]+0.2+len(stable_residues)-0.5 for i in range(10)],np.linspace(0,1,10),'--',c='black')
    ax1.plot([x1[0]+0.2+2*len(stable_residues)-0.5 for i in range(10)],np.linspace(0,1,10),'--',c='black')

    ax1.fill_between(np.arange(round(x1[0])-1,x1[0]+0.2+len(stable_residues)-0.45,0.1),0,1,alpha=0.1,color=colors[0])
    ax1.fill_between(np.arange(x1[0]+0.2+len(stable_residues)-0.5,x1[0]+0.2+2*len(stable_residues)-0.45,0.1),0,1,alpha=0.1,color=colors[1])
    ax1.fill_between(np.arange(x1[0]+0.2+2*len(stable_residues)-0.5,round(x1[-1])+1.05,0.1),0,1,alpha=0.1,color=colors[2])

    ax1.set_xlim([round(x1[0])-0.7,round(x1[-1])+0.7])
    ax1.set_ylim([0.0,1.0])
    ax2.set_ylim([0,-8])
    ax1.set_yticks(np.arange(0,1.1,0.2),np.round(np.arange(0,1.1,0.2),1),**hfont,fontsize=18)
    ax1.set_ylabel('Contact Probability',**hfont,color='blue',fontsize=20)

    ax2.set_yticks(np.arange(-8,1,2),np.arange(-8,1,2),**hfont,fontsize=18)
    ax2.set_ylabel('Interaction Energy (kcal/mol)', **hfont, color='orange',fontsize=20)

    ax1.set_xticks([i for i in range(1,len(stable_residues)*states+1)],residues_ticks,fontsize=16,rotation=90)

    plt.tight_layout()
    plt.savefig('./Main_Figure_4A.png',dpi=300)
    plt.close()


if __name__=='__main__':
    states = 3
    stable_residues = pickle.load(open('CB1_stable_residues.pkl','rb'))
    bar_plot(stable_residues)
