import os
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def bar_plot():
    fig,axs = plt.subplots(1,1,figsize=(15,7))
    systems = list(energies.keys())
    width = 0.4
    
    for i in range(len(Molecules)):
        x1 = i - 0.2
        x2 = i + 0.2
        
        if Molecules[i] in ['AM6','AM8','AM10']:
            axs.bar(x1,energies[systems[0]][0][i],yerr=energies[systems[0]][1][i],width=0.4,color='Orange',alpha=0.7,edgecolor='black',linewidth=1,hatch='O')

        else:
            axs.bar(x1,energies[systems[0]][0][i],yerr=energies[systems[0]][1][i],width=0.4,color='Orange',alpha=0.7,edgecolor='black',linewidth=1,hatch='/')
    
        if Molecules[i] in ['AM3']:
            axs.bar(x2,energies[systems[1]][0][i],yerr=energies[systems[1]][1][i],width=0.4,color='Green',alpha=0.7,edgecolor='black',linewidth=1,hatch='O')

        else:
            axs.bar(x2,energies[systems[1]][0][i],yerr=energies[systems[1]][1][i],width=0.4,color='Green',alpha=0.7,edgecolor='black',linewidth=1,hatch='/')



    axs.set_xlim([-0.5,len(Molecules)-0.5])
    axs.set_ylim([-5,5])

    x = np.arange(0,len(Molecules),1)
    axs.set_xticks(x,labels,**hfont, fontsize=20,rotation=30)
    axs.set_yticks(np.arange(-5,6,2.5),np.round(np.arange(-5,6,2.5),1),**hfont,fontsize=18)

    axs.set_ylabel('$\Delta\Delta G$',**hfont,fontsize=28)
    #plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('./CB_ddg_comparison.png',dpi=300)
    plt.close()





if __name__=='__main__':
    Molecules = ['R-Meth','AM3','AM6','ACPA','AM7','AM8','AM10','AM883','AM12','AM17'] 
    labels = ['120526','346659','162244','418353','150942','421967','359278','274016','74909','151504']

    energies = {'CB1_system_0':[np.zeros(len(Molecules)), np.zeros(len(Molecules))],'CB2_system_2': [np.zeros(len(Molecules)), np.zeros(len(Molecules))]}

    for system in list(energies.keys()):
        protein = system[:3]
        for i,molecule in enumerate(Molecules):
            lig = open('../relative_energy/' + molecule + '.txt','r').readlines()
            rec_com = open('../relative_energy/' + protein + '_' + molecule + '.txt','r').readlines()
            DD_G_lig = float(lig[-1].split()[1])
            DD_G_lig_std = float(lig[-1].split()[3])

            DD_G_com = float(rec_com[-1].split()[1])
            DD_G_com_std = float(rec_com[-1].split()[3])
    
            energies[system][0][i] = DD_G_com - DD_G_lig
            
            energies[system][1][i] = DD_G_com_std + DD_G_lig_std
            print(protein,molecule,np.round(energies[system][0][i],1),np.round(energies[system][1][i],1))
    bar_plot()
    
        
                
