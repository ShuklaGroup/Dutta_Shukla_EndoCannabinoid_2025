import numpy as np 
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def density_plot(ax):
    state_labels = ['U1$_{CB1}^{Mem}$','U2$_{CB1}^{Mem}$']

    for i,state in enumerate(macrostate_names):
        all_density = np.empty([25,3])
        for j in range(3):
            hist,bins = np.histogram(pocket_volume[state]['bt_' +str(j)],bins=25,range=(pocket_volume[state]['min'],pocket_volume[state]['max']),density=True)
            all_density[:,j] = hist

        mean = np.mean(all_density,axis=1)
        std = np.std(all_density,axis=1)

        x_edge = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        ax.plot(x_edge, mean, linewidth=5,label=state_labels[i],color=colors[i])
        ax.fill_between(x_edge, mean-std, mean+std, alpha=0.3, color=colors[i])

    ax.set_ylim([0.0,0.06])
    ax.set_yticks(np.arange(0.0,0.07,0.02),np.round(np.arange(0.0,0.07,0.02),2),fontsize=18,fontname='Helvetica')

    ax.set_xlim([0,100])
    ax.set_xticks(np.arange(0,126,25),np.arange(0,126,25),fontsize=18,fontname='Helvetica')

    ax.set_xlabel('TM4-TM5 Pocket Volume (\AA $^3$)',fontsize=24,fontname='Helvetica')
    ax.set_ylabel('Probability Density',fontsize=24,fontname='Helvetica')

    #plt.legend(fontsize=22)
    plt.savefig('CB1_TM45_pocket_volume.png',dpi=300)
    plt.close()


if __name__ == '__main__':
    macrostate_names = ['Bulk1_position_TM45','Bulk2_position_TM45']
    colors = ['skyblue','lightgreen']

    pocket_volume = {}
    for state in macrostate_names:
        pocket_volume[state] = {}
        min_pocket_volume = 1000
        max_pocket_volume = 0
        
        for i in range(3):
            flag = []
            for j in range(1,1001):
                if os.path.exists('../pb_analysis/CB1_assym_binding_' + state + '_' + str(i) + '_' + str(j) + 'volumes.tabbed.txt'):
                    temp = np.loadtxt('../pb_analysis/CB1_assym_binding_' + state + '_' + str(i) + '_' + str(j) + 'volumes.tabbed.txt')
                    flag.append(temp[1])
                    
                    if temp[1] > max_pocket_volume:
                        max_pocket_volume = temp[1]
                    
                    if temp[1] < min_pocket_volume:
                        min_pocket_volume = temp[1]

            
            pocket_volume[state]['bt_' +str(i)] = flag
            pocket_volume[state]['min'] = min_pocket_volume
            pocket_volume[state]['max'] = max_pocket_volume
            
    fig,ax = plt.subplots(1,1,figsize=(10,7))
    density_plot(ax)
