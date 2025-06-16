import numpy as np 
import os
import matplotlib.pyplot as plt 
from matplotlib import rc
import pickle
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def density_plot(ax):
    for i,protein in enumerate(proteins):
        all_density = np.empty([25,3])
        for j in range(3):
            hist,bins = np.histogram(pocket_volume[protein]['bt_' +str(j)],bins=25,range=(pocket_volume[protein]['min'],pocket_volume[protein]['max']),density=True)
            all_density[:,j] = hist

        mean = np.mean(all_density,axis=1)
        std = np.std(all_density,axis=1)

        x_edge = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        print(pocket_volume[protein]['max'])
        ax.plot(x_edge, mean, linewidth=5, label=proteins[i], color=colors[i])
        ax.fill_between(x_edge, mean-std, mean+std, alpha=0.3, color=colors[i])

    ax.set_ylim([0.0,0.02])
    ax.set_yticks(np.arange(0.0,0.021,0.01),np.round(np.arange(0.0,0.021,0.01),2),fontsize=22,fontname='Helvetica')

    ax.set_xlim([0,300])
    ax.set_xticks(np.arange(0,301,100),np.arange(0,301,100),fontsize=22,fontname='Helvetica')

    ax.set_xlabel('Orthosteric Pocket Volume (\AA $^3$)',fontsize=24,fontname='Helvetica')
    ax.set_ylabel('Probability Density',fontsize=24,fontname='Helvetica')

    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig('./Main_Figure_9B.png',dpi=300)
    plt.close()


if __name__ == '__main__':
    proteins = ['CB1','CB2']
    state = 'B1_ortho'
    colors = ['orange','green']

    pocket_volume = {}
    for protein in proteins:
        pocket_volume[protein] = {}
        min_pocket_volume = 1000
        max_pocket_volume = 0
        
        for i in range(3):
            flag = []
            for j in range(1,1001):
                if os.path.exists('../pb_analysis/' + protein + '_assym_binding_' + state + '_' + str(i) + '_' + str(j) + 'volumes.tabbed.txt'):
                    temp = np.loadtxt('../pb_analysis/' + protein + '_assym_binding_' + state + '_' + str(i) + '_' + str(j) + 'volumes.tabbed.txt')
                    flag.append(temp[1])
                    
                    if temp[1] > max_pocket_volume:
                        max_pocket_volume = temp[1]
                    
                    if temp[1] < min_pocket_volume:
                        min_pocket_volume = temp[1]

            
            pocket_volume[protein]['bt_' +str(i)] = flag
            pocket_volume[protein]['min'] = min_pocket_volume
            pocket_volume[protein]['max'] = max_pocket_volume
    pocket_volume = pickle.load(open('./orthostaric_volume.pkl','rb'))
    fig,ax = plt.subplots(1,1,figsize=(10,7))
    density_plot(ax)
