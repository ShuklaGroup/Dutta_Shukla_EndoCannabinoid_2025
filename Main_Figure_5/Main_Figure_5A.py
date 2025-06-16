import mdtraj as md
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc

hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def distance_calculation(t,nFrames):
    distI = np.empty([nFrames,1])
    dihe = [['C5','C10','C18']]
    for i,d in enumerate(dihe):
        LIG1 = t.topology.select('resname LIG and name ' + d[0])[0]
        LIG2 = t.topology.select('resname LIG and name ' + d[1])[0]
        LIG3 = t.topology.select('resname LIG and name ' + d[2])[0]
        #LIG4 = t.topology.select('resid 324 and name ' + d[3])[0]
        temp = md.compute_angles(t,[[LIG1,LIG2,LIG3]],periodic=False)*180/np.pi
        print(temp.mean(),temp.std())
        distI[:,[i]] = np.where(temp<0,temp+360,temp)

    return distI

def main(protein):
    top = './' + protein + '_assym-strip.prmtop'
    ref = './' + protein + '_assym-strip.pdb'
    trefI = md.load(ref,top=top)
    alpha_ca = trefI.topology.select('name CA')

    totdist = np.empty([1000,1,3])

    for ii in range(3):
        traj = '../macrostate_rep_dcd/' + protein + '_assym_binding_B1_' + str(ii) + '.dcd'
        t = md.load(traj,top=top)
        nFrames=t.n_frames
        t.superpose(trefI,0,alpha_ca,alpha_ca)
        totdist[:,:,ii] = distance_calculation(t,nFrames)

    return totdist

def density_plot(dist):
    fig,ax = plt.subplots(1,1,figsize=(10,7))
    for i,protein in enumerate(proteins):
        all_density = np.empty([25,3])
        for j in range(3):
            hist,bins = np.histogram(dist[:,j,i],bins=25,range=(dist.min(),dist.max()),density=True)
            all_density[:,j] = hist

        mean = np.mean(all_density,axis=1)
        std = np.std(all_density,axis=1)

        x_edge = [(bins[ii]+bins[ii+1])/2 for ii in range(len(bins)-1)]
        ax.plot(x_edge, mean, linewidth=5, label=proteins[i], color=colors[i])
        ax.fill_between(x_edge, mean-std, mean+std, alpha=0.3, color=colors[i])


    ax.set_ylim([0.0,0.02])
    ax.set_yticks(np.arange(0.0,0.021,0.01),np.round(np.arange(0.0,0.021,0.01),2),fontsize=22,fontname='Helvetica')

    ax.set_xlim([0,180])
    ax.set_xticks(np.arange(0,181,90),np.arange(0,181,90),fontsize=22,fontname='Helvetica')

    ax.set_xlabel('Angle (degree)',fontsize=24,fontname='Helvetica')
    ax.set_ylabel('Probability Density',fontsize=24,fontname='Helvetica')

    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig('./Main_Figure_9A.png',dpi=300)
    plt.close()

if __name__=='__main__':
    number_samples = 3
    number_features = 1

    proteins = ['CB1','CB2']
    colors = ['orange','green']

    dist = np.empty([1000,number_features,number_samples,len(proteins)])

    for i,protein in enumerate(proteins):
        dist[:,:,:,i] = main(protein)

    density_plot(dist[:,0,:,:])

