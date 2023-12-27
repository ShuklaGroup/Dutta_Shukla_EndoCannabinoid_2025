import mdtraj as md 
import numpy as np
import os
import pickle
import scipy
import sys
import matplotlib.pyplot as plt 
#import seaborn as sns
import tol_colors as tc
import pyemma
from matplotlib import rc

hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def angle_calculation(t,state,ii,nFrames):
    indices,angles = md.compute_chi1(t)
    distI = np.empty([nFrames,2])
    PHE174 = t.topology.select('resid ' + str(174-89) + ' and name N')[0]

    Nterm = t.topology.select('name CA and resid 14')
    ASPTM2 = t.topology.select('name CA and resid 74')

    temp = angles[:,np.where(indices[:,0] == PHE174)[0]]*180/np.pi
    
    distI[:,[0]] = np.where(temp<0,temp+360,temp)
    distI[:,[1]] = md.compute_distances(t,[[Nterm[0],ASPTM2[0]]],periodic=False)*10

    return distI


def density_plot(distX,distY,ax_histx,ax_histy):
    state_labels = ['M1','A','AB','B']
    #colors = ['blue','red','yellow','purple']

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    for i,state in enumerate(query_macrostates):
        all_density = np.empty([50,3])
        for j in range(3):
            hist,bins = np.histogram(distX[:,j,i],bins=50,range=(distX.min(),distX.max()),density=True)
            all_density[:,j] = hist

        mean = np.mean(all_density,axis=1)
        std = np.std(all_density,axis=1)

        x_edge = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        ax_histx.plot(x_edge, mean, linewidth=5, label=macrostates[state], color=colors[i])
        ax_histx.fill_between(x_edge, mean-std, mean+std, alpha=0.3, color=colors[i])
        
    ax_histx.set_ylim([0.0,0.1])
    ax_histx.set_yticks(np.arange(0.0,0.11,0.05),np.round(np.arange(0.0,0.11,0.05),2),fontsize=18,fontname='Helvetica')


    for i,state in enumerate(query_macrostates):
        all_density = np.empty([50,3])
        for j in range(3):
            hist,bins = np.histogram(distY[:,j,i],bins=50,range=(distY.min(),distY.max()),density=True)
            all_density[:,j] = hist

        mean = np.mean(all_density,axis=1)
        std = np.std(all_density,axis=1)

        x_edge = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        ax_histy.plot(mean, x_edge, linewidth=5,label=macrostates[state],color=colors[i])
        ax_histy.fill_betweenx(x_edge, mean-std, mean+std, alpha=0.3, color=colors[i])

    ax_histy.set_xlim([0.0,1.0])
    ax_histy.set_xticks(np.arange(0.0,1.1,0.5),np.round(np.arange(0.0,1.1,0.5),1),fontsize=18,fontname='Helvetica')


def scatter_density_plot(dist):
    fig = plt.figure(figsize=(12,12))
    Max_energy = 8
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.17, right=0.95, bottom=0.07, top=0.95,
                      wspace=0.1, hspace=0.07)

    state_labels = ['M1','A','AB','B']
    #colors = ['blue','red','yellow','purple']

    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    
    for i,state in enumerate(query_macrostates):
        query_cluster = macro_cluster_list[list(macrostates.keys()).index(state)]
        query_indices = []
        for qc in query_cluster:
            query_indices.append(np.where(dtrajs==qc)[0])

        query_indices = np.concatenate(query_indices)
        w = weights[query_indices]
        X  = txx[query_indices,0]
        Y  = txx[query_indices,1]
        X,Y,G = histogram(w,X,Y,bins=25)
        
        ax.contour(X,Y,G, linewidths=5, levels=[0,2], colors=colors[i],linestyles='dashed')

    X,Y,G = histogram(weights,txx[:,0],txx[:,1],bins=200)

    cd = ax.contourf(X,Y,G, np.linspace(0,Max_energy,Max_energy*5+1),
                     vmin=0.0, vmax=Max_energy,cmap=tc.tol_cmap('rainbow_PuBr'))
    
    cbaxes = fig.add_axes([0.07, 0.07, 0.015, 0.68])
    cbar = fig.colorbar(cd,ticks=range(Max_energy+1),cax = cbaxes)
    cbar.ax.yaxis.set_label_position("left")
    cbar.ax.yaxis.tick_left()
    cbar.ax.set_yticklabels(range(Max_energy+1),fontsize=18)                                        #ticklabels of color bar
    cbar.ax.set_ylabel('Free Energy (Kcal/mol)', labelpad=10,fontname='Helvetica',fontsize=24)

    density_plot(dist[:,0,:,:],dist[:,1,:,:],ax_histx,ax_histy)
    ax_histx.legend(fontsize=18,loc='upper left')
    ax.set_xlim([240,330])
    ax.set_ylim([15,30])

    ax.set_xticks(np.arange(240,331,30),np.arange(240,331,30),fontsize=18,fontname='Helvetica')
    ax.set_yticks(np.arange(15,31,5),np.arange(15,31,5),fontsize=18,fontname='Helvetica')

    ax.set_xlabel('F174$^{2.61}$ angle($\chi_1$)',fontsize=24,fontname='Helvetica')
    ax.set_ylabel('N-term(M103$^{N-term}$) -- TM2(D163$^{2.50}$) (\AA)',fontsize=24,fontname='Helvetica')
    #plt.legend(fontsize=18,loc='lower right')
    plt.tight_layout()
    plt.savefig('./Main_Figure_5A.png',dpi=300)
    plt.close()



def feature_load(tech):
    totdist = pickle.load(open('CB1_binding_conf_feature.pkl','rb'))
    txx = np.concatenate(totdist)
    dtrajs = np.concatenate(pickle.load(open('CB1_assym_cluster_indices.pkl','rb')))

    if tech == 'MSM':
        weights = np.concatenate(pickle.load(open('CB1_msm_weight.pkl','rb')))
    else:
        weights = np.ones(len(txx))/len(txx)

    return txx, dtrajs, weights

def main(state):
    top = './CB1_assym-strip.prmtop'
    ref = './CB1_assym-strip.pdb'
    trefI = md.load(ref,top=top)
    alpha_ca = trefI.topology.select('name CA')

    totdist = np.empty([1000,2,3])
    
    for ii in range(3):
        traj = '../macrostate_rep_dcd/CB1_assym_binding_' + state + '_' + str(ii) + '.dcd'
        t = md.load(traj,top=top)
        nFrames=t.n_frames
        t.superpose(trefI,0,alpha_ca,alpha_ca)
        totdist[:,:,ii] = angle_calculation(t,state,ii,nFrames)

    return totdist

    
def histogram(w,X,Y,bins):
    hist= np.histogram2d(X,Y, bins=[bins,bins],
                     range = [[240,330],[15,30]],
                     density= True,weights=w)

    prob_density = hist[0]
    xedge = hist[1]
    yedge = hist[2]

    x_bin_size = xedge[1]-xedge[0]
    y_bin_size = yedge[1]-yedge[0]

    R = 0.001987
    T = 300

    free_energy = -R*T*np.log(prob_density*x_bin_size*y_bin_size)
    min_free_energy= np.min(free_energy)
    delta_free_energy = free_energy - min_free_energy

    xx = [(xedge[i]+xedge[i+1])/2 for i in range(len(xedge)-1)]
    yy = [(yedge[i]+yedge[i+1])/2 for i in range(len(yedge)-1)]

    X, Y = np.meshgrid(xx,yy)
    
    return X,Y,delta_free_energy.T


if __name__=='__main__':
    number_samples = 3
    number_features = 2
    
    macrostates = {'TM17':'I1$_{CB1}^{TM17}$','Bulk1':'U2$_{CB1}^{Mem}$','Ex':'U$_{CB1}^{Sol}$','Bulk2':'U1$_{CB1}^{Pro}$','AB':'I2$_{CB1}^{TM17}$','B1':'B$_{CB1}^{Pro}$'}
    query_macrostates = ['Bulk1','TM17','AB','B1']
    colors = ['skyblue','salmon','gold','mediumpurple']
    
    macro_cluster_list = pickle.load(open('./CB1_macro_to_micro_indices.pkl','rb')) 
    
    txx, dtrajs, weights = feature_load('MSM')

    dist = np.empty([1000,number_features,number_samples,len(query_macrostates)])

    for i,state in enumerate(query_macrostates):
        dist[:,:,:,i] = main(state)

    scatter_density_plot(dist)
