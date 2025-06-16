import numpy as np
import sys
import matplotlib.pyplot as plt
import mdtraj as md
import pickle
from matplotlib import rc
import matplotlib as mpl
import tol_colors as tc
from sklearn.cluster import KMeans

def helix_cor_cal(protein):
    protein_info = pickle.load(open('./important_proteins_with_resid.pkl','rb'))
    if protein == 'CB1':
        cnr = protein_info['CNR1']['Human']['residue info']
        strip = 88
    else:
        cnr = protein_info['CNR2']['Human']['residue info']
        strip = 0
    helices = ['1x','2x','3x','4x','5x','6x','7x']
    parm = './' +  protein + '_assym-strip.prmtop'
    t = md.load('./' + protein + '_assym-strip.pdb',top=parm)
    temp = []
    start_residue = []
    end_residue = []
    for resid in cnr:
        if resid[2].startswith(helices[0]):
            temp.append(resid[0])
        else:
            if len(temp) > 0:
                start_residue.append(temp[0])
                end_residue.append(temp[-1])
                temp = []
                helices.pop(0)

        if len(helices) == 0:
            break

    extra_cellular = np.empty([7,3])
    intra_cellular = np.empty([7,3])

    for i in range(0,len(start_residue)):
        if i%2==0:
            alpha = t.topology.select('name CA and resid ' + str(start_residue[i]-strip-1))[0]
            extra_cellular[i,:] = t.xyz[0,alpha,:]
            alpha = t.topology.select('name CA and resid ' + str(end_residue[i]-strip-1))[0]
            intra_cellular[i,:] = t.xyz[0,alpha,:]
        else:
            alpha = t.topology.select('name CA and resid ' + str(start_residue[i]-strip-1))[0]
            intra_cellular[i,:] = t.xyz[0,alpha,:]
            alpha = t.topology.select('name CA and resid ' + str(end_residue[i]-strip-1))[0]
            extra_cellular[i,:] = t.xyz[0,alpha,:]

    extra_z = np.mean(extra_cellular[:,2])
    intra_z = np.mean(intra_cellular[:,2]) + 1.5

    return extra_z,intra_z,extra_cellular[:,0],extra_cellular[:,1]


def free_energy_cal(extra_z,intra_z,extra_x,extra_y,protein,weights):
    x = []
    y = []
    new_weights = []
    for dist,wei in zip(totdist,weights):
        cond_1 = np.where(dist[:,2]<extra_z,1,0)
        cond_2 = np.where(dist[:,2]>intra_z,1,0)
        cond = np.multiply(cond_1,cond_2)
        x.append(dist[np.where(cond==1)[0],0])
        y.append(dist[np.where(cond==1)[0],1])
        new_weights.append(wei[np.where(cond==1)[0]])

    x_data = np.concatenate(x)*10
    y_data = np.concatenate(y)*10
    new_weights = np.concatenate(new_weights)

    hfont = {'fontname':'Helvetica','fontweight':'bold'}
    rc('text', usetex=True)
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

    x_bins = 200
    y_bins = 200
    R = 0.001987
    T = 300

    fig_wid = 10
    fig_hig = 7
    cmap = tc.tol_cmap('rainbow_PuBr')

    Max_energy = 8
    x_data_min =  np.min(x_data)
    y_data_min =  np.min(y_data)
    x_data_max =  np.max(x_data)
    y_data_max =  np.max(y_data)


    x_hist_lim_low =  x_data_min -0.5
    y_hist_lim_low =  y_data_min -0.5
    x_hist_lim_high = x_data_max +0.5
    y_hist_lim_high = y_data_max  +0.5

    hist= np.histogram2d(x_data,y_data, bins=[x_bins,y_bins],
                         range = [[x_hist_lim_low,x_hist_lim_high],[y_hist_lim_low,y_hist_lim_high]],density= True,weights=new_weights)

    prob_density = hist[0]
    xedge = hist[1]
    yedge = hist[2]

    x_bin_size = xedge[1]-xedge[0]
    y_bin_size = yedge[1]-yedge[0]

    free_energy = -R*T*np.log(prob_density*x_bin_size*y_bin_size)
    min_free_energy= np.min(free_energy)

    delta_free_energy = free_energy - min_free_energy

    xx = [(xedge[i]+xedge[i+1])/2 for i in range(len(xedge)-1)]
    yy = [(yedge[i]+yedge[i+1])/2 for i in range(len(yedge)-1)]

    X, Y = np.meshgrid(xx,yy)
    fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))

    cd =axs.contourf(X,Y,delta_free_energy.T, np.linspace(0,Max_energy,Max_energy*5+1),
                     vmin=0.0, vmax=Max_energy,cmap=cmap)
    cbar = fig.colorbar(cd,ticks=range(Max_energy+1))
    cbar.ax.set_yticklabels(range(Max_energy+1),fontsize=22)                                        #ticklabels of color bar
    cbar.ax.set_ylabel('Free Energy (Kcal/mol)', labelpad=15, **hfont,fontsize=24)

    plt.contour(X,Y,delta_free_energy.T, linewidths=0.25, levels=[0,1,2,3,4,5,6,7,8], colors='black')
    axs.scatter(extra_x*10,extra_y*10,s=100,marker='*',c='k')

    axs.set_xlim([x_lim[0],x_lim[1]])
    axs.set_ylim([y_lim[0],y_lim[1]])

    axs.set_xticks(range(int(x_lim[0]),int(x_lim[1])+1,15))
    axs.set_xticklabels(range(int(x_lim[0]),int(x_lim[1])+1,15))

    axs.set_yticks(range(int(y_lim[0]),int(y_lim[1])+1,15))
    axs.set_yticklabels(range(int(y_lim[0]),int(y_lim[1])+1,15))

    plt.xlabel('X (\AA)', **hfont,fontsize=24)
    plt.ylabel('Y (\AA)', **hfont,fontsize=24)

    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)

    plt.tight_layout()
    plt.savefig('./Main_Figure_3A.png',transparent=False,dpi =500)

if __name__ == '__main__':
    protein = 'CB1'
    totdist = pickle.load(open('./CB1_lig_com_cor.pkl','rb'))
    weights = pickle.load(open('./CB1_msm_weight.pkl','rb'))

    extra_z,intra_z,extra_x,extra_y = helix_cor_cal(protein)

    x_lim = [-45,45]
    y_lim = [-45,45]

    free_energy_cal(extra_z,intra_z,extra_x,extra_y,protein,weights)

