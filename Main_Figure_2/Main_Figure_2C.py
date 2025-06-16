import numpy as np
import sys
import matplotlib.pyplot as plt
import mdtraj as md
import pickle
from matplotlib import rc
import matplotlib as mpl
import tol_colors as tc
from sklearn.cluster import KMeans

hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

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

def selection(extra_z,intra_z,totcom,mtM,ms):
    dtrajs_ref  = pickle.load(open('./CB1_frame_cluster_indices.pkl','rb'))

    x = []
    y = []
    new_proj = []

    for dist,dj in zip(totcom,dtrajs_ref):
        cond_1 = np.where(dist[:,2]<extra_z,1,0)
        cond_2 = np.where(dist[:,2]>intra_z,1,0)
        cond = np.multiply(cond_1,cond_2)
        x.append(dist[np.where(cond==1)[0],0])
        y.append(dist[np.where(cond==1)[0],1])
        new_proj.append(mtM[dj[np.where(cond==1)[0]]])


    x_data = np.concatenate(x)*10
    y_data = np.concatenate(y)*10
    c_data = np.concatenate(new_proj)

    c_data = np.array([new_macrostate_index[c] for c in c_data])
    scatter_plot(x_data,y_data,c_data,ms)

def scatter_plot(x_data,y_data,c_data,ms):
    colors = ['salmon','skyblue','lightcyan','lightgreen','gold','mediumpurple','orange','black','violet']
    new_colors = [colors[i] for i in list(new_macrostate_index.keys())]
    cmap = mpl.colors.ListedColormap(new_colors)
    fig, axs = plt.subplots(1,1,figsize=(10,7))
    sc = plt.scatter(x_data,y_data,c=c_data,vmin =np.min(c_data),vmax = np.max(c_data),cmap=cmap,s=4)
    cbar = fig.colorbar(sc)
    cbar.ax.set_yticks(np.linspace(0,5,7)[:ms] + np.linspace(0,5,7)[1]/2,['U$_{CB1}^{Sol}$','U1$_{CB1}^{Mem}$','U2$_{CB1}^{Mem}$','I1$_{CB1}^{TM17}$','I2$_{CB1}^{TM17}$','B$_{CB1}^{Pro}$'],fontsize=20)

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
    plt.savefig('./Main_Figure_3C.png',transparent=False,dpi =500)
    plt.close()


if __name__ == '__main__':
    protein = 'CB1'
    ms =  6
    new_macrostate_index = {2:0,3:1,1:2,0:3,4:4,5:5}

    totdist = pickle.load(open('./CB1_lig_com_cor.pkl','rb'))
    weights = pickle.load(open('./CB1_msm_weight.pkl','rb'))
    mtM = pickle.load(open('./CB1_micro_to_macro_index.pkl','rb'))
    extra_z,intra_z,extra_x,extra_y = helix_cor_cal(protein)

    x_lim = [-45,45]
    y_lim = [-45,45]

    selection(extra_z,intra_z,totdist,mtM,ms)


