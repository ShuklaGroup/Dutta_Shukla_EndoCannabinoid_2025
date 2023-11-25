import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import os
from sklearn.preprocessing import MinMaxScaler

def kl_divergence(p, q):
    kl_div = 0
    for i in range(len(p)):
        if not(p[i] == 0 or q[i] == 0):
            kl_div += p[i] * np.log2(p[i]/q[i])
    return kl_div


def kl_calculation(ref_state,mc):
    totdist = []
    query_filename = './CB1_assym_binding_' + mc + '_0_chv.npy'
    txx_query = 1/np.load(query_filename)

    ref_filename = './CB1_assym_binding_' + ref_state + '_0_chv.npy'
    txx_ref = 1/np.load(ref_filename)

    number_of_feature = txx_query.shape[1]

    kl_avg = np.empty([number_of_feature,1])

    for i in range(number_of_feature):
        x_data = txx_ref[:,i]
        y_data = txx_query[:,i]
        minimum = np.min(np.minimum(x_data,y_data))
        maximum = np.max(np.maximum(x_data,y_data))
        bins = np.arange(minimum,maximum,(maximum-minimum)/100)
        xhist,xedges = np.histogram(x_data,bins=bins,density=True)
        yhist,yedges = np.histogram(y_data,bins=bins,density=True)
        x_prob = xhist * np.diff(xedges)
        y_prob = yhist * np.diff(yedges)
        kl_avg[i,0] = (kl_divergence(x_prob,y_prob) + kl_divergence(y_prob,x_prob))/2

    np.save('./CB1_assym_binding_' + ref_state + '_' +  mc + '_kl_avg_score.npy',kl_avg)

    return kl_avg

def pdb_b_factor(r,q,j):
    f = open('./CB1_assym_binding_' + q +'_0.pdb','r')
    fh = open('./CB1_assym_binding_'+ r + '_' + q +'_b.pdb','w')

    data = f.readlines()
    b_values = {}
    thickness_flag = []
    thickness = []
    for i in range(len(weights_per_residue)):
        temp = "{:.2f}".format(weights_per_residue[i,int(j)])
        b_values[str(indexs[i]+1)] = str(temp)
        thickness_flag.append(float(temp)+0.1)

    for i in range(len(thickness_flag)):
        if i < 4:
            thickness.append(thickness_flag[i])
        else:
            thickness.append((thickness_flag[i] + thickness_flag[i-1] + thickness_flag[i-2]+thickness_flag[i-3]+thickness_flag[i-4])/5)
        #print(str(indexs[i]+1),thickness[i],b_values[str(indexs[i]+1)])

    for line in data:
        if line[:4] == 'ATOM' and not (int(line[23:26]) > 323):
            modified_line = line.replace('1.00  0.00','1.00  ' + b_values[str(int(line[23:26]))])
            fh.write(modified_line)
        else:
            fh.write(line)

    f.close()
    fh.close()

    vmd_writer(r,q,thickness,b_values)

def vmd_writer(r,q,thickness,b_values):
    residues = list(b_values.keys())
    f = open('./CB1_assym_binding_' + r + '_' + q + '.vmd','w')
    f.write('mol new ./CB1_assym_binding_' + r + '_' + q + '_b.pdb type pdb first 0 last -1 \n')
    for i,resid in enumerate(residues):
        if i < 4:
            f.write('mol modselect ' + str(i) + ' 0 resid ' + resid + ' \n')
        else:
            f.write('mol modselect ' + str(i) + ' 0 resid ' + resid + ' ' + residues[i-1] + ' ' + residues[i-2] + ' ' + residues[i-3] + ' ' + residues[i-4] + ' \n')
        f.write('mol modmaterial ' + str(i) + ' 0 AOChalky \n')
        f.write('mol modstyle ' + str(i) + ' 0 NewCartoon ' + str(thickness[i]) + ' 10.000000 2.100000 0 \n')
        f.write('mol modcolor ' + str(i) + ' 0 Beta \n')
        if i < 322:
            f.write('mol addrep 0 \n')
    f.write('color scale method BWR \n')
    f.write('display projection Orthographic \n')
    f.write('color Display Background white \n')
    f.write('display depthcue off \n')
    f.write('axes location Off \n')
    f.close()

def closest_atom_dist(state):
    top = './CB1_assym-strip.prmtop'
    traj = '../macrostate_rep_dcd/CB1_assym_binding_' + state + '_0.dcd'
    t = md.load(traj,top=top)
    dist,index = md.compute_contacts(t, contacts='all', scheme='closest-heavy', ignore_nonprotein=True, periodic=False)
    filename = traj.split('.dcd')[0].split('../macrostate_rep_dcd/')[1]
    np.save('./' + filename +'_chv.npy',dist)
    np.save('./' + filename+'_chv_ind.npy',index)

if __name__=='__main__':
    ref_state = 'TM17'
    macrostate_names = ['TM17','AB','B1'] #state B1 represents the B$_{CB1}^{Pro}$; AB represents the I2$_{CB1}^{TM17}$; TM17 represents the I1$_{CB1}^{Pro}$

    for j,mc in enumerate(macrostate_names):
        if not os.path.exists('./important_files/CB1_assym_binding_' + mc + '_0_chv.npy'):
            closest_atom_dist(mc)

    a = np.load('./CB1_assym_binding_TM17_0_chv_ind.npy')
    indexs = np.unique(a)
    weights_per_residue = np.empty([len(indexs),len(macrostate_names)])

    for j,mc in enumerate(macrostate_names):
        kl_file = './CB1_assym_binding_' + ref_state + '_' +  mc + '_kl_avg_score.npy'
        if not os.path.exists(kl_file):
            kl_div = kl_calculation(ref_state,mc)

        kl_div = np.load(kl_file)

        for i in range(len(indexs)):
            position = np.where(a==indexs[i])[0]
            weights_per_residue[i,j] = np.sum(kl_div[position])


    min_max_scaler = MinMaxScaler()
    min_max_scaler.fit(weights_per_residue)
    weights_per_residue = min_max_scaler.transform(weights_per_residue)

    for j,state2 in enumerate(ref_state):
        pdb_b_factor(ref_state,state2,j)


