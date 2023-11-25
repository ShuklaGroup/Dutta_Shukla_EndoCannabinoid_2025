import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import pyemma
import glob
import os

def update_index(A,active_set):
    new_index = []
    for i in A:
        if i in active_set:
            new_index.append(np.where(active_set==i)[0][0])

    return new_index


def TPT_bt(A,B,tech):
    filenames = glob.glob('../MSM_obj/CB1_assym_bt_80_*_msm.pkl') #bootstrapped MSM object with 80% of the data
    timescale = []
    for i,filename in enumerate(filenames):
        msm =  pickle.load(open(filename,'rb'))
        active_set = msm.active_set
        A_new = update_index(A,active_set)
        B_new = update_index(B,active_set)
        tpt = pyemma.msm.tpt(msm,A_new,B_new)
        timescale.append(tpt.mfpt)
        del msm
        del tpt
        del A_new
        del B_new
        del active_set
    return np.mean(timescale),np.std(timescale)


if __name__=='__main__':
    protein = 'CB1'
    tech = 'MSM'

    macro_cluster_list = pickle.load(open('./CB1_macro_to_micro_indices.pkl','rb'))

    macrostate_names = {'TM17':'I1$_{CB1}^{TM17}$','Bulk1':'U2$_{CB1}^{Mem}$','Ex':'U$_{CB1}^{Sol}$','Bulk2':'U1$_{CB1}^{Pro}$','AB':'I2$_{CB1}^{TM17}$','B1':'B$_{CB1}^{Pro}$'} 
    count = 0

    for ind1,mc in enumerate(list(macrostate_names.keys())):
        temp = list(macrostate_names.keys()).copy()
        temp.remove(mc)
        for tp in temp:
            ind2 = list(macrostate_names.keys()).index(tp)
            mean,std = TPT_bt(macro_cluster_list[ind1],macro_cluster_list[ind2],tech)
            print(macrostate_names[mc] + ' to ' + macrostate_names[tp] + ' MFPT: ' + str(mean/10000) + ' +/- ' + str(std/10000))




