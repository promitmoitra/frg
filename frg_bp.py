import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas as pd
import pathlib

ids = [i.stem for i in pathlib.Path().glob('./data/[0-9]*')]
for i in [ids[-1]]:
    subid = int(i)
    df = pd.read_csv('./hssm/sub_data/'+i+'_hddm_norm.csv')
    df.drop(["max_cz","max_c3","max_c4"],axis=1,inplace=True)
    channel = 'CZ'

##    stress_flag = 'pre';stress=1
##    stress_flag = 'post';stress=2
##    tt_flag = 'short';tt=5
##    tt_flag = 'long';tt=20
##    subid_filt = df["subj_idx"]==subid

    stress_filt = df["condition"]==stress
    tt_filt = df["travel_time"]==tt
    leave_filt = df["response"]==0

    patch_leave_idx = pd.Index(df[leave_filt].index)
    patch_entry_idx = pd.Index([df.index[0]]+list(patch_leave_idx+1))
    patch_trial_idx = [list(range(i,j+1)) for (i,j) in list(zip(patch_entry_idx,patch_leave_idx))]
##    for idxs in patch_trial_idx:
##        if len(df.loc[idxs]["travel_time"].unique())>1:
##            df.loc[idxs,"travel_time"] = df.loc[idxs[-1],"travel_time"]
##            df.at[idxs,"travel_time"] = df.loc[idxs[-1]]["travel_time"]
    
    stress_flag = 'pre';stress=1;tt_flag = 'short';tt=5
    stay_lock_res = sio.loadmat('./data/'+str(subid)+'/'+str(subid)+'_'+channel+'_'+stress_flag+tt_flag+'_stay_lock_res.mat')
    leave_lock_res = sio.loadmat('./data/'+str(subid)+'/'+str(subid)+'_'+channel+'_'+stress_flag+tt_flag+'_leave_lock_res.mat')
    patch_trials = sio.loadmat('./data/'+str(subid)+'/'+str(subid)+'_'+stress_flag+tt_flag+'_patch_trials.mat');patch_trials = patch_trials['patch_trials'][0]
    num_patches = len(patch_trials)
    
    trial_idx = dict(zip(range(num_patches),patch_trials-1))    ##-1 because python is 0-indexed
    for key in trial_idx:
        trial_idx[key]=trial_idx[key][0]
    

