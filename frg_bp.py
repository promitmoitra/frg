import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

subid = 7873

channel = 'CZ'

stress_flag = 'pre'
##stress_flag = 'post'

tt_flag = 'short'
##tt_flag = 'long'

stay_lock_res = sio.loadmat(str(subid)+'_'+channel+'_'+stress_flag+tt_flag+'_stay_lock_res.mat')
leave_lock_res = sio.loadmat(str(subid)+'_'+channel+'_'+stress_flag+tt_flag+'_leave_lock_res.mat')

patch_trials = sio.loadmat(str(subid)+'_patch_trials.mat');patch_trials = patch_trials['patch_trials'][0]
num_patches = len(patch_trials)

trial_idx = dict(zip(range(num_patches),patch_trials-1))    ##-1 because python is 0-indexed
for key in trial_idx:
    trial_idx[key]=trial_idx[key][0]


