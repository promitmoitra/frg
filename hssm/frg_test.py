import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("frg_norm_hssm.csv")
##subid=7873
##subid_filt = df["subj_idx"]==subid
trial_time = 3
stress_dict = {1:'pre',2:'post'}
tt_dict = {5:'short',20:'long'}
leave_resp_filt = df["response"]==0

def cum_rew_func(df,subid,stress_cond,tt_cond,ax,c):
    subid_filt = df["subj_idx"]==subid
    stress_filt = df["condition"]==stress_cond
    tt_filt = df["travel_time"]==tt_cond
    patch_leave_idx = pd.Index(df[subid_filt & stress_filt & tt_filt & leave_resp_filt].index)
    patch_entry_idx = pd.Index([df[subid_filt & stress_filt & tt_filt].index[0]]+list(patch_leave_idx+1))
    patch_trial_idx = [list(range(i,j+1)) for (i,j) in list(zip(patch_entry_idx,patch_leave_idx))]
#    rel_patch_trial_idx = [np.array(i)-i[0] for i in patch_trial_idx]
    rel_leave_idxs = list(patch_leave_idx - patch_entry_idx[0])

##    for idxs in patch_trial_idx:
##        patch_num = patch_trial_idx.index(idxs)
##        if len(df.loc[idxs]["travel_time"].unique())>1:
##            patch_trial_idx.remove(idxs)

    env_cum_rew = []
    patch_cum_rew = []
    for idxs in patch_trial_idx:
        patch_num = patch_trial_idx.index(idxs)
        tot_rew = np.array(df.loc[idxs][:]["tot_reward"])
        patch_rew = tot_rew-tot_rew[0]
#        env_cum_rew.append(tot_rew)
        env_cum_rew+=list(tot_rew)
        patch_cum_rew.append(patch_rew)
##        ax.plot(range(len(idxs)-1),np.diff(patch_rew),'--'+c)#,label="Patch {}".format(patch_num))
##        ax.legend()
##    env_cum_rew = np.concatenate([i for i in env_cum_rew])
##    ax.plot(env_cum_rew-env_cum_rew[0],'-'+c,label=str(subid)+' '+str(stress_dict[stress])+str(tt_dict[tt]))
    x_data = np.array(range(0,trial_time*(len(env_cum_rew)),trial_time))
    x_data[rel_leave_idxs:] = x_data[rel_leave_idxs:]*(tt_cond/trial_time)
    print(x_data)
    ax.plot(x_data[1:],np.diff(env_cum_rew)/trial_time,'-'+c,label=str(subid)+' '+str(stress_dict[stress])+str(tt_dict[tt]))
    return env_cum_rew,patch_cum_rew,ax


fig,ax=plt.subplots()

subid = 7873; 

stress = 1; tt = 5
preshort_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C0')

stress = 1; tt = 20
prelong_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C1')

##stress = 2; tt = 5
##postshort_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C2')
##
##stress = 2; tt = 20
##postlong_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C3')
##
##subid = 45194; 
##
##stress = 1; tt = 5
##preshort_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C4')
##
##stress = 1; tt = 20
##prelong_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C5')
##
##stress = 2; tt = 5
##postshort_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C6')
##
##stress = 2; tt = 20
##postlong_env_rew,_,ax = cum_rew_func(df,subid,stress,tt,ax,'C7')

##handles,_ = ax.get_legend_handles_labels()
##post_idx = len(handles)
##long_patch_len = max([len(i) for i in all_rew])
##num_patches = len(all_rew)
##all_cum_rew = np.nan*np.ones((num_patches,long_patch_len))
##
####stay_lock
##for i in range(num_patches):
##    all_cum_rew[i] = list(all_rew[i])+(long_patch_len-len(all_rew[i]))*[np.nan]
##
##pre_mean_rew = np.nanmean(all_cum_rew,0)

##post_env_rew,post_patch_rew,ax = cum_rew_func(df,subid,2,ax,'o')
##handles,_ = ax.get_legend_handles_labels()
###ax.legend([handles[0],handles[post_idx]],["Pre-stress","Post-stress"])
##
##long_patch_len = max([len(i) for i in all_rew])
##num_patches = len(all_rew)
##all_cum_rew = np.nan*np.ones((num_patches,long_patch_len))
##
####stay_lock
##for i in range(num_patches):
##    all_cum_rew[i] = list(all_rew[i])+(long_patch_len-len(all_rew[i]))*[np.nan]
##
##post_mean_rew = np.nanmean(all_cum_rew,0)
##ax.plot(pre_mean_rew,'oC0')
##ax.plot(post_mean_rew,'+C3')

ax.legend()
plt.show()

