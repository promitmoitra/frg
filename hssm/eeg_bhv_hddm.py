import numpy as np
import scipy.io as sio
import csv
import pathlib
import pandas as pd
from sklearn.preprocessing import StandardScaler as stds
from sklearn.preprocessing import MinMaxScaler as mms
'''
'''
def read_bhv_mat(fname,pre=True):
    counts = {}
    field_names = []
    with open("fields_bhv.txt",'r') as f:
        fields = csv.reader(f)
        for field in fields:
            field_names.append(field[0])
    field_names = field_names
    matfile = sio.loadmat(fname)
    if pre:
        prepost = matfile.get('prestress')
    else:
        prepost = matfile.get('poststress')
    field_vals = [prepost[0,0][i] for i in range(len(prepost[0,0]))]
    data = dict(zip(field_names,field_vals))
    return data

def concat_csv_to_df(filepaths):
    df = pd.DataFrame()
    for f in filepaths:
        tmp = pd.read_csv(f)
        df = pd.concat([df,tmp])
    return df

def extract_from_bhv(subid,bhvdict,preflag):
    bhv_response = []
    bhv_reward = []
    bhv_tot_reward = []
    bhv_rt = []
    if preflag and subid in [42656,48276]:
        rew_idx = 1
    else:
        rew_idx = 0
    for idx in range(len(bhvdict["CodeNumbers"][0])):
        codes = bhvdict["CodeNumbers"][0,idx].flatten()
        times = bhvdict["CodeTimes"][0,idx].flatten()
        if 5 in codes:
            bhv_response.append(1)
            bhv_reward.append(bhvdict["UserVars"][0,idx][rew_idx][0][0])
            bhv_tot_reward.append(bhvdict["UserVars"][0,idx][rew_idx+1][0][0])
            rt = times[codes==5]-times[codes==4]
            bhv_rt.append(rt[0])
        elif 8 in codes:
            bhv_response.append(0)
            bhv_reward.append(0)
            try:
                bhv_tot_reward.append(bhv_tot_reward[-1])
            except IndexError:
                bhv_tot_reward.append(0)
            rt = times[codes==8]-times[codes==4]
            bhv_rt.append(rt[0])
    return bhv_response, bhv_reward, bhv_tot_reward, bhv_rt    

if __name__=="__main__":
    bhv_data_path = '../frg/foraging/Neuroflow/'
    l0_data_path = './l0_eeg_hddm/raw/'
    all_ids = [int(i.stem) for i in pathlib.Path().glob(bhv_data_path+'[0-9]*')]
    all_ids.sort()
    bad_ids = [37532, 38058, 39862, 42125, 43543, 46037, 47744, 47801, 48238, 48278]    
    sub_ids = np.setdiff1d(all_ids,bad_ids)

    for sub in sub_ids:
        print('='*50,sub,'='*50)
        fpath = list(pathlib.Path().glob(bhv_data_path+str(sub)+'/*?oraging*/'+str(sub)+'.mat'))[0]
        pre_bhv = read_bhv_mat(fpath,pre=True)
        post_bhv = read_bhv_mat(fpath,pre=False)

        fpath = list(pathlib.Path().glob(l0_data_path+str(sub)+'*_hddm.csv'))[0]
        df_eeg = pd.read_csv(fpath)

##  ############################################################################################

        """calculate response, reward and rt from bhv trials"""
        pre_bhv_response,pre_bhv_rew,pre_bhv_tot_rew,pre_bhv_rt = extract_from_bhv(sub,pre_bhv,1)
        post_bhv_response,post_bhv_rew,post_bhv_tot_rew,post_bhv_rt = extract_from_bhv(sub,post_bhv,0)

        sub_bhv_response = pre_bhv_response + post_bhv_response
        sub_bhv_rew = pre_bhv_rew+post_bhv_rew
        sub_bhv_tot_rew = pre_bhv_tot_rew+post_bhv_tot_rew
        sub_bhv_rt = pre_bhv_rt + post_bhv_rt
        sub_bhv_rt = (np.array(sub_bhv_rt)-1000)/1000

        if len(df_eeg)==len(sub_bhv_response):
            if np.equal(df_eeg["response"].to_numpy(),np.array(sub_bhv_response)).all():
                print('#'*10,'Check','#'*10)
            else:
                print('misaligned')
        else:
            print('unequal len')
        
        resp = df_eeg["response"].to_numpy()
        leave_idx = np.where(np.array(resp)==0)[0]
    
        start = 0
        for idx in leave_idx:
            end = idx
            num_patch_trials = len(resp[start:end])
            resp[start:end] = np.array(range(num_patch_trials,0,-1))
            start = end+1

        df_eeg["response"] = sub_bhv_response[:len(df_eeg)]

##      Normalizing:
##        mmscale = mms(); stdscale = stds()
##        resp = mmscale.fit_transform(resp.reshape(-1,1)).flatten()
##        sub_bhv_rew = mmscale.fit_transform(np.array(sub_bhv_rew[:len(df_eeg)]).reshape(-1,1)).flatten()
##        sub_bhv_tot_rew = mmscale.fit_transform(np.array(sub_bhv_tot_rew[:len(df_eeg)]).reshape(-1,1)).flatten()
        
        df_eeg.insert(loc=4,column='n_from_leave',value=resp)
        df_eeg.insert(loc=5,column='reward',value=sub_bhv_rew)
        df_eeg.insert(loc=6,column='tot_reward',value=sub_bhv_tot_rew)
        df_eeg.insert(loc=7,column='rt',value=sub_bhv_rt[:len(df_eeg)])

        col_name = "reward"#n_from_leave
        df_eeg["stim"] = df_eeg[col_name]
        df_eeg.loc[df_eeg[col_name]==0.0,"stim"]="leave"
        df_eeg.loc[(df_eeg[col_name]>0.0) & (df_eeg[col_name]<=0.3),"stim"] = "low"##"late"
        df_eeg.loc[(df_eeg[col_name]>0.3) & (df_eeg[col_name]<=0.6),"stim"] = "med"##"mid"
        df_eeg.loc[(df_eeg[col_name]>0.6) & (df_eeg[col_name]<=1.0),"stim"] = "high"##"early"

##        df_eeg['max_cz'] = stdscale.fit_transform(df_eeg['max_cz'].to_numpy().reshape(-1,1)).flatten()
##        df_eeg['max_c3'] = stdscale.fit_transform(df_eeg['max_c3'].to_numpy().reshape(-1,1)).flatten()
##        df_eeg['max_c4'] = stdscale.fit_transform(df_eeg['max_c4'].to_numpy().reshape(-1,1)).flatten()

        df_eeg.to_csv(f"./l1_eeg_bhv_hddm/raw/sub_data/{sub}_hddm_norm.csv",index=False)
###############################################################################################################################

    fpaths = pathlib.Path().glob('./l1_eeg_bhv_hddm/raw/sub_data/*_hddm_norm.csv')
    df_norm = concat_csv_to_df(fpaths)
    df_norm.to_csv("./l1_eeg_bhv_hddm/raw/foraging_hddm_norm.csv",index=False)
