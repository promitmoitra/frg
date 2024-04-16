import numpy as np
import scipy.io as sio
import csv
import pathlib
import pandas as pd
'''
Next attempt:   Read mat file, check pre/poststress["CodeNumbers"][0,trial_number] for 3:'M', 4:'N', 5:'O', 8:'R' (roughly!)
                Extract trial numbers and response from above, to compare with df_eeg
Update:         Done!!
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

def insert_row(df,row_num,row_val):
    df_top = df[0:row_num].copy()
    df_bottom = df[row_num:].copy()
    df_top.loc[row_num] = row_val
    df_out = pd.concat([df_top,df_bottom])
    df_out.index = [*range(df_out.shape[0])]
    return df_out

def concat_csv_to_df(filepaths):
    df = pd.DataFrame()
    for f in filepaths:
        tmp = pd.read_csv(f)
        df = pd.concat([df,tmp])
    return df

if __name__=="__main__":
    data_path = '../Mrugsen/foraging/Neuroflow/'
    all_ids = [int(i.stem) for i in pathlib.Path().glob(data_path+'[0-9]*')]
    all_ids.sort()
    bad_ids = [37532, 38058, 39862, 42125, 43543, 46037, 47744, 47801, 48238, 48278]
    sub_ids = np.setdiff1d(all_ids,bad_ids)

    bhv_response = []
    bhv_rt = []
    df_comp = pd.read_csv("nflw_eeg_comp.csv")

    for sub in sub_ids:
        sub_start_idx = list(df_comp[df_comp["subject"]==sub].index)[0]
        sub_len = len(df_comp[df_comp["subject"]==sub])

        fpath = list(pathlib.Path().glob(data_path+str(sub)+'/*?oraging*/'+str(sub)+'.mat'))[0]
        pre_bhv = read_bhv_mat(fpath,pre=True)
        post_bhv = read_bhv_mat(fpath,pre=False)

        fpath = list(pathlib.Path().glob(str(sub)+'_hddm.csv'))[0]
        df_eeg = pd.read_csv(fpath)
        eeg_pre_num_trial = len(df_eeg[(df_eeg["subj_idx"]==sub) & (df_eeg["condition"]==1)])
        eeg_post_num_trial = len(df_eeg[(df_eeg["subj_idx"]==sub) & (df_eeg["condition"]==2)])
        eeg_num_trial = eeg_pre_num_trial + eeg_post_num_trial

##############################################################################################################################

        """calculate num_trials from bhv raw data"""
##        bhv_pre_num_trial = 0  
##        for idx in range(len(pre_bhv["CodeNumbers"][0])):
##            if 3 in pre_bhv["CodeNumbers"][0,idx].flatten():
##                bhv_pre_num_trial += 1
##
##        bhv_post_num_trial = 0
##        for idx in range(len(post_bhv["CodeNumbers"][0])):
##            if 3 in post_bhv["CodeNumbers"][0,idx].flatten():
##                bhv_post_num_trial += 1
##
##        bhv_num_trial = bhv_pre_num_trial + bhv_post_num_trial

##        print(sub,": ")
##        print("pre_bhv: ",bhv_pre_num_trial)
##        print("pre_eeg: ",eeg_pre_num_trial)
##        print("post_bhv: ",bhv_post_num_trial)
##        print("post_eeg: ",eeg_post_num_trial)

        """calculate response and rt from bhv trials"""
        bhv_pre_response = []
        bhv_pre_rt = []
        for idx in range(len(pre_bhv["CodeNumbers"][0])):
            codes = pre_bhv["CodeNumbers"][0,idx].flatten()
            times = pre_bhv["CodeTimes"][0,idx].flatten()
            if 5 in codes:
                bhv_pre_response.append(1)
                rt = times[codes==5]-times[codes==4]
                bhv_pre_rt.append(rt[0])
            elif 8 in codes:
                bhv_pre_response.append(0)
                rt = times[codes==8]-times[codes==4]
                bhv_pre_rt.append(rt[0])

        bhv_post_response = []
        bhv_post_rt = []
        for idx in range(len(post_bhv["CodeNumbers"][0])):
            codes = post_bhv["CodeNumbers"][0,idx].flatten()
            times = post_bhv["CodeTimes"][0,idx].flatten()
            if 5 in codes:
                bhv_post_response.append(1)
                rt = times[codes==5]-times[codes==4]
                bhv_post_rt.append(rt[0])
            elif 8 in codes:
                bhv_post_response.append(0)
                rt = times[codes==8]-times[codes==4]
                bhv_post_rt.append(rt[0])

        sub_bhv_response = bhv_pre_response + bhv_post_response
        sub_bhv_rt = bhv_pre_rt + bhv_post_rt
        sub_bhv_rt = np.array(sub_bhv_rt)/1000
        sub_bhv_rt = list(sub_bhv_rt)

##      adjust loaded df_comp to acco new columns
        len_diff = len(sub_bhv_response) - sub_len 
        if len_diff<0:
            sub_bhv_response += [np.nan]*len_diff
            sub_bhv_rt += [np.nan]*len_diff
        elif len_diff>0:
            for i in range(len_diff):
                new_sub_len = len(df_comp[df_comp["subject"]==sub])
                df_comp = insert_row(df_comp,sub_start_idx+new_sub_len,[np.nan]*4+[sub])
        print(len(df_comp[df_comp["subject"]==sub]),len(sub_bhv_response))
        bhv_response += sub_bhv_response
        bhv_rt += sub_bhv_rt

##      compensate for missing sub 31730
##    sub_start_idx = list(df_comp[df_comp["subject"]==31730].index)[0]
##    sub_len = len(df_comp[df_comp["subject"]==31730])
##    bhv_response[sub_start_idx:sub_start_idx] = [np.nan]*sub_len
##    bhv_rt[sub_start_idx:sub_start_idx] = [np.nan]*sub_len
        
    df_comp.insert(loc=2, column='response_bhv', value=bhv_response)
    df_comp.insert(loc=5, column='rt_bhv', value=bhv_rt)
##    df_comp.to_csv("nflw_eeg_bhv_2.csv")
###############################################################################################################################

##        bhv_pre_num_trial = sum(pre_bhv['TrialNumber'].flatten() == np.array(range(len(pre_bhv['TrialNumber'].flatten())))+1)
##        bhv_post_num_trial = sum(post_bhv['TrialNumber'].flatten() == np.array(range(len(post_bhv['TrialNumber'].flatten())))+1)

##        df_comp = pd.read_csv("nflw_eeg_comp.csv")
##        nflw_num_trial = sum(np.logical_not(pd.isna(df_comp[df_comp["subject"]==sub]["response_nflw"])))
##        eeg_num_trial = sum(np.logical_not(pd.isna(df_comp[df_comp["subject"]==sub]["response_eeg"])))
##        print("bhv_raw:\t",bhv_num_trial," trials")
##        print("bhv_nflw:\t",nflw_num_trial," trials")
##        print("eeg:\t\t",eeg_num_trial," trials\n")
##        print("Mismatch: ",bhv_num_trial-eeg_num_trial,"\n")
