%% Preamble:
%%% Load EEGLAB: DO NOT add to MATLAB path - causes problems with BioSig
%%% and add data path
%%% NOTE: eeglab2024.2.1 unable to load chanlocs from
clear;clc;
% rmpath_pat('NBTpublic')
eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
wrk_dir = '/home/decision_lab/work/github/frg/';
dir_sep = '/';

% eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1\";
% wrk_dir = "C:\Users\promitmoitra\Documents\GitHub\frg\";
% dir_sep = '\';

data_path = [char(wrk_dir),'data'];
%%
cd(eeglab_dir)
eeglab nogui;
cd(wrk_dir)
%%
% badids = [37532 38058 39862 43543 45528 47801 48278];
% %%%Assuming each file is in a folder named <subid>
% edf_fnames = dir('**/*.edf');
% data_path_dirs = split({edf_fnames(:).folder},dir_sep);
% subids = str2num(char(data_path_dirs(:,:,end)));
% subids = setdiff(subids,badids);

subids = [31730];%,47204,47131,48238,47324,43000];
%%
for idx = 1:length(subids)
%%    
    clearvars -except idx subids data_path %data_path_dirs
    subid = subids(idx);
    fprintf([repmat('=',1,80),'\n',char(string(subid)),'\n',repmat('=',1,80),'\n'])
    
    data_file = dir(fullfile(data_path,string(subid),'*.edf'));
    data_file_path = fullfile(data_file(end).folder,data_file(end).name);
    
    origchanlocs = readlocs([data_path,'/Statnet_F3F4FCz.ced']);
    
    EEG = pop_biosig(data_file_path,'channels',1:19);
    EEG = pop_select(EEG,'rmchannel',{'EKG2'});
    
    EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
    EEG = pop_chanedit(EEG,'load',[data_path,'/Statnet_F3F4FCz.ced']);
    EEG = pop_reref(EEG,{'A1','A2'});
    EEG = pop_interp(EEG,origchanlocs);
    EEG = pop_reref(EEG,[]);

    EEG = frg_bandpass_clean(EEG);
%%
    all_events = {EEG.event.type};
    if size(find(ismember(all_events,'TRIGGER EVENT A')),2)>0
        pre_start = find(ismember(all_events,'TRIGGER EVENT A')); pre_start = pre_start(1);
        pre_start_time = EEG.event(pre_start(1)).latency/EEG.srate;
    else
        pre_start = 1; pre_start_time=0;
        fprintf([repmat('=',1,80),'\n','NO TRIGGER EVENT A FOUND!!!\n',repmat('=',1,80),'\n'])
    end
    
    if subid==42125
        block_end_marker = 'condition 20';
    elseif find(ismember([46037,47744],subid))
        block_end_marker = 'condition 21';
    else
        block_end_marker = 'TRIGGER EVENT Z';
    end
    
    if subid==7873
        pre_end = find(strcmp(all_events,'TRIGGER EVENT B'));
    else
        pre_end = find(strcmp(all_events(pre_start:end),block_end_marker));
    end
    pre_end = pre_end(1)+pre_start;
    pre_end_time = EEG.event(pre_end).latency/EEG.srate;
    pre_timestamp = [pre_start_time,pre_end_time];
    
    pre_game_start = find(strcmp(all_events,'TRIGGER EVENT B'));
    pre_game_start_time = EEG.event(pre_game_start).latency/EEG.srate;
    
    pre_game_end = find(strcmp(all_events(pre_game_start:end),block_end_marker));
    pre_game_end = pre_game_end(1)+pre_game_start;
    pre_game_end_time = EEG.event(pre_game_end).latency/EEG.srate;
    
    pre_game_events = all_events(pre_game_start:pre_game_end);
    
    pre_game_short_end = find(strcmp(pre_game_events,'TRIGGER EVENT T'));
    pre_game_short_end = pre_game_short_end(end);
    pre_game_long_end = find(strcmp(pre_game_events,'TRIGGER EVENT U'));
    pre_game_long_end = pre_game_long_end(end);
    
    pre_game_short_end_time = EEG.event(pre_game_start+pre_game_short_end-1).latency/EEG.srate;
    pre_game_long_end_time = EEG.event(pre_game_start+pre_game_long_end-1).latency/EEG.srate;
    
    if pre_game_short_end<pre_game_long_end
        pre_short_events = all_events(pre_game_start:pre_game_short_end);
        pre_short_timestamp = [pre_game_start_time pre_game_short_end_time];
        
        pre_long_events = all_events(pre_game_short_end:pre_game_long_end);
        pre_long_timestamp = [pre_game_short_end_time pre_game_long_end_time];
    else
        pre_long_events = all_events(pre_game_start:pre_game_long_end);
        pre_long_timestamp = [pre_game_start_time pre_game_long_end_time];
        
        pre_short_events = all_events(pre_game_long_end:pre_game_short_end);
        pre_short_timestamp = [pre_game_long_end_time pre_game_short_end_time];
    end
    
    post_game_start = find(strcmp(all_events,'TRIGGER EVENT I'));
    post_game_start_time = EEG.event(post_game_start).latency/EEG.srate;
    
    post_game_end = find(strcmp(all_events(post_game_start:end),block_end_marker));
    post_game_end = post_game_end(1)+post_game_start;
    post_game_end_time = EEG.event(post_game_end).latency/EEG.srate;
    
    post_game_events = all_events(post_game_start:post_game_end);
    
    post_game_short_end = find(strcmp(post_game_events,'TRIGGER EVENT T'));
    post_game_short_end = post_game_short_end(end);
    post_game_long_end = find(strcmp(post_game_events,'TRIGGER EVENT U'));
    post_game_long_end = post_game_long_end(end);
    
    post_game_short_end_time = EEG.event(post_game_start+post_game_short_end-1).latency/EEG.srate;
    post_game_long_end_time = EEG.event(post_game_start+post_game_long_end-1).latency/EEG.srate;
    
    if post_game_short_end<post_game_long_end
        post_short_events = all_events(post_game_start:post_game_short_end);
        post_short_timestamp = [post_game_start_time post_game_short_end_time];
        
        post_long_events = all_events(post_game_short_end:post_game_long_end);
        post_long_timestamp = [post_game_short_end_time post_game_long_end_time];
    else
        post_long_events = all_events(post_game_start:post_game_long_end);
        post_long_timestamp = [post_game_start_time post_game_long_end_time];
        
        post_short_events = all_events(post_game_long_end:post_game_short_end);
        post_short_timestamp = [post_game_long_end_time post_game_short_end_time];
    end

    preshort_EEG = pop_select(EEG,'time',pre_short_timestamp);
    preshort_EEG = frg_bandpass_clean(preshort_EEG);
    preshort_EEG = pop_saveset(preshort_EEG,'filename',[num2str(subid) '_preshort.set'],...
                                'filepath',[data_path '/' num2str(subid)],'savemode','onefile');

    prelong_EEG = pop_select(EEG,'time',pre_long_timestamp);
    prelong_EEG = frg_bandpass_clean(prelong_EEG);
    prelong_EEG = pop_saveset(prelong_EEG,'filename',[num2str(subid) '_prelong.set'],...
                                'filepath',[data_path '/' num2str(subid)],'savemode','onefile');

    postshort_EEG = pop_select(EEG,'time',post_short_timestamp);
    postshort_EEG = frg_bandpass_clean(postshort_EEG);
    postshort_EEG = pop_saveset(postshort_EEG,'filename',[num2str(subid) '_postshort.set'],...
                                'filepath',[data_path '/' num2str(subid)],'savemode','onefile');

    postlong_EEG = pop_select(EEG,'time',post_long_timestamp);
    postlong_EEG = frg_bandpass_clean(postlong_EEG);
    postlong_EEG = pop_saveset(postlong_EEG,'filename',[num2str(subid) '_postlong.set'],...
                                'filepath',[data_path '/' num2str(subid)],'savemode','onefile');

%%
    [epoch_data,early,mid,late,leave] = frg_epoch(preshort_EEG,'TRIGGER EVENT N',[-1 2],1);
    [epoch_data,early,mid,late,leave] = frg_epoch(prelong_EEG,'TRIGGER EVENT N',[-1 2],1);
    [epoch_data,early,mid,late,leave] = frg_epoch(postshort_EEG,'TRIGGER EVENT N',[-1 2],1);
    [epoch_data,early,mid,late,leave] = frg_epoch(postlong_EEG,'TRIGGER EVENT N',[-1 2],1);

%%
% early = pop_loadset('filename', '31730_preshort_early_FZ-CZ_theta-gamma.set','filepath','/home/decision_lab/work/github/frg/data/pac_data')
% mid = pop_loadset('filename', '31730_preshort_mid_FZ-CZ_theta-gamma.set','filepath','/home/decision_lab/work/github/frg/data/pac_data')
% late = pop_loadset('filename', '31730_preshort_late_FZ-CZ_theta-gamma.set','filepath','/home/decision_lab/work/github/frg/data/pac_data')
% leave = pop_loadset('filename', '31730_preshort_leave_FZ-CZ_theta-gamma.set','filepath','/home/decision_lab/work/github/frg/data/pac_data')
% 
%     phase_chan_idx = find(ismember({epoch_data.chanlocs.labels},{'FZ'}));
%     amp_chan_idx = find(ismember({epoch_data.chanlocs.labels},{'CZ'}));
%     pac_dir = [data_path '/pac_data/'];
% 
%     tic
%     early = pop_pac(early,'Channels',[4 8],[30 90],[phase_chan_idx],[amp_chan_idx],...
%                     'method','ermipac','nboot',200,'alpha',[],...
%                     'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
% %     early = pop_saveset(early,'filename',[num2str(subid) '_preshort_early_FZ-CZ_theta-gamma.set'],...
% %                         'filepath',pac_dir,'savemode','onefile');
%     
%     mid = pop_pac(mid,'Channels',[4 8],[30 90],[phase_chan_idx],[amp_chan_idx],...
%                     'method','ermipac','nboot',200,'alpha',[],...
%                     'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
% %     mid = pop_saveset(mid,'filename',[num2str(subid) '_preshort_mid_FZ-CZ_theta-gamma.set'],...
% %                         'filepath',pac_dir,'savemode','onefile');
%     
%     late = pop_pac(late,'Channels',[4 8],[30 90],[phase_chan_idx],[amp_chan_idx],...
%                     'method','ermipac','nboot',200,'alpha',[],...
%                     'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
% %     late = pop_saveset(late,'filename',[num2str(subid) '_preshort_late_FZ-CZ_theta-gamma.set'],...
% %                         'filepath',pac_dir,'savemode','onefile');
%     
%     leave = pop_pac(leave,'Channels',[4 8],[30 90],[phase_chan_idx],[amp_chan_idx],...
%                     'method','ermipac','nboot',200,'alpha',[],...
%                     'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
% %     leave = pop_saveset(leave,'filename',[num2str(subid) '_preshort_leave_FZ-CZ_theta-gamma.set'],...
% %                         'filepath',pac_dir,'savemode','onefile');
%     toc
    
end

function eeg = frg_bandpass_clean(eeg)
    freqs = [0.5 120];
    wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,eeg.srate,df);
    eeg = pop_firws(eeg,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);
    
    eeg = clean_data_with_zapline_plus_eeglab_wrapper(eeg,struct('noisefreqs',60, ...
                                       'chunkLength',0,'adaptiveNremove',true, ...
                                       'fixedNremove',1,'plotResults',0));
end

function [ep_dat,early,mid,late,leave] = frg_epoch(raw_eeg,event_str,epoch_lims,saveflag)
    ep_dat = pop_epoch(raw_eeg,{event_str},epoch_lims);
    trigger_event = {'TRIGGER EVENT R'};

    %Explain how these two lines find leave indexes
    trial_events = cellfun(@cell2mat,{ep_dat.epoch(:).eventtype},'UniformOutput',false);
    leave_idxs = find(contains(trial_events,trigger_event));
    
    entry_idxs = [1 leave_idxs(1:end-1)+1];
    patch_trials = arrayfun(@(f,g) (f:g),entry_idxs,leave_idxs,'UniformOutput',false);
    patch_len = cell2mat(cellfun(@length,patch_trials,'UniformOutput',false));

    %Explain how these lines categorize patch trials as early, mid, late and leave
    patch_trial_cat = arrayfun(@(f) discretize((1:f)/f,[0 0.34 0.67 0.99 1]),patch_len,'UniformOutput',false);
    flat_cat = cell2mat(patch_trial_cat);
    flat_cat = [flat_cat zeros(1,ep_dat.trials-size(flat_cat,2))];
    flat_cat = mat2cell(flat_cat,size(flat_cat,1),ones(1,size(flat_cat,2)));
    [ep_dat.epoch.cat_latency] = flat_cat{:};
    early = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==1));
    mid = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==2));
    late = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==3));
    leave = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==4));
    if saveflag
        early = pop_saveset(early,'filename',[num2str(subid) '_preshort_early.set'],...
                            'filepath',[data_path '/' num2str(subid)],'savemode','onefile');
        mid = pop_saveset(mid,'filename',[num2str(subid) '_preshort_mid.set'],...
                            'filepath',[data_path '/' num2str(subid)],'savemode','onefile');
        late = pop_saveset(late,'filename',[num2str(subid) '_preshort_late.set'],...
                            'filepath',[data_path '/' num2str(subid)],'savemode','onefile');
        leave = pop_saveset(leave,'filename',[num2str(subid) '_preshort_leave.set'],...
                            'filepath',[data_path '/' num2str(subid)],'savemode','onefile');
    end
end