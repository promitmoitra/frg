%% Preamble:
%%% Load EEGLAB: DO NOT add to MATLAB path - causes problems with BioSig
%%% and add data path
%%% NOTE: eeglab2024.2.1 unable to load chanlocs from EDF.
clear;clc;
runtime = struct();
eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
wrk_dir = '/home/decision_lab/work/github/frg/';
dir_sep = '/';

% eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1\";
% wrk_dir = "C:\Users\promitmoitra\Documents\GitHub\frg\";
% dir_sep = '\';

global data_path; data_path = [char(wrk_dir),'data'];
%%
cd(eeglab_dir)
eeglab nogui;
cd(wrk_dir)
%%
badids = [37532 38058 39862 43543 45528 47801 48278];
%%%Marker issues with [42125,46037,47744,48238]
%%%Assuming each file is in a folder named <subid>
sub_dirs = dir(data_path);
subids = str2double({sub_dirs(~isnan(str2double({sub_dirs(:).name}))).name});
subids = setdiff(subids,badids);

% subids = [31730,47204,47131,48238,47324,43000];
% subids = [42949]
loop_start = tic;

for idx = 1:length(subids)
%%
    SAVEFLAG = 0;
    clearvars -except idx subids data_path loop_start runtime SAVEFLAG
    subid_start = tic;
    global subid; subid = subids(idx);
    fprintf([repmat('=',1,80),'\n',num2str(subid),'\n',repmat('=',1,80),'\n'])
    
    data_file = dir(fullfile(data_path,num2str(subid),'*.edf'));
    data_file_path = fullfile(data_file(end).folder,data_file(end).name);
    
    origchanlocs = readlocs([data_path,'/Statnet_F3F4FCz.ced']);
    if isfile([data_path '/raw/' num2str(subid) '.set'])
	    fprintf([num2str(subid) ' processed. Loading next...\n'])
	    continue
    else
    EEG = pop_biosig(data_file_path,'channels',1:19);
    EEG = pop_select(EEG,'rmchannel',{'EKG2'});
    
    %%%Implicitly setting FCz as online reference (during recording)
    EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
    EEG = pop_chanedit(EEG,'load',[data_path,'/Statnet_F3F4FCz.ced']);
%%%%%%%%%%%%%%%%

%% Check marker issues: NOTE: CHECK "TRIGGER EVENT I" IN BADIDS
    bad_event_idx = find(ismember({EEG.event.type},{'condition 20','condition 21'}));
    for ix = bad_event_idx
        try
            msg = {EEG.event(ix-5:ix+5).type};
        catch
            try
                msg = {EEG.event(ix:ix+5).type};
            catch
                msg = {EEG.event(ix-5:ix).type};
            end
        end
        if ~isempty(bad_event_idx)
            runtime.(['subid_' num2str(subid) '_bad_events']) = bad_event_idx;
            runtime.(['subid_' num2str(subid) '_bad_event_' num2str(ix)]) = msg;
        end
    end
%% Corrections and time stamps using marker:
% 
%     switch subid
%     case 43000
%         EEG.urevent(912).type = 'TRIGGER EVENT N';
%         EEG.event(912).type = 'TRIGGER EVENT N';
% 
%     case 42125
%         [EEG.urevent(contains({EEG.urevent.type},'condition 20')).type] = deal('TRIGGER EVENT Z');
%         [EEG.event(contains({EEG.event.type},'condition 20')).type] = deal('TRIGGER EVENT Z');
% 
%     case 46037
%         [EEG.urevent(contains({EEG.urevent.type},'condition 21')).type] = deal('TRIGGER EVENT Z');
%         [EEG.event(contains({EEG.event.type},'condition 21')).type] = deal('TRIGGER EVENT Z');
% 
%     case 47744
%         [EEG.urevent(contains({EEG.urevent.type},'condition 21')).type] = deal('TRIGGER EVENT Z');
%         [EEG.event(contains({EEG.event.type},'condition 21')).type] = deal('TRIGGER EVENT Z');
% 
%     case 48031
%         EEG.urevent(878).type = 'TRIGGER EVENT N';
%         EEG.event(878).type = 'TRIGGER EVENT N';
% 
%     case 48276
%         EEG.urevent(50).type = 'TRIGGER EVENT N';
%         EEG.event(50).type = 'TRIGGER EVENT N';
% 
%     case 48238
%         EEG.urevent(161).type = 'TRIGGER EVENT R';
%         EEG.event(161).type = 'TRIGGER EVENT R';
%     end
% 
%     all_events = {EEG.event.type};
%     if subid == 7873
%         pre_start = 1; pre_start_time=0;
%         pre_end = find(contains(all_events,'TRIGGER EVENT B'));
%     else
%         pre_start = find(ismember(all_events,'TRIGGER EVENT A')); pre_start = pre_start(1);
%         pre_start_time = EEG.event(pre_start(1)).latency/EEG.srate;
%         pre_end = find(contains(all_events(pre_start:end),'TRIGGER EVENT Z'));pre_end = pre_end(1)+pre_start;
%     end
% 
%     pre_end_time = EEG.event(pre_end).latency/EEG.srate;
%     pre_timestamp = [pre_start_time,pre_end_time];
%     
%     pre_game_start = find(strcmp(all_events,'TRIGGER EVENT B'));
%     pre_game_start_time = EEG.event(pre_game_start).latency/EEG.srate;
%     
%     pre_game_end = find(strcmp(all_events(pre_game_start:end),'TRIGGER EVENT Z'));
%     pre_game_end = pre_game_end(1)+pre_game_start;
%     pre_game_end_time = EEG.event(pre_game_end).latency/EEG.srate;
%     
%     pre_game_events = all_events(pre_game_start:pre_game_end);
%     
%     pre_game_short_end = find(strcmp(pre_game_events,'TRIGGER EVENT T'));
%     pre_game_short_end = pre_game_short_end(end);
%     pre_game_long_end = find(strcmp(pre_game_events,'TRIGGER EVENT U'));
%     pre_game_long_end = pre_game_long_end(end);
%     
%     pre_game_short_end_time = EEG.event(pre_game_start+pre_game_short_end-1).latency/EEG.srate;
%     pre_game_long_end_time = EEG.event(pre_game_start+pre_game_long_end-1).latency/EEG.srate;
%     
%     if pre_game_short_end<pre_game_long_end
%         pre_short_events = all_events(pre_game_start:pre_game_short_end);
%         pre_short_timestamp = [pre_game_start_time pre_game_short_end_time];
%         
%         pre_long_events = all_events(pre_game_short_end:pre_game_long_end);
%         pre_long_timestamp = [pre_game_short_end_time pre_game_long_end_time];
%     else
%         pre_long_events = all_events(pre_game_start:pre_game_long_end);
%         pre_long_timestamp = [pre_game_start_time pre_game_long_end_time];
%         
%         pre_short_events = all_events(pre_game_long_end:pre_game_short_end);
%         pre_short_timestamp = [pre_game_long_end_time pre_game_short_end_time];
%     end
%     
%     post_game_start = find(strcmp(all_events,'TRIGGER EVENT I'));
%     post_game_start_time = EEG.event(post_game_start).latency/EEG.srate;
%     
%     post_game_end = find(strcmp(all_events(post_game_start:end),'TRIGGER EVENT Z'));
%     post_game_end = post_game_end(1)+post_game_start;
%     post_game_end_time = EEG.event(post_game_end).latency/EEG.srate;
%     
%     post_game_events = all_events(post_game_start:post_game_end);
%     
%     post_game_short_end = find(strcmp(post_game_events,'TRIGGER EVENT T'));
%     post_game_short_end = post_game_short_end(end);
%     post_game_long_end = find(strcmp(post_game_events,'TRIGGER EVENT U'));
%     post_game_long_end = post_game_long_end(end);
%     
%     post_game_short_end_time = EEG.event(post_game_start+post_game_short_end-1).latency/EEG.srate;
%     post_game_long_end_time = EEG.event(post_game_start+post_game_long_end-1).latency/EEG.srate;
%     
%     if post_game_short_end<post_game_long_end
%         post_short_events = all_events(post_game_start:post_game_short_end);
%         post_short_timestamp = [post_game_start_time post_game_short_end_time];
%         
%         post_long_events = all_events(post_game_short_end:post_game_long_end);
%         post_long_timestamp = [post_game_short_end_time post_game_long_end_time];
%     else
%         post_long_events = all_events(post_game_start:post_game_long_end);
%         post_long_timestamp = [post_game_start_time post_game_long_end_time];
%         
%         post_short_events = all_events(post_game_long_end:post_game_short_end);
%         post_short_timestamp = [post_game_long_end_time post_game_short_end_time];
%     end
% %%%%%%%%%%%%%%%%
% %%
%     %%%Rereferencing to mastoids
%     EEG = pop_reref(EEG,{'A1','A2'});
%     EEG = pop_interp(EEG,origchanlocs);
%     raw = EEG;
%     
%     %%%ZapLine at 60Hz, FIR bandpass [0.5 120] and ICA
%     thresh = [0 0; 0.9 1; 0.9 1; 0.9 1; 0.9 1; 0.9 1; 0 0];
%     EEG = frg_clean_filt_ica(EEG,thresh);
%     
%     %%%Rereferencing to average
%     raw = pop_reref(raw,[]);
%     EEG = pop_reref(EEG,[]);
%     if SAVEFLAG
%         raw = pop_saveset(raw,'filename',[num2str(subid) '.set'],...
%                               'filepath',[data_path '/raw/'],'savemode','onefile');
%         EEG = pop_saveset(EEG,'filename',[num2str(subid) '_ica.set'],...
%                               'filepath',[data_path '/ica/'],'savemode','onefile');
%     end
% %% OLD CORRECTIONS
% %     all_events = {EEG.event.type};
% %     if size(find(ismember(all_events,'TRIGGER EVENT A')),2)>0
% %         pre_start = find(ismember(all_events,'TRIGGER EVENT A')); pre_start = pre_start(1);
% %         pre_start_time = EEG.event(pre_start(1)).latency/EEG.srate;
% %     else
% %         pre_start = 1; pre_start_time=0;
% %         fprintf([repmat('=',1,80),'\n','NO TRIGGER EVENT A FOUND!!!\n',repmat('=',1,80),'\n'])
% %     end
% %     
% %     if subid==42125
% %         block_end_marker = 'condition 20';
% %     elseif find(ismember([46037,47744],subid))
% %         block_end_marker = 'condition 21';
% %     else
% %         block_end_marker = 'TRIGGER EVENT Z';
% %     end
% %     
% %     if subid==7873
% %         pre_end = find(strcmp(all_events,'TRIGGER EVENT B'));
% %     else
% %         pre_end = find(strcmp(all_events(pre_start:end),block_end_marker));
% %     end
% %     pre_end = pre_end(1)+pre_start;
% %     pre_end_time = EEG.event(pre_end).latency/EEG.srate;
% %     pre_timestamp = [pre_start_time,pre_end_time];
% %     
% %     pre_game_start = find(strcmp(all_events,'TRIGGER EVENT B'));
% %     pre_game_start_time = EEG.event(pre_game_start).latency/EEG.srate;
% %     
% %     pre_game_end = find(strcmp(all_events(pre_game_start:end),block_end_marker));
% %     pre_game_end = pre_game_end(1)+pre_game_start;
% %     pre_game_end_time = EEG.event(pre_game_end).latency/EEG.srate;
% %     
% %     pre_game_events = all_events(pre_game_start:pre_game_end);
% %     
% %     pre_game_short_end = find(strcmp(pre_game_events,'TRIGGER EVENT T'));
% %     pre_game_short_end = pre_game_short_end(end);
% %     pre_game_long_end = find(strcmp(pre_game_events,'TRIGGER EVENT U'));
% %     pre_game_long_end = pre_game_long_end(end);
% %     
% %     pre_game_short_end_time = EEG.event(pre_game_start+pre_game_short_end-1).latency/EEG.srate;
% %     pre_game_long_end_time = EEG.event(pre_game_start+pre_game_long_end-1).latency/EEG.srate;
% %     
% %     if pre_game_short_end<pre_game_long_end
% %         pre_short_events = all_events(pre_game_start:pre_game_short_end);
% %         pre_short_timestamp = [pre_game_start_time pre_game_short_end_time];
% %         
% %         pre_long_events = all_events(pre_game_short_end:pre_game_long_end);
% %         pre_long_timestamp = [pre_game_short_end_time pre_game_long_end_time];
% %     else
% %         pre_long_events = all_events(pre_game_start:pre_game_long_end);
% %         pre_long_timestamp = [pre_game_start_time pre_game_long_end_time];
% %         
% %         pre_short_events = all_events(pre_game_long_end:pre_game_short_end);
% %         pre_short_timestamp = [pre_game_long_end_time pre_game_short_end_time];
% %     end
% %     
% %     post_game_start = find(strcmp(all_events,'TRIGGER EVENT I'));
% %     post_game_start_time = EEG.event(post_game_start).latency/EEG.srate;
% %     
% %     post_game_end = find(strcmp(all_events(post_game_start:end),block_end_marker));
% %     post_game_end = post_game_end(1)+post_game_start;
% %     post_game_end_time = EEG.event(post_game_end).latency/EEG.srate;
% %     
% %     post_game_events = all_events(post_game_start:post_game_end);
% %     
% %     post_game_short_end = find(strcmp(post_game_events,'TRIGGER EVENT T'));
% %     post_game_short_end = post_game_short_end(end);
% %     post_game_long_end = find(strcmp(post_game_events,'TRIGGER EVENT U'));
% %     post_game_long_end = post_game_long_end(end);
% %     
% %     post_game_short_end_time = EEG.event(post_game_start+post_game_short_end-1).latency/EEG.srate;
% %     post_game_long_end_time = EEG.event(post_game_start+post_game_long_end-1).latency/EEG.srate;
% %     
% %     if post_game_short_end<post_game_long_end
% %         post_short_events = all_events(post_game_start:post_game_short_end);
% %         post_short_timestamp = [post_game_start_time post_game_short_end_time];
% %         
% %         post_long_events = all_events(post_game_short_end:post_game_long_end);
% %         post_long_timestamp = [post_game_short_end_time post_game_long_end_time];
% %     else
% %         post_long_events = all_events(post_game_start:post_game_long_end);
% %         post_long_timestamp = [post_game_start_time post_game_long_end_time];
% %         
% %         post_short_events = all_events(post_game_long_end:post_game_short_end);
% %         post_short_timestamp = [post_game_long_end_time post_game_short_end_time];
% %     end
% %%%%%%%%%%%%%%%%
% %% Segment and save raw
%     preshort_raw = pop_select(raw,'time',pre_short_timestamp);
%     prelong_raw = pop_select(raw,'time',pre_long_timestamp);
%     postshort_raw = pop_select(raw,'time',post_short_timestamp);
%     postlong_raw = pop_select(raw,'time',post_long_timestamp);
%     if SAVEFLAG
%         preshort_raw = pop_saveset(preshort_raw,'filename',[num2str(subid) '_preshort.set'],...
%                                   'filepath',[data_path '/raw_seg/'],'savemode','onefile');
%         prelong_raw = pop_saveset(prelong_raw,'filename',[num2str(subid) '_prelong.set'],...
%                                     'filepath',[data_path '/raw_seg/'],'savemode','onefile');
%         postshort_raw = pop_saveset(postshort_raw,'filename',[num2str(subid) '_postshort.set'],...
%                                     'filepath',[data_path '/raw_seg/'],'savemode','onefile');
%         postlong_raw = pop_saveset(postlong_raw,'filename',[num2str(subid) '_postlong.set'],...
%                                     'filepath',[data_path '/raw_seg/'],'savemode','onefile');
%     end
% %%%%%%%%%%%%%%%%
% %% Segment and save clean_filt_ica
%     preshort_EEG = pop_select(EEG,'time',pre_short_timestamp);
%     prelong_EEG = pop_select(EEG,'time',pre_long_timestamp);
%     postshort_EEG = pop_select(EEG,'time',post_short_timestamp);
%     postlong_EEG = pop_select(EEG,'time',post_long_timestamp);
% %     SAVEFLAG=0;
%     if SAVEFLAG
%         preshort_EEG = pop_saveset(preshort_EEG,'filename',[num2str(subid) '_preshort.set'],...
%                                     'filepath',[data_path '/ica_seg/'],'savemode','onefile');
%         prelong_EEG = pop_saveset(prelong_EEG,'filename',[num2str(subid) '_prelong.set'],...
%                                     'filepath',[data_path '/ica_seg/'],'savemode','onefile');
%         postshort_EEG = pop_saveset(postshort_EEG,'filename',[num2str(subid) '_postshort.set'],...
%                                     'filepath',[data_path '/ica_seg/'],'savemode','onefile');
%         postlong_EEG = pop_saveset(postlong_EEG,'filename',[num2str(subid) '_postlong.set'],...
%                                     'filepath',[data_path '/ica_seg/'],'savemode','onefile');
%     end
% %%%%%%%%%%%%%%%%
% %%  Epoch and save early, mid, late and leave
% %     SAVEFLAG = 0;
%     [epoch_data,early,mid,late,leave] = frg_epoch(preshort_EEG,'TRIGGER EVENT N',[-1 2],SAVEFLAG);
%     [epoch_data,early,mid,late,leave] = frg_epoch(prelong_EEG,'TRIGGER EVENT N',[-1 2],SAVEFLAG);
%     [epoch_data,early,mid,late,leave] = frg_epoch(postshort_EEG,'TRIGGER EVENT N',[-1 2],SAVEFLAG);
%     [epoch_data,early,mid,late,leave] = frg_epoch(postlong_EEG,'TRIGGER EVENT N',[-1 2],SAVEFLAG);
% %%%%%%%%%%%%%%%%
% 
% %%
%     subid_end = toc(subid_start);
%     runtime.(['subid_' num2str(subid)]) = subid_end;
    end
end    
loop_end = toc(loop_start);
runtime.loop = loop_end;
% save('frg_marker_issues.mat','-struct',"runtime");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%==Function Definitions==%%%%%%%%%%%%%%%%%%%%%%%%%

function eeg = frg_clean_filt_ica(eeg,thresh)
    eeg = clean_data_with_zapline_plus_eeglab_wrapper(eeg,struct('noisefreqs',60, ...
                                       'chunkLength',0,'adaptiveNremove',true, ...
                                       'fixedNremove',1,'plotResults',0));

    freqs = [2 120];
    wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,eeg.srate,df);
    eeg_ica = pop_firws(eeg,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);
    
    freqs = [0.5 120];
    wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,eeg.srate,df);
    eeg = pop_firws(eeg,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);

    eeg_fields = fieldnames(eeg);
    eeg_icafields = eeg_fields(contains(fieldnames(eeg),'ica'));
    if sum(cellfun(@(f) isempty(eeg.(f)),eeg_icafields))
%         eeg_ica = pop_runica(eeg_ica,'icatype','runica','extended',1,'interrupt','off');
        eeg_ica = pop_runica(eeg_ica,'icatype','fastica');
        for idx=1:size(eeg_icafields,1)
            eeg.(char(eeg_icafields(idx))) = eeg_ica.(char(eeg_icafields(idx)));
        end
    end
%            [ Brain,Muscle, Eye, Heart, Line,Channel,Other ]
%     thresh = [0 0; 0.67 1; 0.67 1; 0 0; 0.8 1; 0.8 1; 0 0];
    eeg = pop_iclabel(eeg,'Default');
    eeg = pop_icflag(eeg,thresh);
    rejected_comps = find(eeg.reject.gcompreject > 0);
    eeg = pop_subcomp(eeg,rejected_comps);
end

function [ep_dat,early,mid,late,leave] = frg_epoch(eeg_dat,event_str,epoch_lims,saveflag)
    global subid;
    state_cond = inputname(1); prefix = state_cond(1:strfind(state_cond,'_')-1);

    %Explain how leave indexes are directly found from urevents
    leave_S_uridxs = [eeg_dat.event(find(contains({eeg_dat.event.type},'TRIGGER EVENT S'))).urevent];
    leave_R_uridxs = [eeg_dat.event(find(contains({eeg_dat.event.type},'TRIGGER EVENT R'))).urevent];

    if length(leave_R_uridxs)>length(leave_S_uridxs)
        leave_N_uridxs = leave_R_uridxs-1;
    else
        leave_N_uridxs = leave_S_uridxs-2;
    end

    ep_dat = pop_epoch(eeg_dat,{event_str},epoch_lims);
    ep_trial_uridxs = cellfun(@cell2mat, {ep_dat.epoch(:).eventurevent},'UniformOutput',false);
    ep_leave = cellfun(@(x) ismember(x,leave_N_uridxs),ep_trial_uridxs,'UniformOutput',false);
    leave_idxs = find(cellfun(@(x) sum(x),ep_leave));

    %OLD AND BUGGY LEAVE IDX FINDER - not all epochs contain R
%     %trial_events = cellfun(@cell2mat,{ep_dat.epoch(:).eventtype},'UniformOutput',false);
%     %leave_idxs = find(contains(trial_events,'TRIGGER EVENT R'));

    entry_idxs = [1 leave_idxs(1:end-1)+1];
    patch_trials = arrayfun(@(f,g) (f:g),entry_idxs,leave_idxs,'UniformOutput',false);
    patch_len = cell2mat(cellfun(@length,patch_trials,'UniformOutput',false));
    
    %Explain how these lines categorize patch trials as early, mid, late and leave
    patch_trial_cat = arrayfun(@(f) discretize((1:f)/f,[0 0.34 0.67 0.99 1]),...
                                patch_len,'UniformOutput',false);
    patch_trial_cat([find(patch_len==2)]) = {[1 4]};
    patch_trial_cat([find(patch_len==3)]) = {[1 3 4]};
    flat_cat = cell2mat(patch_trial_cat);
    flat_cat = [flat_cat zeros(1,ep_dat.trials-size(flat_cat,2))];
    flat_cat = mat2cell(flat_cat,size(flat_cat,1),ones(1,size(flat_cat,2)));
    [ep_dat.epoch.cat_latency] = flat_cat{:};
    try
        early = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==1));
    catch
        early = struct();
    end

    try
        mid = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==2));
    catch
        mid = struct();
    end

    try
        late = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==3));
    catch
        late = struct();
    end

    try
        leave = pop_select(ep_dat,'trial',find([ep_dat.epoch(:).cat_latency]==4));
    catch
        leave = struct();
    end

    if saveflag
        if ~isempty(fieldnames(early))
            early = pop_saveset(early,'filename',[num2str(subid) '_' prefix '_early.set'],...
                                'filepath',[data_path '/epoch_data/'],'savemode','onefile');
        end
        if ~isempty(fieldnames(mid))
            mid = pop_saveset(mid,'filename',[num2str(subid) '_' prefix '_mid.set'],...
                                'filepath',[data_path '/epoch_data/'],'savemode','onefile');
        end
        if ~isempty(fieldnames(late))
            late = pop_saveset(late,'filename',[num2str(subid) '_' prefix '_late.set'],...
                                'filepath',[data_path '/epoch_data/'],'savemode','onefile');
        end
        if ~isempty(fieldnames(leave))
            leave = pop_saveset(leave,'filename',[num2str(subid) '_' prefix '_leave.set'],...
                                'filepath',[data_path '/epoch_data/'],'savemode','onefile');
        end
    end
end