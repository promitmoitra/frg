%% Notes:
% +. Expand timestamp extraction section (224-284) to all behavioral variables.

%% Preamble:
%%% Load EEGLAB: DO NOT add to MATLAB path - causes problems with BioSig
%%% and add data path
clear;clc;
current_dir = pwd;
% cd "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1";
cd '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
eeglab;close;
cd(current_dir)
data_path = strcat(current_dir,'/data/');
badids = [37532 38058 39862 43543 45528 47801 48278];
processed_ids = dir([data_path,'spec_data/*.mat']);
proc_fnames = char({processed_ids(:).name});
for idx = 1:size(proc_fnames,1)
    proc_id = char(regexp(proc_fnames(idx,:),'[0-9]','match'))';
    badids(end+1) = str2num(proc_id);
end    
subids = readmatrix(fullfile(data_path,'subids.txt')); subids = setdiff(subids,badids);

%%%%%-New code, no need to read in subids.txt-%%%%%
% edf_fnames = dir('**/*.edf');
% data_path_dirs = split({edf_fnames(:).folder},'/');
% subids = str2num(char(data_path_dirs(:,:,end))); subids = setdiff(subids,badids);
%%%%%%%%%%%%%%%%%%%%
% bp_table = table();

% subids = [7873];
% subid = 31730;
%% Main loop begins:

%%%Comment out for single sub: Set up loop vars
for idx = 1:length(subids)
clearvars -global -except idx subids data_path %data_path_dirs %bp_table
global subid; subid = subids(idx);

fprintf([repmat('=',1,80),'\n',char(string(subid)),'\n',repmat('=',1,80)])

%%% Load data and channel locations:
data_file = dir(fullfile(data_path,string(subid),'*.edf'));
data_file_path = fullfile(data_file(end).folder,data_file(end).name);
% f_idx = find(ismember(data_path_dirs(:,:,end),num2str(subid)));
% data_file_path = fullfile(edf_fnames(f_idx).folder,edf_fnames(f_idx).name)

origchanlocs = readlocs('./data/Statnet_F3F4FCz.ced');

%%Comment out for single sub: Check existing set file
% if isfile([pwd '/data/clean_set/' num2str(subid) '_preproc.set'])
%     disp([num2str(subid) ': Already preprocessed'])
%     continue
% else    
%     disp(['Preprocessing:' num2str(subid)]);
% end    

EEG = pop_biosig(data_file_path,'channels',1:19);
EEG = pop_select(EEG,'rmchannel',{'EKG2'});

%Online reference was FCz, which is not there in the data. Adding an
%empty channel at index 19 (same as chanloc file) before loading the channel locations.
EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
EEG = pop_chanedit(EEG,'load','./data/Statnet_F3F4FCz.ced');

% fcz_idx = find(ismember({origchanlocs(:).labels},'FCz'));
% fcz = origchanlocs(fcz_idx); fcz.idx = fcz_idx;
% fcz = orderfields(fcz,[find(ismember(fieldnames(fcz),{'idx'})),1:size(fieldnames(fcz),1)-1]);
% fcz.a=0;fcz.b=0;
% EEG = pop_chanedit(EEG,'add',struct2cell(fcz));%,'load','./data/Statnet_F3F4FCz.ced','setref',{1:EEG.nbchan,fcz.labels});

%Re-referencing to mastoids A1 A2, interpolating FCz (17 in data, 19 in chanloc), and
%re-referencing to averef
EEG = pop_reref(EEG,{'A1','A2'});
EEG = pop_interp(EEG,origchanlocs);
EEG = pop_reref(EEG,[]);

%%
% freqs = [0.5 120];
% wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,EEG.srate,df);
% EEG = pop_firws(EEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);
% 
% EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('noisefreqs',60, ...
%                                    'chunkLength',0,'adaptiveNremove',true, ...
%                                    'fixedNremove',1,'plotResults',0));

% EEG_asr = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.1, ...
%     'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
%     'WindowCriterion',0.5,'BurstRejection','off','Distance','Euclidian', ...
%     'WindowCriterionTolerances',[-Inf 10]);
% 
% EEG_asr = pop_interp(EEG_asr,origchanlocs);
%
% % EEG_asr.etc=EEG.etc;
% % vis_artifacts(EEG_asr,EEG)
% %
%%%%%%%%
% [seg_delta_power,seg_theta_power,seg_alpha_power,seg_lbeta_power,...
%     seg_hbeta_power,seg_lgamma_power] = seg_chan_bp_abs(EEG,[0 300],sel_chans)
% [eeg,specdata] = seg_chan_bp_specpar(EEG,[0 300],sel_chans);
% figure;hold on
% ax=gca;ax.XScale='log';ax.YScale='log';
% mark = ['r','g','b','k'];
% 
% f_hz = 1:0.25:126;
% delta = find(f_hz>0 & f_hz<=2);
% theta = find(f_hz>4 & f_hz<=7);
% alpha = find(f_hz>7 & f_hz<=13);
% lbeta = find(f_hz>13 & f_hz<=20);
% hbeta = find(f_hz>20 & f_hz<=30);
% lgamma= find(f_hz>30 & f_hz<=45);
% 
% band = alpha;%[delta,theta,alpha,lbeta,hbeta,lgamma];
% for i = 1:4
%     loglog(f_hz(band),specdata(i,band),mark(i),'DisplayName',eeg.chanlocs(i).labels)
% %     P = polyfit(log(f_hz(band)),log(specdata(i,band)),1);
% %     spec_fit = f_hz(band).^P(1)*exp(1)^P(2);
% %     loglog(f_hz(band),spec_fit,[mark(i),'.-'],'HandleVisibility', 'off')
% %     loglog(f_hz(1:13),cell2mat(eeg.etc.FOOOF_results(i)).ap_fit(1:13),mark(i))
% end
% % legend('-DynamicLegend')
% 
% P(1)
%%%%%%%%%%
% bp = [];
% for i=1:eeg.nbchan
%     bp = [bp;cell2mat(eeg.etc.FOOOF_results(i)).bandpowers'];
% end    
% mean(bp,1)

%% Tests for filters,asr,ica... etc:
%
% asr = pop_select(EEG,'time',[pre_start_time 300]);
% 
% sel_chans = {'FP1','FP2','F7','F8'};
% chan_idx = cellfun(@(x) ismember(x,sel_chans),{asr.chanlocs.labels},'UniformOutput',0);
% disp('Selecting channels:');disp({asr.chanlocs.labels})
% asr = pop_select(asr,'channel',find(cell2mat(chan_idx)));
%
% asr.etc = rmfield(asr.etc,["clean_channel_mask","clean_sample_mask","clean_drifts_kernel"]);
% vis_artifacts(asr,filtd);
%
%%% Filtering and cleaning options:
%
% % EEGLAB built-in
% EEG = pop_eegfilt(EEG,0.5,80);%,'firtype','fir1');
%
% %Windowed FIR
% freqs = [0.5 80];
% wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,EEG.srate,df);
% EEG = pop_firws(EEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);
%
% %Zapline
% EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('noisefreqs',60, ...
%                                                  'chunkLength',0,'adaptiveNremove',true, ...
%                                                  'fixedNremove',1,'plotResults',0));
%
% %Artifact Subspace Reconstruction
% EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.1, ...
%     'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
%     'WindowCriterion',0.5,'BurstRejection','off','Distance','Euclidian', ...
%     'WindowCriterionTolerances',[-Inf 10]);
%
%%% ICA: try amica
%
%Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
% thresholds = [0.7 1;0 0.2;0 0.2;0 0.2;0 0.1;0 0.1;0 1];
% EEG = pop_runica(EEG,'icatype','runica','extended',1,'interrupt','off');
% EEG = iclabel(EEG);
% EEG = pop_icflag(EEG,thresholds);
% reconsEEG = pop_subcomp(EEG,EEG.reject.gcompreject);
% 
% sel_chans = {'FP1','FP2','F7','F8'};
% chan_idx = cellfun(@(x) ismember(x,sel_chans),{reconsEEG.chanlocs.labels},'UniformOutput',0);
% reconsEEG = pop_select(reconsEEG,'channel',find(cell2mat(chan_idx)));%,'time',[pre_start_time/1000 250]);
% disp('Selecting channels:');disp({reconsEEG.chanlocs.labels})
% 
% EEG = pop_select(EEG,'channel',find(cell2mat(chan_idx)));%,'time',[pre_start_time/1000 250]);
% vis_artifacts(reconsEEG,EEG)

%% Spectral analysis tests:
%
% sel_chans = {'FP1','FP2','F7','F8'};
% chan_idx = cellfun(@(x) ismember(x,sel_chans),{EEG.chanlocs.labels},'UniformOutput',0);
% pre = pop_select(EEG,'time',[pre_start_time 300],'channel',find(cell2mat(chan_idx)));
% % pre = pop_select(EEG,);
% 
% % figure; hold on
% % %EEG_in,type_proc,channel_num,tlimits,cycles
% [spectra_db,f_hz,~,~,~] = spectopo(pre.data,0,EEG.srate, ...
%                                        'freq',[8 10 12],'freqrange',[0.5 60], ...
%                                        'chanlocs',pre.chanlocs,'plot','off');%, ...
% %                                        'plot','on','electrodes','labels');
% % hold off 
% 
% % for i = 1:4
% % figure; hold on
% % [ersp itc powbase times frequencies] = pop_newtimef(EEG,1,i,[],[1 0.5], ...
% %                                     'plottype','image','plotersp','on','plotitc','off');
% % tftopo(ersp,times,frequencies)
% % hold off
% % figure; hold on
% % [ersp itc powbase times frequencies] = pop_newtimef(reconsEEG,1,i,[],[1 0.5], ...
% %                                     'plottype','image','plotersp','on','plotitc','off');
% % % hold off
% % end
% 
% % Define band frequency ranges:
% delta_hz = find(f_hz>1 & f_hz<4);
% theta_hz = find(f_hz>4 & f_hz<7);
% alpha_hz = find(f_hz>8 & f_hz<13);
% lbeta_hz = find(f_hz>13 & f_hz<20);
% hbeta_hz = find(f_hz>20 & f_hz<30);
% lgamma_hz = find(f_hz>30 & f_hz<45);
% 
% % Compute band power
% pre_delta_power = 10^(mean(spectra_db(delta_hz))/10);
% pre_theta_power = 10^(mean(spectra_db(theta_hz))/10);
% pre_alpha_power = 10^(mean(spectra_db(alpha_hz))/10);
% pre_lbeta_power = 10^(mean(spectra_db(lbeta_hz))/10);
% pre_hbeta_power = 10^(mean(spectra_db(hbeta_hz))/10);
% pre_lgamma_power=10^(mean(spectra_db(lgamma_hz))/10);
% 
% T=table(subid,pre_delta_power,pre_theta_power,pre_alpha_power,pre_lbeta_power,pre_hbeta_power,pre_lgamma_power)

%% Find pre and post stress, short and long timestamps:

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

pre_game_EEG = pop_select(EEG,'time',[pre_game_start_time pre_game_end_time]);
pre_game_events = {pre_game_EEG.event.type};

pre_game_short_end = find(strcmp(pre_game_events,'TRIGGER EVENT T'));
pre_game_short_end = pre_game_short_end(end);
pre_game_long_end = find(strcmp(pre_game_events,'TRIGGER EVENT U'));
pre_game_long_end = pre_game_long_end(end);

pre_game_short_end_time = pre_game_EEG.event(pre_game_short_end).latency/EEG.srate;
pre_game_long_end_time = pre_game_EEG.event(pre_game_long_end).latency/EEG.srate;

if pre_game_short_end<pre_game_long_end
    pre_short_timestamp = [pre_game_start_time...
                           pre_game_start_time+pre_game_short_end_time];
    pre_long_timestamp = [pre_game_start_time+pre_game_short_end_time...
                          pre_game_start_time+pre_game_long_end_time];
else
    pre_long_timestamp = [pre_game_start_time...
                          pre_game_start_time+pre_game_long_end_time];
    pre_short_timestamp = [pre_game_start_time+pre_game_long_end_time...
                           pre_game_start_time+pre_game_short_end_time];
end

post_game_start = find(strcmp(all_events,'TRIGGER EVENT I'));
post_game_start_time = EEG.event(post_game_start).latency/EEG.srate;

post_game_end = find(strcmp(all_events(post_game_start:end),block_end_marker));
post_game_end = post_game_end(1)+post_game_start;
post_game_end_time = EEG.event(post_game_end).latency/EEG.srate;

post_game_EEG = pop_select(EEG,'time',[post_game_start_time post_game_end_time]);
post_game_events = {post_game_EEG.event.type};

post_game_short_end = find(strcmp(post_game_events,'TRIGGER EVENT T'));
post_game_short_end = post_game_short_end(end);
post_game_long_end = find(strcmp(post_game_events,'TRIGGER EVENT U'));
post_game_long_end = post_game_long_end(end);

post_game_short_end_time = post_game_EEG.event(post_game_short_end).latency/EEG.srate;
post_game_long_end_time = post_game_EEG.event(post_game_long_end).latency/EEG.srate;

if post_game_short_end<post_game_long_end
    post_short_timestamp = [post_game_start_time...
                            post_game_start_time+post_game_short_end_time];
    post_long_timestamp = [post_game_start_time+post_game_short_end_time...
                           post_game_start_time+post_game_long_end_time];
else
    post_long_timestamp = [post_game_start_time...
                           post_game_start_time+post_game_long_end_time];
    post_short_timestamp = [post_game_start_time+post_game_long_end_time...
                            post_game_start_time+post_game_short_end_time];
end

%% Test E/I Variability plots:
% band = [1 100];chan_idx=13;
% 
% pre_EEG = pop_select(EEG,'time',[0 300]);t_epoch=60;
% pre_EEG = eeg_regepochs(pre_EEG,'recurrence',t_epoch,'limits',[0 t_epoch],'rmbase',NaN);
% % pre_EEG = pop_rejepoch(pre_EEG,[1 0 0 0],0);
% [pre_ersp,pre_exponents,pre_offsets,pre_times] = ap_exp_ts(pre_EEG,chan_idx,band);
% figure;hold on;
% hist(pre_exponents)
% % plot(pre_timepoints,pre_exponents);
% % plot(timepoints,offsets);
% % hold off;
% 
% pre_short_EEG = pop_select(EEG,'time',pre_short_timestamp);
% % sprintf('Trials: %d',sum(ismember({pre_short_EEG.event.type},{'TRIGGER EVENT N'})))
% s_idx = find(ismember({pre_short_EEG.event.type},{'TRIGGER EVENT S'}));
% % sprintf('Num patches: %d',size(s_idx,2))
% % n_patch_idx = [1,s_idx];
% % for idx=2:size(n_patch_idx,2)
% %     sprintf('Trials in patch %d: %d',idx-1,sum(ismember({pre_short_EEG.event(n_patch_idx(idx-1):n_patch_idx(idx)).type},{'TRIGGER EVENT N'})))
% % end
% s_timestamp = [pre_short_EEG.event(s_idx).latency]/pre_short_EEG.srate;
% short_patch_durations = [s_timestamp(1),diff(s_timestamp)];
% pre_short_epoch = pop_epoch(pre_short_EEG,{'TRIGGER EVENT S'},[-min(short_patch_durations)+0.5 0]);
% sprintf('Num patches (after epoching): %d',pre_short_epoch.trials)
% [pre_short_ersp,pre_short_exponents,pre_short_offsets,pre_short_timepoints] = ap_exp_ts(pre_short_epoch,chan_idx,band);
% % short_trial_times = nan(pre_short_epoch.trials,12);
% % for idx = 1:pre_short_epoch.trials
% %     T_idx = find(ismember(pre_short_epoch.epoch(idx).eventtype,{'TRIGGER EVENT T'}));
% %     short_trial_times(idx,1:size(T_idx,2)) = cell2mat(pre_short_epoch.epoch(idx).eventlatency(T_idx))/1000;
% % end
% % mean_short_trial_times = mean(short_trial_times,1,'omitnan');
% % 
% % figure;hold on;
% % hist(short_exponents,10)
% % plot(short_exponents);
% % xline(short_trial_times(:,3),'--');
% % xline(mean_short_trial_times(3));
% % % plot(timepoints,offsets);
% % hold off;
% 
% pre_long_EEG = pop_select(EEG,'time',pre_long_timestamp);
% % sprintf('Trials: %d',sum(ismember({pre_long_EEG.event.type},{'TRIGGER EVENT N'})))
% s_idx = find(ismember({pre_long_EEG.event.type},{'TRIGGER EVENT S'}));
% % sprintf('Num patches: %d',size(s_idx,2))
% % n_patch_idx = [1,s_idx];
% % for idx=2:size(n_patch_idx,2)
% %     sprintf('Trials in patch %d: %d',idx-1,sum(ismember({pre_long_EEG.event(n_patch_idx(idx-1):n_patch_idx(idx)).type},{'TRIGGER EVENT N'})))
% % end
% s_timestamp = [pre_long_EEG.event(s_idx).latency]/pre_long_EEG.srate;
% long_patch_durations = [s_timestamp(1),diff(s_timestamp)];
% pre_long_epoch = pop_epoch(pre_long_EEG,{'TRIGGER EVENT S'},[-min(long_patch_durations)+0.5 0]);
% % sprintf('Num patches (after epoching): %d',pre_long_epoch.trials)
% [pre_long_ersp,pre_long_exponents,pre_long_offsets,pre_long_timepoints] = ap_exp_ts(pre_long_epoch,chan_idx,band);
% % long_trial_times = nan(pre_long_epoch.trials,12);
% % for idx = 1:pre_long_epoch.trials
% %     U_idx = find(ismember(pre_long_epoch.epoch(idx).eventtype,{'TRIGGER EVENT U'}));
% %     long_trial_times(idx,1:size(U_idx,2)) = cell2mat(pre_long_epoch.epoch(idx).eventlatency(U_idx))/1000;
% % end
% % mean_long_trial_times = mean(long_trial_times,1,'omitnan');
% % 
% % figure;hold on;
% % hist(long_exponents,10)
% % plot(long_exponents);
% % xline(long_trial_times(:,3),'--');
% % xline(mean_long_trial_times(3));
% % % plot(timepoints,offsets);
%%
bands = [[1 100];[1 10];[10 100]];
sel_chans = {'FCz','FZ','CZ'};
chan_id_fn = cellfun(@(x) ismember(x,sel_chans),{EEG.chanlocs.labels},'UniformOutput',0); %%Completely unnecessary
chan_idxs = find(cell2mat(chan_id_fn)); chan_idx = chan_idxs(2);
fprintf([repmat('=',1,80),'\n',sel_chans(chan_idx),'\n',repmat('=',1,80),'\n'])

%%
% for chan_idx = chan_idxs
%     for idx = 1%:3
band = bands(1,:);
fprintf('Processing pre task resting state...\n')
pre_EEG = pop_select(EEG,'time',pre_timestamp);t_epoch=60;
pre_EEG = eeg_regepochs(pre_EEG,'recurrence',t_epoch,'limits',[0 t_epoch],'rmbase',NaN);
% pre_EEG = pop_rejepoch(pre_EEG,[1 0 0 0],0);
[pre_times,pre_freqs,pre_ersp,pre_exponents,pre_offsets,...
 pre_bp,pre_ap_bp] = ap_exp_ts(pre_EEG,chan_idx,band);

fprintf('Processing pre stress short TT resting state...\n')
pre_short_EEG = pop_select(EEG,'time',pre_short_timestamp);
s_idx = find(ismember({pre_short_EEG.event.type},{'TRIGGER EVENT S'}));
s_timestamp = [pre_short_EEG.event(s_idx).latency]/pre_short_EEG.srate;
short_patch_durations = [s_timestamp(1),diff(s_timestamp)];
pre_short_epoch = pop_epoch(pre_short_EEG,{'TRIGGER EVENT S'},[-min(short_patch_durations)+0.5 0]);
[pre_short_times,pre_short_freqs,pre_short_ersp,pre_short_exponents,pre_short_offsets,...
 pre_short_bp,pre_short_ap_bp] = ap_exp_ts(pre_short_epoch,chan_idx,band);

fprintf('Processing pre stress long TT resting state...\n')
pre_long_EEG = pop_select(EEG,'time',pre_long_timestamp);
s_idx = find(ismember({pre_long_EEG.event.type},{'TRIGGER EVENT S'}));
s_timestamp = [pre_long_EEG.event(s_idx).latency]/pre_long_EEG.srate;
long_patch_durations = [s_timestamp(1),diff(s_timestamp)];
pre_long_epoch = pop_epoch(pre_long_EEG,{'TRIGGER EVENT S'},[-min(long_patch_durations)+0.5 0]);
[pre_long_times,pre_long_freqs,pre_long_ersp,pre_long_exponents,pre_long_offsets,...
 pre_long_bp,pre_long_ap_bp] = ap_exp_ts(pre_long_epoch,chan_idx,band);

fprintf('Processing post stress short TT resting state...\n')
post_short_EEG = pop_select(EEG,'time',post_short_timestamp);
s_idx = find(ismember({post_short_EEG.event.type},{'TRIGGER EVENT S'}));
s_timestamp = [post_short_EEG.event(s_idx).latency]/post_short_EEG.srate;
post_short_patch_durations = [s_timestamp(1),diff(s_timestamp)];
post_short_epoch = pop_epoch(post_short_EEG,{'TRIGGER EVENT S'},[-min(post_short_patch_durations)+0.5 0]);
[post_short_times,post_short_freqs,post_short_ersp,post_short_exponents,post_short_offsets,...
 post_short_bp,post_short_ap_bp] = ap_exp_ts(post_short_epoch,chan_idx,band);

fprintf('Processing post stress long TT resting state...\n')
post_long_EEG = pop_select(EEG,'time',post_long_timestamp);
s_idx = find(ismember({post_long_EEG.event.type},{'TRIGGER EVENT S'}));
s_timestamp = [post_long_EEG.event(s_idx).latency]/post_long_EEG.srate;
post_long_patch_durations = [s_timestamp(1),diff(s_timestamp)];
post_long_epoch = pop_epoch(post_long_EEG,{'TRIGGER EVENT S'},[-min(post_long_patch_durations)+0.5 0]);
[post_long_times,post_long_freqs,post_long_ersp,post_long_exponents,post_long_offsets,...
 post_long_bp,post_long_ap_bp] = ap_exp_ts(post_long_epoch,chan_idx,band);

%%
state_conditions = {'pre','pre_short','pre_long','post_short','post_long'};
spec_vars = {'times','freqs','ersp','exponents','offsets'};
freq_bands = {'delta','theta','alpha','lbeta','hbeta','lgamma','hgamma'};
% sc_idx = find(ismember(state_conditions,'pre'));
% freq_band_idx = find(ismember(freq_bands,'lgamma'));

spec_data = struct();
for sc_idx = 1:size(state_conditions,2)
    state_condition = state_conditions(sc_idx);
    for sv_idx = 1:size(spec_vars,2)
        sv = spec_vars(sv_idx);
        spec_data.(char(state_condition)).(char(sv)) = ...
            eval(sprintf('%s_%s',char(state_condition),char(sv)));
    end
    for fb_idx = 1:size(freq_bands,2)
        freq_band = freq_bands(fb_idx);
        spec_data.(char(state_condition)).([char(freq_band),'_bp']) = ...
                eval(sprintf('%s_bp(%d,:)',char(state_condition),fb_idx));
        spec_data.(char(state_condition)).([char(freq_band),'_ap_bp']) = ...
                eval(sprintf('%s_ap_bp(%d,:)',char(state_condition),fb_idx));
        
%         [pre_exp_pdf,pre_exp] = ksdensity(pre_exponents);
%         [pre_short_exp_pdf,pre_short_exp] = ksdensity(pre_short_exponents);
%         [pre_long_exp_pdf,pre_long_exp] = ksdensity(pre_long_exponents);
%         [post_short_exp_pdf,post_short_exp] = ksdensity(post_short_exponents);
%         [post_long_exp_pdf,post_long_exp] = ksdensity(post_long_exponents);
%
%         pre_stat = sprintf('mean: %1.2f\nstd: %1.2f',mean(pre_exponents,'omitnan'),std(pre_exponents,'omitnan'));
%         pre_short_stat = sprintf('mean: %1.2f\nstd: %1.2f',mean(pre_short_exponents,'omitnan'),std(pre_short_exponents,'omitnan'));
%         pre_long_stat = sprintf('mean: %1.2f\nstd: %1.2f',mean(pre_long_exponents,'omitnan'),std(pre_long_exponents,'omitnan'));
%         post_short_stat = sprintf('mean: %1.2f\nstd: %1.2f',mean(post_short_exponents,'omitnan'),std(post_short_exponents,'omitnan'));
%         post_long_stat = sprintf('mean: %1.2f\nstd: %1.2f',mean(post_long_exponents,'omitnan'),std(post_long_exponents,'omitnan'));
        
%         [pre_bp_pdf,pre_bp] = ksdensity(pre_bandpowers(freq_band_idx,:));
%         [pre_short_bp_pdf,pre_short_bp] = ksdensity(pre_short_bandpowers(freq_band_idx,:));
%         [pre_long_bp_pdf,pre_long_bp] = ksdensity(pre_long_bandpowers(freq_band_idx,:));
%         [post_short_bp_pdf,post_short_bp] = ksdensity(post_short_bandpowers(freq_band_idx,:));
%         [post_long_bp_pdf,post_long_bp] = ksdensity(post_long_bandpowers(freq_band_idx,:));

%         figure;hold on;ax=gca;c=ax.ColorOrder;close;
%         figure;hold on;
%         bin_width = 1;

%         histogram(pre_exponents,'BinWidth',0.15,'Normalization','pdf',...
%                   'FaceColor',c(1,:),'FaceAlpha',0.7,'HandleVisibility','off','Visible','off')
%         histogram(pre_short_exponents,'BinWidth',0.15,'Normalization','pdf',...
%                   'FaceColor',c(2,:),'FaceAlpha',0.6,'HandleVisibility','off','Visible','on')
%         histogram(pre_long_exponents,'BinWidth',0.15,'Normalization','pdf',...
%                   'FaceColor',c(3,:),'FaceAlpha',0.5,'HandleVisibility','off','Visible','on')
%         histogram(post_short_exponents,'BinWidth',0.15,'Normalization','pdf',...
%                   'FaceColor',c(4,:),'FaceAlpha',0.4,'HandleVisibility','off','Visible','on')
%         histogram(post_long_exponents,'BinWidth',0.15,'Normalization','pdf',...
%                   'FaceColor',c(5,:),'FaceAlpha',0.3,'HandleVisibility','off','Visible','on')

%         histogram(pre_bandpowers(freq_band_idx,:),'BinWidth',bin_width,'Normalization','pdf',...
%                   'FaceColor',c(1,:),'FaceAlpha',0.7,...
%                   'DisplayName','Pre Game','HandleVisibility','on','Visible','off')
% 
%         histogram(pre_short_bandpowers(freq_band_idx,:),'BinWidth',bin_width,'Normalization','pdf',...
%                   'FaceColor',c(2,:),'FaceAlpha',0.6,...
%                   'DisplayName','Pre Stress Short TT','HandleVisibility','on','Visible','on')
% 
%         histogram(pre_long_bandpowers(freq_band_idx,:),'BinWidth',bin_width,'Normalization','pdf',...
%                   'FaceColor',c(3,:),'FaceAlpha',0.5,...
%                   'DisplayName','Pre Stress Long TT','HandleVisibility','on','Visible','off')
% 
%         histogram(post_short_bandpowers(freq_band_idx,:),'BinWidth',bin_width,'Normalization','pdf',...
%                   'FaceColor',c(4,:),'FaceAlpha',0.4,...
%                   'DisplayName','Post Stress Short TT','HandleVisibility','on','Visible','on')
% 
%         histogram(post_long_bandpowers(freq_band_idx,:),'BinWidth',bin_width,'Normalization','pdf',...
%                   'FaceColor',c(5,:),'FaceAlpha',0.3,...
%                   'DisplayName','Post Stress Long TT','HandleVisibility','on','Visible','off')

%         area(pre_bp_pdf,pre_bp,'EdgeColor',c(1,:),'LineStyle','-',...
%             'FaceColor',c(1,:),'FaceAlpha',0.7,'DisplayName','Pre Game',...
%             'Visible','on');
%         
%         area(pre_short_bp_pdf,pre_short_bp,'EdgeColor',c(2,:),'LineStyle','-',...
%             'FaceColor',c(2,:),'FaceAlpha',0.6,'DisplayName','Pre Stress Short TT',...
%             'Visible','on')
%         
%         area(pre_long_bp_pdf,pre_long_bp,'EdgeColor',c(3,:),'LineStyle','-',...
%             'FaceColor',c(3,:),'FaceAlpha',0.5,'DisplayName','Pre Stress Long TT',...
%             'Visible','on')
%         
%         area(post_short_bp_pdf,post_short_bp,'EdgeColor',c(4,:),'LineStyle','-',...
%             'FaceColor',c(4,:),'FaceAlpha',0.4,'DisplayName','Post Stress Short TT',...
%             'Visible','on')
%         
%         area(post_long_bp_pdf,post_long_bp,'EdgeColor',c(5,:),'LineStyle','-',...
%             'FaceColor',c(5,:),'FaceAlpha',0.3,'DisplayName','Post Stress Long TT',...
%             'Visible','on')
        
%         legend();
%         t_txt_exp = sprintf('Sub: %d \t Channel: %s \t Band: [%s]',...
%                         subid,EEG.chanlocs(chan_idx).labels,join(string(band),' '));
%         title(t_txt_exp)
%         t_txt_bp = sprintf('Sub: %d \t Channel: %s \t Band: %s',...
%                         subid,EEG.chanlocs(chan_idx).labels,char(freq_bands(freq_band_idx)));
%         title(t_txt_bp)
%         hold off;
        
%         exp_fname = sprintf('%d_%s_%s_%s',subid,EEG.chanlocs(chan_idx).labels,join(string(band),'-'),char(state_condition));
%         bp_fname = sprintf('%d_%s_%s_bp',subid,EEG.chanlocs(chan_idx).labels,char(state_condition));
%         saveas(gcf,['./fig/' fname,'.svg'])
    end
end
%     end
% %     close all;
% end
% % end
%% Save processed data:
%%Comment out for single sub
% % pop_saveset(EEG,'filename',[num2str(subid) '_preproc'],'filepath',[pwd '/data/clean_set/'],'savemode','onefile');
% global delta_hz; global theta_hz; global alpha_hz; global lbeta_hz; global hbeta_hz; global lgamma_hz
% delta_hz = find(f_hz>1 & f_hz<4); theta_hz = find(f_hz>4 & f_hz<7); alpha_hz = find(f_hz>8 & f_hz<13);
% lbeta_hz = find(f_hz>13 & f_hz<20); hbeta_hz = find(f_hz>20 & f_hz<30); lgamma_hz = find(f_hz>30 & f_hz<45);
% subT=frg_extract_bp(subid,EEG,pre_start_time,short_timestamp,long_timestamp);
% bp_table = [bp_table;subT];
save(sprintf('./data/spec_data/spec_%d_%s',subid,EEG.chanlocs(chan_idx).labels),'spec_data')
end
% writetable(bp_table,'bandpowers.csv')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%==Function Definitions==%%%%%%%%%%%%%%%%%%%%%%%%%

function [times,frequencies,ersp,exponents,offsets,bp,ap_bp] = ap_exp_ts(eeg,chan_idx,band)
    [ersp,~,~,times,frequencies,~,~,~] = ...
    pop_newtimef(eeg,1,chan_idx,[eeg.xmin eeg.xmax]*1000,[1 0.6],...
                'freqs',[1 120],'baseline',NaN,'basenorm','off',...
                'plotersp','off','plotitc','off','scale','abs'); %'scale','log'

    ersp_abs = ersp;%10.^(ersp./10);
    lo_f = band(1); hi_f = band(2);
    half_win = 0.75; %%Sliding window!

    exponents = []; offsets = [];
    bp = []; ap_bp = [];
%     figure;ax=gca;ax.XScale='log';ax.YScale='log';hold on;
    for tpoint=times./1000
        x = frequencies(frequencies>lo_f & frequencies<hi_f);
        y = mean(ersp_abs((frequencies>lo_f & frequencies<hi_f),...
                        times>(tpoint-half_win)*1000 & times<(tpoint+half_win)*1000),2);
%         loglog(x,y)
        fooof_msg = fprintf('FOOOFing... Timepoint:%.4f\n',tpoint);
        fooof_res = fooof(x,y,[x(1),x(end)],struct(),0);
        offsets(end+1) = fooof_res.aperiodic_params(1);
        exponents(end+1) = fooof_res.aperiodic_params(2);
        bp(:,end+1) = fooof_res.bandpowers;
        ap_bp(:,end+1) = fooof_res.ap_bandpowers;
        fprintf(repmat('\b', 1, (fooof_msg)));
        
%         x = log(frequencies(frequencies>lo_f & frequencies<hi_f));
%         y = log(mean(ersp_abs(find(frequencies>lo_f & frequencies<hi_f),...
%                         find(times>(tpoint-0.75)*1000 & times<(tpoint+0.75)*1000)),2));
%         [P,S] = polyfit(x,y,1);
%         R_squared = 1 - (S.normr/norm(y - mean(y)))^2
%         exponents(end+1)=P(1);offsets(end+1)=P(2);
%         spec_fit = frequencies(frequencies>lo_f & frequencies<hi_f).^P(1)*exp(1)^P(2);
    end
    sprintf('Done!')
end

% function [seg_delta_power,seg_theta_power,seg_alpha_power,seg_lbeta_power,...
%     seg_hbeta_power,seg_lgamma_power] = seg_chan_bp_abs(eeg,time_stamps,sel_chan)
%     eeg = pop_select(eeg,'time',time_stamps,'channel',sel_chan);
%     [spectra_db,f_hz,~,~,~] = spectopo(eeg.data,0,eeg.srate, ...
%                                            'freq',[8 10 12],'freqrange',[0.5 60], ...
%                                            'chanlocs',eeg.chanlocs,'plot','off');
%     delta_hz = f_hz(f_hz>1 & f_hz<4);theta_hz = f_hz(f_hz>4 & f_hz<7);alpha_hz = f_hz(f_hz>8 & f_hz<13);
%     lbeta_hz = f_hz(f_hz>13 & f_hz<20);hbeta_hz = f_hz(f_hz>20 & f_hz<30);lgamma_hz = f_hz(f_hz>30 & f_hz<45);
%     seg_delta_power = mean(10.^(spectra_db(:,delta_hz)/10),'all');
%     seg_theta_power = mean(10.^(spectra_db(:,theta_hz)/10),'all');
%     seg_alpha_power = mean(10.^(spectra_db(:,alpha_hz)/10),'all');
%     seg_lbeta_power = mean(10.^(spectra_db(:,lbeta_hz)/10),'all');
%     seg_hbeta_power = mean(10.^(spectra_db(:,hbeta_hz)/10),'all');
%     seg_lgamma_power=mean(10.^(spectra_db(:,lgamma_hz)/10),'all');
% end
%
% function [eeg,specdata] = seg_chan_bp_specpar(eeg,time_stamps,sel_chan)
%     chan_id = cellfun(@(x) ismember(x,sel_chan),{eeg.chanlocs.labels},'UniformOutput',0);chan_idxs = find(cell2mat(chan_id));
%     eeg = pop_select(eeg,'time',time_stamps,'channel',sel_chan);
%     eeg = eeg_fooof(eeg,'channel',1:eeg.nbchan,[1 4]);%,[eeg.xmin*1000,eeg.xmax*1000],100,struct('aperiodic_mode','knee'));
%     [spec_db,f_hz,~,~,~] = spectopo(eeg.data,0,eeg.srate, ...
%                                        'freq',[10],'freqrange',[1 4], ...
%                                        'winsize',1000,'overlap',250,...
%                                        'chanlocs',eeg.chanlocs,'plot','off');
% 
%     specdata = arrayfun(@(y) 10^(y/10), spec_db);%specfreqs = f_hz';
% %     specpar_out = fooof(specfreqs,specdata(1,:),[1 60],struct(),0);
% end
% 
% function T = frg_extract_bp(subid,eeg,pre_start_time,short_timestamps,long_timestamps)
%     sel_chans = {'FP1','FP2','F7','F8'};
%     chan_idx = cellfun(@(x) ismember(x,sel_chans),{eeg.chanlocs.labels},'UniformOutput',0);
%     
%     pre = pop_select(eeg,'time',[pre_start_time 300],'channel',find(cell2mat(chan_idx)));
%     
%     [spectra_db,f_hz,~,~,~] = spectopo(pre.data,0,eeg.srate, ...
%                                            'freq',[8 10 12],'freqrange',[0.5 60], ...
%                                            'chanlocs',pre.chanlocs,'plot','off');
%     
% %     delta_hz = find(f_hz>1 & f_hz<4);
% %     theta_hz = find(f_hz>4 & f_hz<7);
% %     alpha_hz = find(f_hz>8 & f_hz<13);
% %     lbeta_hz = find(f_hz>13 & f_hz<20);
% %     hbeta_hz = find(f_hz>20 & f_hz<30);
% %     lgamma_hz = find(f_hz>30 & f_hz<45);
%     
%     pre_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
%     pre_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
%     pre_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
%     pre_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
%     pre_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
%     pre_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);
% 
%     short = pop_select(eeg,'time',short_timestamps,'channel',find(cell2mat(chan_idx)));
% 
%     [spectra_db,f_hz,~,~,~] = spectopo(short.data,0,eeg.srate, ...
%                                        'freq',[8 10 12],'freqrange',[0.5 60], ...
%                                        'chanlocs',pre.chanlocs,'plot','off');
%     
%     short_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
%     short_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
%     short_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
%     short_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
%     short_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
%     short_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);
%     
%     long = pop_select(eeg,'time',long_timestamps,'channel',find(cell2mat(chan_idx)));
% 
%     [spectra_db,f_hz,~,~,~] = spectopo(long.data,0,eeg.srate, ...
%                                        'freq',[8 10 12],'freqrange',[0.5 60], ...
%                                        'chanlocs',pre.chanlocs,'plot','off');
%     
%     long_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
%     long_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
%     long_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
%     long_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
%     long_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
%     long_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);
% 
%     T=table(subid,pre_delta_power,pre_theta_power,pre_alpha_power, ...
%                   pre_lbeta_power,pre_hbeta_power,pre_lgamma_power,...
%                   short_delta_power,short_theta_power,short_alpha_power,...
%                   short_lbeta_power,short_hbeta_power,short_lgamma_power,...
%                   long_delta_power,long_theta_power,long_alpha_power,...
%                   long_lbeta_power,long_hbeta_power,long_lgamma_power);
% end