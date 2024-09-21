%% Load and reref EEG:
clear;clc;
current_dir = pwd;
% cd "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1";
% cd "C:\Users\Peeusa\OneDrive\Documents\MATLAB\eeglab2023.1";
% cd '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
eeglab;close;
cd(current_dir)
% origchanlocs = readlocs('./live_amp_32.ced');

%%
data_path = 'F:\phd works\experimental data\analysis\final data analysis\eeg_analysis\data\Loop_data\';
fname = 'p2_import.set';
% set_files = dir('*.set');
% fnames = {set_files(:).name};
% 
% for idx = 1:size(fnames,2)
% fname = cell2mat(fnames(idx));
% p_id = fname(1:find(ismember(fname,'_'))-1);

%%
EEG = pop_loadset(fname);
%%
% %add fcz and reref by ch. nos. [29 30]
% EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
% EEG = pop_chanedit(EEG,'load','./live_amp_32.ced');
% EEG = pop_reref(EEG,[29 30]);
% EEG = pop_interp(EEG,origchanlocs);
% EEG = pop_reref(EEG,[]);

%% Bandpass, line noise cleaning and ASR:
%zapline,asr etc.
freqs = [0.5 120];
wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,EEG.srate,df);
EEG = pop_firws(EEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);

EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('noisefreqs',60, ...
                                   'chunkLength',0,'adaptiveNremove',true, ...
                                   'fixedNremove',1,'plotResults',0));

% filtd = pop_select(EEG,'time',[pre_start_time 300]);
% 
% sel_chans = {'FP1','FP2','F7','F8'};
% chan_idx = cellfun(@(x) ismember(x,sel_chans),{filtd.chanlocs.labels},'UniformOutput',0);
% disp('Selecting channels - filtered:');disp({filtd.chanlocs.labels})
% filtd = pop_select(filtd,'channel',find(cell2mat(chan_idx)));

% EEG_asr = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.1, ...
%     'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
%     'WindowCriterion',0.5,'BurstRejection','off','Distance','Euclidian', ...
%     'WindowCriterionTolerances',[-Inf 10]);
% 
% EEG_asr = pop_interp(EEG_asr,origchanlocs);

%% Segment short and long envs:
all_events = {EEG.event.type};
leave_idxs = find(ismember(all_events,'N'));
travel_times=zeros(size(leave_idxs));
for idx = [1:size(leave_idxs,2)]
    travel_times(idx) = EEG.event(leave_idxs(idx)+1).latency/EEG.srate - EEG.event(leave_idxs(idx)).latency/EEG.srate;
end

short_env_leave_idxs = leave_idxs(travel_times<10);
long_env_leave_idxs = leave_idxs(travel_times>10);

short_env_start_time = EEG.event(1).latency/EEG.srate;
short_env_end_time = EEG.event(short_env_leave_idxs(end)).latency/EEG.srate;
short_env_timestamps = [short_env_start_time short_env_end_time];

% Since short env always precedes long env...
long_env_start_time = EEG.event(short_env_leave_idxs(end)+1).latency/EEG.srate;
long_env_end_time = EEG.event(leave_idxs(end)).latency/EEG.srate;
long_env_timestamps = [long_env_start_time long_env_end_time];

EEG_short = pop_select(EEG,'time',short_env_timestamps);
EEG_long = pop_select(EEG,'time',long_env_timestamps);
%%
