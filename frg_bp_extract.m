%% Notes:
% 1. Use clean_rawdata asr for bad channel rejection right at the top and
% multiple times
% 2. Avg ref before interp and ica

%% Load EEGLAB: DO NOT add to MATLAB path - causes problems with BioSig
clear;clc;
current_dir = pwd;
% cd "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1";
cd '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
eeglab;close;
cd(current_dir)

%% Path to data:
data_path = strcat(current_dir,'/data/');
subid = 7873;

%%Comment out for single sub: Exclude badids
% badids = [37532 38058 39862 42125 43543 45194 45528 46037 47678 47744 47801 48238 48278];
% subids = readmatrix(fullfile(data_path,'subids.txt')); subids = setdiff(subids,badids);
% bp_table = table();
%% Main loop begins:

%Comment out for single sub: Set up loop vars
% for idx = 1:length(subids)
% clearvars -except idx subids data_path bp_table
% global subid; subid = subids(idx); 

repelem('=',80)
disp(subid)
% sprintf('%d of %d',idx,length(subids))
repelem('=',80)

%%% Load data and channel locations:
data_file = dir(fullfile(data_path,string(subid),'*.edf'));
data_file_path = fullfile(data_file(end).folder,data_file(end).name);
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

% all_events = {EEG.event.type};
% if size(find(ismember(all_events,'TRIGGER EVENT A')),2)>0
%     pre_start = find(ismember(all_events,'TRIGGER EVENT A'));
%     pre_start_time = EEG.event(pre_start(1)).latency/EEG.srate;
%     sprintf('Start time: %0.2g seconds',pre_start_time)
% else
%     pre_start_time=0;
% end
pre_start_time=0;
% EEG = pop_select(EEG,'time',[pre_start_time 300]);

%Online reference was FCz, which is not there in the data. Adding an
%empty channel at index 19 (same as chanloc file) before loading the channel locations.
EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
EEG = pop_chanedit(EEG,'load','./data/Statnet_F3F4FCz.ced');

%Re-referencing to mastoids A1 A2, interpolating FCz (17 in data, 19 in chanloc), and
%re-referencing to averef
EEG = pop_reref(EEG,{'A1','A2'});
EEG = pop_interp(EEG,origchanlocs);
EEG = pop_reref(EEG,[]);
% raw = EEG;%pop_select(EEG,'time',[pre_start_time 300]);

freqs = [0.5 80];
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

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.1, ...
    'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
    'WindowCriterion',0.5,'BurstRejection','off','Distance','Euclidian', ...
    'WindowCriterionTolerances',[-Inf 10]);

EEG = pop_interp(EEG,origchanlocs);

sel_chans = {'FP1','FP2','F7','F8'};
%%
% [seg_delta_power,seg_theta_power,seg_alpha_power,seg_lbeta_power,...
%     seg_hbeta_power,seg_lgamma_power] = seg_chan_bp_abs(EEG,[0 300],sel_chans)
eeg = seg_chan_bp_specpar(EEG,[0 300],sel_chans);
bp = [];
for i=1:eeg.nbchan
    bp = [bp;cell2mat(eeg.etc.FOOOF_results(i)).bandpowers']
end    

[seg_delta_power,seg_theta_power,seg_alpha_power,seg_lbeta_power,...
seg_hbeta_power,seg_lgamma_power] = seg_chan_bp_abs(EEG,[0 300],sel_chans)
%% Check filter,asr,ica... etc:
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
%
%% Spectral analysis:
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
%
%% Trim short and long timestamps:
%
% all_events = {EEG.event.type};
% game_start = find(strcmp(all_events,'TRIGGER EVENT B'));
% game_end = find(strcmp(all_events(game_start:end),'TRIGGER EVENT Z'));
% game_end = game_end(1)+game_start;
% game_start_time = EEG.event(game_start).latency/EEG.srate;
% game_end_time = EEG.event(game_end).latency/EEG.srate;
% 
% game_EEG = pop_select(EEG,'time',[game_start_time game_end_time]);
% game_events = {game_EEG.event.type};
% 
% game_short_end = find(strcmp(game_events,'TRIGGER EVENT T'));game_short_end = game_short_end(end);
% game_long_end = find(strcmp(game_events,'TRIGGER EVENT U'));game_long_end = game_long_end(end);
% 
% game_short_end_time = game_EEG.event(game_short_end).latency/EEG.srate;
% game_long_end_time = game_EEG.event(game_long_end).latency/EEG.srate;
% 
% if game_short_end<game_long_end
%     short_timestamp = [game_start_time game_start_time+game_short_end_time];
%     long_timestamp = [game_start_time+game_short_end_time game_start_time+game_long_end_time];
% else
%     long_timestamp = [game_start_time game_start_time+game_long_end_time];
%     short_timestamp = [game_start_time+game_long_end_time game_start_time+game_short_end_time];
% end
%% Save preprocessed data:
%%Comment out for single sub
% pop_saveset(EEG,'filename',[num2str(subid) '_preproc'],'filepath',[pwd '/data/clean_set/'],'savemode','onefile');
% global delta_hz; global theta_hz; global alpha_hz; global lbeta_hz; global hbeta_hz; global lgamma_hz
% delta_hz = find(f_hz>1 & f_hz<4); theta_hz = find(f_hz>4 & f_hz<7); alpha_hz = find(f_hz>8 & f_hz<13);
% lbeta_hz = find(f_hz>13 & f_hz<20); hbeta_hz = find(f_hz>20 & f_hz<30); lgamma_hz = find(f_hz>30 & f_hz<45);
% subT=frg_extract_bp(subid,EEG,pre_start_time,short_timestamp,long_timestamp);
% bp_table = [bp_table;subT];
% end
% writetable(bp_table,'bandpowers.csv')
%%
function [seg_delta_power,seg_theta_power,seg_alpha_power,seg_lbeta_power,...
    seg_hbeta_power,seg_lgamma_power] = seg_chan_bp_abs(eeg,time_stamps,sel_chan)
    eeg = pop_select(eeg,'time',time_stamps,'channel',sel_chan);
    [spectra_db,f_hz,~,~,~] = spectopo(eeg.data,0,eeg.srate, ...
                                           'freq',[8 10 12],'freqrange',[0.5 60], ...
                                           'chanlocs',eeg.chanlocs,'plot','off');
    delta_hz = f_hz(f_hz>1 & f_hz<4);theta_hz = f_hz(f_hz>4 & f_hz<7);alpha_hz = f_hz(f_hz>8 & f_hz<13);
    lbeta_hz = f_hz(f_hz>13 & f_hz<20);hbeta_hz = f_hz(f_hz>20 & f_hz<30);lgamma_hz = f_hz(f_hz>30 & f_hz<45);
    seg_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
    seg_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
    seg_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
    seg_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
    seg_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
    seg_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);
end

function eeg = seg_chan_bp_specpar(eeg,time_stamps,sel_chan)
    chan_id = cellfun(@(x) ismember(x,sel_chan),{eeg.chanlocs.labels},'UniformOutput',0);chan_idxs = find(cell2mat(chan_id));
    eeg = pop_select(eeg,'time',time_stamps,'channel',sel_chan);

%     eeg = pop_select(eeg,'time',time_stamps);
    eeg = eeg_fooof(eeg,'channel',1:eeg.nbchan,[1 50]);
%     eeg = pop_select(eeg,'channel',sel_chan);

%     [spec_db,f_hz,~,~,~] = spectopo(eeg.data,0,eeg.srate, ...
%                                        'freq',[8 10 12],'freqrange',[0.5 60], ...
%                                        'chanlocs',eeg.chanlocs,'plot','off');
% 
%     specdata = arrayfun(@(y) 10^(y/10), spec_db);specfreqs = f_hz';
%     specpar_out = fooof(specfreqs,specdata(1,:),[1 60],struct(),0);
end

function T = frg_extract_bp(subid,eeg,pre_start_time,short_timestamps,long_timestamps)
    sel_chans = {'FP1','FP2','F7','F8'};
    chan_idx = cellfun(@(x) ismember(x,sel_chans),{eeg.chanlocs.labels},'UniformOutput',0);
    
    pre = pop_select(eeg,'time',[pre_start_time 300],'channel',find(cell2mat(chan_idx)));
    
    [spectra_db,f_hz,~,~,~] = spectopo(pre.data,0,eeg.srate, ...
                                           'freq',[8 10 12],'freqrange',[0.5 60], ...
                                           'chanlocs',pre.chanlocs,'plot','off');
    
%     delta_hz = find(f_hz>1 & f_hz<4);
%     theta_hz = find(f_hz>4 & f_hz<7);
%     alpha_hz = find(f_hz>8 & f_hz<13);
%     lbeta_hz = find(f_hz>13 & f_hz<20);
%     hbeta_hz = find(f_hz>20 & f_hz<30);
%     lgamma_hz = find(f_hz>30 & f_hz<45);
    
    pre_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
    pre_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
    pre_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
    pre_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
    pre_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
    pre_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);

    short = pop_select(eeg,'time',short_timestamps,'channel',find(cell2mat(chan_idx)));

    [spectra_db,f_hz,~,~,~] = spectopo(short.data,0,eeg.srate, ...
                                       'freq',[8 10 12],'freqrange',[0.5 60], ...
                                       'chanlocs',pre.chanlocs,'plot','off');
    
    short_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
    short_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
    short_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
    short_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
    short_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
    short_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);
    
    long = pop_select(eeg,'time',long_timestamps,'channel',find(cell2mat(chan_idx)));

    [spectra_db,f_hz,~,~,~] = spectopo(long.data,0,eeg.srate, ...
                                       'freq',[8 10 12],'freqrange',[0.5 60], ...
                                       'chanlocs',pre.chanlocs,'plot','off');
    
    long_delta_power = 10^(mean(spectra_db(:,delta_hz),'all')/10);
    long_theta_power = 10^(mean(spectra_db(:,theta_hz),'all')/10);
    long_alpha_power = 10^(mean(spectra_db(:,alpha_hz),'all')/10);
    long_lbeta_power = 10^(mean(spectra_db(:,lbeta_hz),'all')/10);
    long_hbeta_power = 10^(mean(spectra_db(:,hbeta_hz),'all')/10);
    long_lgamma_power=10^(mean(spectra_db(:,lgamma_hz),'all')/10);

    T=table(subid,pre_delta_power,pre_theta_power,pre_alpha_power, ...
                  pre_lbeta_power,pre_hbeta_power,pre_lgamma_power,...
                  short_delta_power,short_theta_power,short_alpha_power,...
                  short_lbeta_power,short_hbeta_power,short_lgamma_power,...
                  long_delta_power,long_theta_power,long_alpha_power,...
                  long_lbeta_power,long_hbeta_power,long_lgamma_power);
end