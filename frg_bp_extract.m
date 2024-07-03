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
subid = 15089

%%Comment out for single sub: Exclude badids
% badids = [37532 38058 39862 42125 43543 45194 46037 47678 47744 47801 48238 48278];
% subids = readmatrix(fullfile(data_path,'subids.txt')); subids = setdiff(subids,badids);

%% Main loop begins:
%%Comment out for single sub: Set up loop vars
% for idx = 1:length(subids)
% clearvars -except idx subids data_path
% global subid; subid = subids(idx); 

%% Load data and trim:
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

all_events = {EEG.event.type};
if size(find(ismember(all_events,'TRIGGER EVENT A')),2)>0
    pre_start = find(ismember(all_events,'TRIGGER EVENT A'));
    pre_start_time = EEG.event(pre_start(1)).latency/EEG.srate;
else
    pre_start_time=0;
end
EEG = pop_select(EEG,'time',[pre_start_time/1000 300]);

%Online reference was FCz, which is not there in the data. Adding an
%empty channel at index 19 (same as chanloc file) before loading the channel locations.
EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
EEG = pop_chanedit(EEG,'load','./data/Statnet_F3F4FCz.ced');

%Re-referencing wrt mastoids A1 A2, interpolating FCz (19)
EEG = pop_reref(EEG,{'A1','A2'});EEG = pop_interp(EEG,[origchanlocs(19)]);

%% Filtering and cleaning:
freqs = [0.5 80];
wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,EEG.srate,df);
EEG = pop_firws(EEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);

% EEG = pop_eegfilt(EEG,0.5,80);

EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('noisefreqs',60, ...
                                                 'chunkLength',0,'adaptiveNremove',true, ...
                                                 'fixedNremove',1,'plotResults',0));
% freqs = [0.5 80];
% wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,EEG.srate,df);
% EEG = pop_firws(EEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);

% EEG = pop_eegfilt(EEG,0.5,80);%,'firtype','fir1');

%% ICA:
%Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
thresholds = [0.7 1;0 0.2;0 0.2;0 0.2;0 0.1;0 0.1;0 1];
EEG = pop_runica(EEG,'icatype','runica','extended',1,'interrupt','off');
EEG = iclabel(EEG);
EEG = pop_icflag(EEG,thresholds);
reconsEEG = pop_subcomp(EEG,EEG.reject.gcompreject);

%%
sel_chans = {'FP1','FP2','F7','F8'};
chan_idx = cellfun(@(x) ismember(x,sel_chans),{reconsEEG.chanlocs.labels},'UniformOutput',0);
reconsEEG = pop_select(reconsEEG,'channel',find(cell2mat(chan_idx)));%,'time',[pre_start_time/1000 250]);
disp('Selecting channels:');disp({reconsEEG.chanlocs.labels})

EEG = pop_select(EEG,'channel',find(cell2mat(chan_idx)));%,'time',[pre_start_time/1000 250]);
% vis_artifacts(reconsEEG,EEG)

%% Spectral analysis:
% [spectra_db,f_hz,~,~,~] = pop_spectopo(reconsEEG,1,[0 reconsEEG.xmax*1000],'EEG', ...
%                                        'freq',[8 10 12],'freqrange',[0.5 60],'chanlocs', ...
%                                        reconsEEG.chanlocs,'plot','on','electrodes','labels');

figure; hold on
[spectra_db,f_hz,~,~,~] = spectopo(EEG.data,0,EEG.srate, ...
                                       'freq',[8 10 12],'freqrange',[0.5 60],'chanlocs', ...
                                       EEG.chanlocs,'plot','on');
hold off 
figure; hold on
[spectra_db,f_hz,~,~,~] = spectopo(reconsEEG.data,0,reconsEEG.srate, ...
                                       'freq',[8 10 12],'freqrange',[0.5 60],'chanlocs', ...
                                       reconsEEG.chanlocs,'plot','on');
hold off

% %EEG_in,type_proc,channel_num,tlimits,cycles

% for i = 1:4
% figure; hold on
% [ersp itc powbase times frequencies] = pop_newtimef(EEG,1,i,[],[1 0.5], ...
%                                     'plottype','image','plotersp','on','plotitc','off');
% tftopo(ersp,times,frequencies)
% hold off
% figure; hold on
% [ersp itc powbase times frequencies] = pop_newtimef(reconsEEG,1,i,[],[1 0.5], ...
%                                     'plottype','image','plotersp','on','plotitc','off');
% % hold off
% end

% %Data ranges
% hz8 = find(f_hz>7 & f_hz<9);
% hz7 = find(f_hz>6 & f_hz<8);
% hz6 = find(f_hz>5 & f_hz<7);
% hz5 = find(f_hz>4 & f_hz<6);
% hz4 = find(f_hz>3 & f_hz<5);
% 
% 
% % Compute Power
% hz8powerFz = 10^(mean(spectra_db(hz8))/10);
% hz7powerFz = 10^(mean(spectra_db(hz7))/10);
% hz6powerFz = 10^(mean(spectra_db(hz6))/10);
% hz5powerFz = 10^(mean(spectra_db(hz5))/10);
% hz4powerFz = 10^(mean(spectra_db(hz4))/10);
%% Save preprocessed data:
%%Comment out for single sub
% pop_saveset(EEG,'filename',[num2str(subid) '_preproc'],'filepath',[pwd '/data/clean_set/'],'savemode','onefile');
% end