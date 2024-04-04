%% Foraging EEG analysis:
cd 'C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1';
clear;clc;eeglab;close;
cd 'C:\Users\promitmoitra\Documents\GitHub\frg\'
% data_path = '/home/decision_lab/work/frg/foraging/Neuroflow/';
% badids = [37532 38058 39862 42125 43543 45194 46037 47678 47744 47801 48238 48278];
% %%% badids are [37532 38058 39862 42125 43543 '45194' 46037 '47678' 47744 47801 48238 48278];
% subids = readmatrix(fullfile(data_path,'subids.txt')); subids = setdiff(subids,badids);
% %%% sub 47801 has no poststress long travel time markers TRIGGER EVENT U
% 
subid = 15083;%subids(subids==7873);

%%% Main loop
% for idx = 1:length(subids)
% clearvars -except idx subids data_path
% subid = subids(idx);disp(subid);
% data_file = dir(fullfile(data_path,string(subid),'*EEG','*.edf'));
% data_file_path = fullfile(data_file(end).folder,data_file(end).name);

data_file_path = ['./data/',num2str(subid),'.edf'];
EEG = pop_biosig(data_file_path,'channels',1:19);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data, rereference, seperate pre and post stress blocks:
EEG = pop_select(EEG,'nochannel',10);EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
EEG = pop_chanedit(EEG,'load','./data/Statnet_F3F4FCz.ced');
origchanlocs = readlocs('./data/Statnet_F3F4FCz.ced');

EEG = pop_reref(EEG,19);EEG = pop_reref(EEG,[2 5]);

all_events = {EEG.event.type};
pre_start = find(strcmp(all_events,'TRIGGER EVENT B'));
pre_end = find(strcmp(all_events(pre_start:end),'TRIGGER EVENT Z')); pre_end = pre_end(1)+pre_start;
pre_start_time = EEG.event(pre_start).latency/EEG.srate; pre_end_time = EEG.event(pre_end).latency/EEG.srate;
preEEG = pop_select(EEG,'time',[pre_start_time pre_end_time]);

post_start = find(strcmp(all_events,'TRIGGER EVENT I'));
post_end = find(strcmp(all_events(post_start:end),'TRIGGER EVENT Z')); post_end = post_end(1)+post_start;
post_start_time = EEG.event(post_start).latency/EEG.srate; post_end_time = EEG.event(post_end).latency/EEG.srate;
postEEG = pop_select(EEG,'time',[post_start_time post_end_time]);

%% Exploratory analysis and plots
%%%'condition 2x' check: If urevent after N is neither O nor R, mark epoch
% ur_n = find(strcmp({EEG.urevent.type},'TRIGGER EVENT N'));
% find(~ismember({EEG.urevent(ur_n+1).type},{'TRIGGER EVENT O','TRIGGER EVENT R'}))

%%% Plot stay and leave trial ERPs:
% channel = 'CZ';
% ch_id = find(strcmp({EEG.chanlocs.labels},{channel})); 
% pre = 1;
% if pre
%     pretitle = strcat(string(subid),' pre stress: ', channel);
% else
%     posttitle = strcat(string(subid),' post stress: ', channel);
% end
%
% % if pre
%     preEEG = pop_epoch(preEEG,{'TRIGGER EVENT N'},[-3,3]);
%     % all_N = eeg_getepochevent(preEEG,'type',{'TRIGGER EVENT N'},'fieldname','urevent');
%     leave = zeros(1,preEEG.trials);
%     for i = 1:preEEG.trials
%         e = preEEG.epoch(i).eventtype;
%         if sum(strcmp(e,'TRIGGER EVENT O')) == 0 %|| sum(strcmp(e,'condition 23')) ~= 0
%             leave(i) = 1;
%         end
%     end
%     stay_trials = find(~leave);
%     % leave_trials = find(leave);
%     prestayEEG = pop_select(preEEG,'trial',stay_trials);
%     figure; pop_erpimage(prestayEEG,1,[ch_id],[],strcat(pretitle,' N-O'),10,1,{'TRIGGER EVENT O'},[],'latency','yerplabel','\muV','erp','on','cbar','on')
% % else
%     postEEG = pop_epoch(postEEG,{'TRIGGER EVENT N'},[-3,3]);
%     % all_N = eeg_getepochevent(postEEG,'type',{'TRIGGER EVENT N'},'fieldname','urevent');
%     leave = zeros(1,postEEG.trials);
%     for i = 1:postEEG.trials
%         e = postEEG.epoch(i).eventtype;
%         if sum(strcmp(e,'TRIGGER EVENT O')) == 0 %|| sum(strcmp(e,'condition 23')) ~= 0
%             leave(i) = 1;
%         end
%     end
%     stay_trials = find(~leave);
%     % leave_trials = find(leave);
%     poststayEEG = pop_select(postEEG,'trial',stay_trials);
%     figure; pop_erpimage(poststayEEG,1,[ch_id],[],strcat(posttitle,' N-O'),10,1,{'TRIGGER EVENT O'},[],'latency','yerplabel','\muV','erp','on','cbar','on')

%%% Example code for ICA and various exploratory plots:
% 
% EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','off');
% % preEEG = pop_iclabel(preEEG);
% ALLEEG = [EEG];CURRENTSET = 1;
% [ALLEEG,EEG,CURRENTSET] = processMARA(ALLEEG,EEG,CURRENTSET);
% EEG = pop_selectcomps_MARA(EEG);
% reconsEEG = pop_subcomp(EEG);
% % pop_saveset(reconsEEG,'filename','filt_icasub_7873_post','filepath',pwd,'savemode','onefile');
%
%
% figure; pop_spectopo(preEEG, 1,[0 preEEG.xmax*1e3], 'EEG','freq',[6 10 22],'freqrange',[0.5 60],'electrodes','on');
% figure; pop_spectopo(postEEG,1,[0 postEEG.xmax*1e3],'EEG','freq',[6 10 22],'freqrange',[0.5 60],'electrodes','on');
%
% preEEG = pop_epoch(preEEG,{'TRIGGER EVENT N'},[-3,3]);
% figure; pop_erpimage(preEEG,1,[14],[],'7873 pre cz',10,1,{'TRIGGER EVENT M'},[],'latency','yerplabel','\muV','erp','on','cbar','on')
% 
% figure; pop_erpimage(preEEG,1,[ch_id],[],strcat(pretitle,' O'),10,1,{'TRIGGER EVENT P'},[],'latency','yerplabel','\muV','erp','on','cbar','on')
% saveas(gcf,strcat('/home/decision_lab/work/figs/',string(subid),'_erptrials_pre.png'));
%
% figure; pop_timtopo(preEEG, [-3000  1996], [-500 0 300], '','plotchans',[2 3 14]);
% saveas(gcf,strcat('/home/decision_lab/work/figs/',string(subid),'_erp_pre.png'))
%
% postEEG = pop_epoch(postEEG,{'TRIGGER EVENT N'},[-3,3]);
%
% figure; pop_erpimage(postEEG,1,[14],[],'7873 post cz',10,1,{'TRIGGER EVENT M'},[],'latency','yerplabel','\muV','erp','on','cbar','on')
% figure; pop_erpimage(postEEG,1,[ch_id],[],strcat(posttitle,' O'),10,1,{'TRIGGER EVENT P'},[],'latency','yerplabel','\muV','erp','on','cbar','on')
% saveas(gcf,strcat('/home/decision_lab/work/figs/',string(subid),'_erptrials_post.png'));
% figure; pop_timtopo(postEEG, [-3000  1996], [-500 0 300], '','plotchans',[2 3 14]);
% saveas(gcf,strcat('/home/decision_lab/work/figs/',string(subid),'_erp_post.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Fit ERPs to Fourier series (Order 4) using fit_erp
% % pre_erp = mean(mean(preEEG_epoch.data([2 3 14],:,:),3),1);
% % post_erp = mean(mean(postEEG_epoch.data([2 3 14],:,:),3),1);
% % plot(preEEG_epoch.times, pre_erp); hold on;
% % plot(postEEG_epoch.times, post_erp);
% % for chidx = 1:5
% epoch_data = preEEG_epoch; stress_flag = 'pre';
% init_time = min(epoch_data.times);fin_time = max(epoch_data.times);
% % init_time = 400;fin_time = 1000;
% 
% figure;
% [x,y,y_fit,ph,lh,pl,ll,trange] = fit_erp(epoch_data,init_time,fin_time);
% plot(x,y,'DisplayName','pre'); hold on; plot(x,y_fit,'DisplayName','prefit'); hold on; 
% plot(ll,-pl,'ko','HandleVisibility','off'); hold on; plot(lh,ph,'ko','HandleVisibility','off');hold on;
% 
% title(strcat('resp_locked_erp_',stress_flag,': ',string(subid),trange),'Interpreter','none');%,' Channel: ',epoch_data.chanlocs(chidx).labels));
% hold off;legend;
% % saveas(gcf,strcat('/home/decision_lab/work/figs/resp',trange,'/',stress_flag,'/',string(subid),'_erp',trange,'.png'))%'_',epoch_data.chanlocs(chidx).labels,'_erp.png'))
% % close;
% 
% epoch_data = postEEG_epoch; stress_flag = 'post';
% figure;
% [x,y,y_fit,ph,lh,pl,ll,trange] = fit_erp(epoch_data,init_time,fin_time);
% plot(x,y,'DisplayName','post'); hold on; plot(x,y_fit,'DisplayName','postfit'); hold on; 
% plot(ll,-pl,'ko','HandleVisibility','off'); hold on; plot(lh,ph,'ko','HandleVisibility','off');
% 
% title(strcat('resp_locked_erp_',stress_flag,': ',string(subid),trange),'Interpreter','none');%,' Channel: ',epoch_data.chanlocs(chidx).labels));
% hold off;legend;
% % saveas(gcf,strcat('/home/decision_lab/work/figs/resp',trange,'/',stress_flag,'/',string(subid),'_erp',trange,'.png'))%'_',epoch_data.chanlocs(chidx).labels,'_erp.png'))
% % close;
% % end

% % ERP plots (erpimage - single frequency multi channel):
% freq_rng = [15 30];
% figure('Position', [10 10 900 600]); 
% epoch_data = preEEG_epoch; stress_flag = 'pre';
% % epoch_data = postEEG_epoch; stress_flag = 'post';
% pop_erpimage(epoch_data,1,chan_set,[],strcat(string(subid), ...
%             ': ',stress_flag,{' '},lock_flag,'\_lock\_cpp\_psd'),5,0,{'TRIGGER EVENT O'},[] ,'latency', ...
%             'erp','on','erpalpha',[0.05],'cbar','on','yerplabel','\muV','caxis',[-10 10], ...
%             'nosort','on','coher',cat(2,freq_rng,[0.05]),'plotamps','on','showwin','on')
% %             'nosort','off','phasesort',cat(2,[+500 0],freq_rng),
% % saveas(gcf,strcat('/home/decision_lab/work/figs/resp_locked_ersp/',stress_flag,'/',string(subid),'_ersp.png'))
% % close;
% 
% figure('Position', [10 10 900 600]);
% epoch_data = postEEG_epoch; stress_flag = 'post';
% % epoch_data = preEEG_epoch; stress_flag = 'pre';
% pop_erpimage(epoch_data,1,chan_set,[],strcat(string(subid), ...
%             ': ',stress_flag,{' '},lock_flag,'\_lock\_cpp\_psd'),5,0,{'TRIGGER EVENT O'},[] ,'latency', ...
%             'erp',1,'erpalpha',[0.05],'cbar','on','yerplabel','\muV','caxis',[-10 10], ...    
%              'nosort','on','coher',cat(2,freq_rng,[0.05]),'plotamps','on','showwin','on')
% %             'nosort','off','ampsort',cat(2,[+1000 0],freq_rng),
% % saveas(gcf,strcat('/home/decision_lab/work/figs/resp_locked_ersp/',stress_flag,'/',string(subid),'_ersp.png'))
% % close;

% % Time-Frequency plots (pop_newtimef - single channel full spectrum)
% figure('Position',[73 615 1842 347])
% tlo = tiledlayout(1,5,'TileSpacing','loose','Padding','compact','Position', [0.15,0.15,0.5,0.5]);
% leave_trial_idx = get_leave_idx(postEEG_epoch); stress_flag = 'post';
% leave_trial_idx(2:end+1)=leave_trial_idx(1:end); leave_trial_idx(1)=0;
% min_patch_len = min(min(diff(leave_trial_idx)-1),4); % % By default, plot from leave-4 to leave, 
%                                                      % % unless there's a patch with less than 3 stay trials
% for i = min_patch_len:-1:0
%     epoch_data = postEEG_epoch; 
%     leave_trial_idx = get_leave_idx(epoch_data);
%     epoch_data.data = epoch_data.data(:,:,leave_trial_idx-i);
%     ax = nexttile(tlo);
%     
% %   [ersp, itc, powbase, times, frequencies,~,~,tfdata] = ...
%     pop_newtimef(epoch_data,1,chan_idx,[],[1 0.5], ...
%                 'plottype','image','plotersp','off','plotitc','on', ...
%                 'plotphase','off','plotphaseonly','off','erspmax',3, ...
%                 'baseline',baseline,'basenorm','off','trialbase','full', ...
%                 'freqs',freq_rng,'padratio',32);%,'title',fig_title);%, ...
% %                 'pcontour','off','alpha',0.01);
%     fname = sprintf('/home/decision_lab/work/figs/itc/%s/%s/%d_itc.png',lock_flag,stress_flag,subid);
%     saveas(gcf,fname);
% end
% close;

% % Time-Frequency plots (newtimef - relative to condition (left (stay1) or right (leave) alligned))
% epoch_data = preEEG_epoch; stress_flag = 'pre';
% epoch_data = postEEG_epoch; stress_flag = 'post';
% freq_rng = [1 50];
% chan = 'CZ';chan_idx = find(strcmp({EEG.chanlocs.labels},{chan}));
% caption_txt = char(strcat(string(subid),'_',stress_flag,'-stress_', ...
%               lock_flag,'-locked_',chan));
% 
% leave_trial_idx = get_leave_idx(epoch_data);
%%% get_leave_idx is called again inside the loop because first idx is set 
%%% to 0 manually, to calculate the min_stay_len.
% first_stay_idx = [1 leave_trial_idx(1:end-1)+1];
% leave_trial_idx(2:end+1)=leave_trial_idx(1:end); leave_trial_idx(1)=0;
% min_stay_len = min(min(diff(leave_trial_idx)-1),4); 
%%% By default, plot from leave-4 to leave, unless there's a patch with 
%%% less than 3 stay trials
%
% for i = min_stay_len:-1:0
%     leave_trial_idx = get_leave_idx(epoch_data);
%     first_stay_data = epoch_data.data(chan_idx,:,first_stay_idx); %cond = 'Stay 1';
%     leave_data = epoch_data.data(chan_idx,:,leave_trial_idx); %cond = 'Leave';
% 
%     stay_alligned_data = epoch_data.data(chan_idx,:,first_stay_idx+i); %fig_title = {[strcat('Stay 1 + ',string(i))],cond};
%     leave_alligned_data = epoch_data.data(chan_idx,:,leave_trial_idx-i); %fig_title = {[strcat('Leave - ',string(i))],cond};
% 
%     [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] =...
%     newtimef({leave_alligned_data leave_data},size(leave_data,2),...
%               epoch_trange*1000,preEEG_epoch.srate,[3 0.5],...
%               'plotitc','off',...
%               'title',fig_title,'caption',caption_txt,...
%               'scale','abs','baseline',NaN,'basenorm','off','commonbase','off','trialbase','off');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean line noise and bandpass filter:
preEEG = clean_data_with_zapline_plus_eeglab_wrapper(preEEG,struct('noisefreqs',60,'chunkLength',0,'adaptiveNremove',true, ...
                                                    'fixedNremove',1,'plotResults',0));
postEEG = clean_data_with_zapline_plus_eeglab_wrapper(postEEG,struct('noisefreqs',60,'chunkLength',0,'adaptiveNremove',true, ...
                                                    'fixedNremove',1,'plotResults',0));

freqs = [1 100];
wtype = 'hamming'; df = 1; m = pop_firwsord(wtype,EEG.srate,df);
preEEG = pop_firws(preEEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);
postEEG = pop_firws(postEEG,'wtype',wtype,'ftype','bandpass','fcutoff',freqs,'forder',m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract travel-time blocks from markers:
pre_events = {preEEG.event.type};
pre_short_end = find(strcmp(pre_events,'TRIGGER EVENT T'));pre_short_end = pre_short_end(end);
pre_long_end = find(strcmp(pre_events,'TRIGGER EVENT U'));pre_long_end = pre_long_end(end);
post_events = {postEEG.event.type};
post_short_end = find(strcmp(post_events,'TRIGGER EVENT T'));post_short_end = post_short_end(end);
post_long_end = find(strcmp(post_events,'TRIGGER EVENT U'));post_long_end = post_long_end(end);

pre_short_end_time = preEEG.event(pre_short_end).latency/EEG.srate;
pre_long_end_time = preEEG.event(pre_long_end).latency/EEG.srate;
post_short_end_time = postEEG.event(post_short_end).latency/EEG.srate;
post_long_end_time = postEEG.event(post_long_end).latency/EEG.srate;

if pre_short_end<pre_long_end
    preshortEEG = pop_select(preEEG,'time',[0 pre_short_end_time+15]);
    prelongEEG = pop_select(preEEG,'time',[pre_short_end_time+15 pre_long_end_time]);
else
    prelongEEG = pop_select(preEEG,'time',[0 pre_long_end_time+15]);
    preshortEEG = pop_select(preEEG,'time',[pre_long_end_time+15 pre_short_end_time]);
end

if post_short_end<post_long_end
    postshortEEG = pop_select(postEEG,'time',[0 post_short_end_time+15]);
    postlongEEG = pop_select(postEEG,'time',[post_short_end_time+15 post_long_end_time]);
else
    postlongEEG = pop_select(postEEG,'time',[0 post_long_end_time+15]);
    postshortEEG = pop_select(postEEG,'time',[post_long_end_time+15 post_short_end_time]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Epoch locked to stimulus (Trigger Event N) or response (Trigger Events [O, R]):
lock_event = {'TRIGGER EVENT N'}; lock_flag = 'stim'; 
epoch_trange = [-1 2]; %baseline = [-1 0]*1000;

% lock_event = {'TRIGGER EVENT O','TRIGGER EVENT R'}; lock_flag = 'resp';
% epoch_trange = [-2 1]; resp_baseline = [-2 -1]*1000;

preEEG_epoch = pop_epoch(preEEG,lock_event,epoch_trange);
postEEG_epoch = pop_epoch(postEEG,lock_event,epoch_trange);

preshort_epoch = pop_epoch(preshortEEG,lock_event,epoch_trange);
prelong_epoch = pop_epoch(prelongEEG,lock_event,epoch_trange);
postshort_epoch = pop_epoch(postshortEEG,lock_event,epoch_trange);
postlong_epoch = pop_epoch(postlongEEG,lock_event,epoch_trange);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SpecParam
% % epoch_data = preEEG_epoch; stress_flag='pre'; tt_flag='';
% % epoch_data = postEEG_epoch; stress_flag='post'; tt_flag='';
% 
% epoch_data = preshort_epoch;   stress_flag='pre'; tt_flag='short';
% % epoch_data = prelong_epoch;  stress_flag='pre'; tt_flag='long';
% % epoch_data = postshort_epoch;stress_flag='post'; tt_flag='short';
% % epoch_data = postlong_epoch; stress_flag='post'; tt_flag='long';
% 
% channel = 'CZ';cz_idx = find(strcmp({EEG.chanlocs.labels},{channel}));
% % channel = 'C3';c3_idx = find(strcmp({EEG.chanlocs.labels},{channel}));
% % channel = 'C4';c4_idx = find(strcmp({EEG.chanlocs.labels},{channel})); 
% % channel = 'FP1';fp1_idx = find(strcmp({EEG.chanlocs.labels},{channel}));
% % channel = 'FP2';fp2_idx = find(strcmp({EEG.chanlocs.labels},{channel}));
% % sma_chan_set = [c3_idx cz_idx c4_idx];
% % fp_chan_set = [fp1_idx fp2_idx];
% 
% chan_idx = cz_idx;
% 
% leave_trial_idxs = get_leave_idx(epoch_data);
% first_stay_idxs = [1 leave_trial_idxs(1:end-1)+1];
% 
% patch_trial_idxs = arrayfun(@(f,g) (f:g),first_stay_idxs,leave_trial_idxs,'UniformOutput',false);
% num_patches = size(patch_trial_idxs,2);
% patch_trial_len = cell2mat(cellfun(@length,patch_trial_idxs,'UniformOutput',false));
% 
% %%%TODO: Add code to fooof chan_sets
% stay_lock_trial_idxs = reshape(cell2mat(cellfun(@(x) [x(1:end-1) nan(1,max(patch_trial_len)-size(x(1:end-1),2)-1)], ...
%                                     patch_trial_idxs,'UniformOutput',false)),max(patch_trial_len)-1,num_patches)';
% leave_lock_trial_idxs = reshape(cell2mat(cellfun(@(x) [nan(1,max(patch_trial_len)-size(x,2)) x], ...
%                                     patch_trial_idxs,'UniformOutput',false)),max(patch_trial_len),num_patches)';
% 
% fooof_settings = struct('peak_width_limits',[2 10],'max_n_peaks',inf,'min_peak_height',0, ...
%                         'peak_threshold',2,'aperiodic_mode','fixed','verbose',false);
% 
% stay_lock_slopes = nan(num_patches+1,max(patch_trial_len)-1);
% stay_lock_delta_bp  = nan(num_patches+1,max(patch_trial_len)-1);
% stay_lock_theta_bp  = nan(num_patches+1,max(patch_trial_len)-1);
% stay_lock_alpha_bp  = nan(num_patches+1,max(patch_trial_len)-1);
% stay_lock_beta_bp   = nan(num_patches+1,max(patch_trial_len)-1);
% stay_lock_gamma_bp  = nan(num_patches+1,max(patch_trial_len)-1);
% 
% leave_lock_slopes = nan(num_patches+1,max(patch_trial_len));
% leave_lock_delta_bp  = nan(num_patches+1,max(patch_trial_len));
% leave_lock_theta_bp  = nan(num_patches+1,max(patch_trial_len));
% leave_lock_alpha_bp  = nan(num_patches+1,max(patch_trial_len));
% leave_lock_beta_bp   = nan(num_patches+1,max(patch_trial_len));
% leave_lock_gamma_bp  = nan(num_patches+1,max(patch_trial_len));
% 
% %alligned stay_lock and leave_lock 
% for idx=1:max(patch_trial_len)
%     if idx ~= max(patch_trial_len)
%         eval(strcat("stay_",num2str(idx),"= pop_select(epoch_data,'trial',stay_lock_trial_idxs(~isnan(stay_lock_trial_idxs(:,",num2str(idx),")),",num2str(idx),"));"));
%         eval(strcat("stay_",num2str(idx),".xmin=-1;stay_",num2str(idx),".xmax=+1.9960;"))
%         eval(strcat("stay_lock_fooof_prestim = eeg_fooof(stay_",num2str(idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"));
%     %     eval(strcat("stay_lock_fooof_prestim = eeg_fooof(stay_",num2str(idx),",'channel',[1:epoch_data.nbchan],[-1 0]*1000,100,freqs,fooof_settings);"));
%         fooof_prestim = cell2mat(stay_lock_fooof_prestim.etc.FOOOF_results(chan_idx));
%     %     stay_lock_slopes(1,idx) = fooof_prestim(chan_idx).aperiodic_params(end);
%         stay_lock_slopes(end,idx) = fooof_prestim.aperiodic_params(end);
%         stay_lock_delta_bp(end,idx) = mean(fooof_prestim.bandpowers(1),'omitnan');
%         stay_lock_theta_bp(end,idx) = mean(fooof_prestim.bandpowers(2),'omitnan');
%         stay_lock_alpha_bp(end,idx) = mean(fooof_prestim.bandpowers(3),'omitnan');
%         stay_lock_beta_bp( end,idx) = mean(fooof_prestim.bandpowers(4),'omitnan');
%         stay_lock_gamma_bp(end,idx) = mean(fooof_prestim.bandpowers(5),'omitnan');
%     end
%     eval(strcat("leave_",num2str(max(patch_trial_len)-idx),"= pop_select(epoch_data,'trial',leave_lock_trial_idxs(~isnan(leave_lock_trial_idxs(:,", ...
%                 num2str(idx),")),",num2str(idx),"));"));
%     eval(strcat("leave_",num2str(max(patch_trial_len)-idx),".xmin=-1;leave_",num2str(max(patch_trial_len)-idx),".xmax=1.9960;"))
%     eval(strcat("leave_lock_fooof_prestim = eeg_fooof(leave_",num2str(max(patch_trial_len)-idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"));
% %     eval(strcat("leave_lock_fooof_prestim = eeg_fooof(leave_",num2str(max(patch_trial_len)-idx),",'channel',[1:epoch_data.nbchan],[-1 0]*1000,100,freqs,fooof_settings);"));
%     fooof_prestim = cell2mat(leave_lock_fooof_prestim.etc.FOOOF_results(chan_idx));
% %     leave_lock_slopes(1,idx) = fooof_prestim(chan_idx).aperiodic_params(end);
%     leave_lock_slopes(end,end-idx+1) = fooof_prestim.aperiodic_params(end);
%     leave_lock_delta_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(1),'omitnan');
%     leave_lock_theta_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(2),'omitnan');
%     leave_lock_alpha_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(3),'omitnan');
%     leave_lock_beta_bp( end,end-idx+1) = mean(fooof_prestim.bandpowers(4),'omitnan');
%     leave_lock_gamma_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(5),'omitnan');
% end
% %process trial-wise
% for patch_idx=1:num_patches
%     first_stay_idx=1;leave_trial_idx=patch_trial_len(patch_idx);
%     patch = pop_select(epoch_data,'trial',patch_trial_idxs{patch_idx});
%     patch.xmin=-1;patch.xmax=1.996;
% 
%     for idx = first_stay_idx:leave_trial_idx
%         if idx~=leave_trial_idx
%             eval(strcat("stay_",num2str(patch_idx),"_",num2str(idx),"= pop_select(patch,'trial',idx);"));
%             eval(strcat("stay_",num2str(patch_idx),"_",num2str(idx),".xmin=-1;stay_",num2str(patch_idx),"_",num2str(idx),".xmax=1.9660;"));
%             eval(strcat("stay_fooof_prestim = eeg_fooof(stay_",num2str(patch_idx),"_",num2str(idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"))
%             fooof_prestim = cell2mat(stay_fooof_prestim.etc.FOOOF_results(chan_idx));
%             stay_lock_slopes(patch_idx,idx) = fooof_prestim.aperiodic_params(end);
%             stay_lock_delta_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(1),'omitnan');
%             stay_lock_theta_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(2),'omitnan');
%             stay_lock_alpha_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(3),'omitnan');
%             stay_lock_beta_bp(  patch_idx,idx) = mean(fooof_prestim.bandpowers(4),'omitnan');
%             stay_lock_gamma_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(5),'omitnan');
%         end
%         eval(strcat("leave_",num2str(patch_idx),"_",num2str(idx-1),"= pop_select(patch,'trial',leave_trial_idx-idx+1);"));
%         eval(strcat("leave_",num2str(patch_idx),"_",num2str(idx-1),".xmin=-1;leave_",num2str(patch_idx),"_",num2str(idx-1),".xmax=1.9660;"));
%         eval(strcat("leave_fooof_prestim = eeg_fooof(leave_",num2str(patch_idx),"_",num2str(idx-1),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"))
%         fooof_prestim = cell2mat(leave_fooof_prestim.etc.FOOOF_results(chan_idx));
%         leave_lock_slopes(patch_idx,end-idx+1) = fooof_prestim.aperiodic_params(end);
%         leave_lock_delta_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(1),'omitnan');
%         leave_lock_theta_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(2),'omitnan');
%         leave_lock_alpha_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(3),'omitnan');
%         leave_lock_beta_bp(  patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(4),'omitnan');
%         leave_lock_gamma_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(5),'omitnan');        
%     end
%     % legend;
%     % ax=gca;ax.XScale='log';ax.YScale='log';
% end
% fprintf('SpecParam Done!\n')

% stress_flag='pre'; tt_flag='short';
% stress_flag='pre'; tt_flag='long';
% stress_flag='post'; tt_flag='short';
stress_flag='post'; tt_flag='long';

eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
channel = 'CZ';cz_idx = find(strcmp({EEG.chanlocs.labels},{channel}));
chan_idx = cz_idx;

leave_trial_idxs = get_leave_idx(epoch_data);
first_stay_idxs = [1 leave_trial_idxs(1:end-1)+1];
patch_trial_idxs = arrayfun(@(f,g) (f:g),first_stay_idxs,leave_trial_idxs,'UniformOutput',false);
num_patches = size(patch_trial_idxs,2);
patch_trial_len = cell2mat(cellfun(@length,patch_trial_idxs,'UniformOutput',false));

[stay_lock_res, leave_lock_res] = specparam(epoch_data,chan_idx,freqs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select, trim and concatenate data for plotting
% data_flag = "slopes"; ylab = "Aperiodic exponent";
% data_flag = "delta_bp"; ylab = "Delta band power";
% data_flag = "theta_bp"; ylab = "Theta band power";
% data_flag = "alpha_bp"; ylab = "Alpha band power";
% data_flag = "beta_bp"; ylab = "Beta band power";
% data_flag = "gamma_bp"; ylab = "Gamma band power";
%%% Trimming trials which are present in less than a fraction of patches
% eval(strcat("stay_lock_data = stay_lock_res.",data_flag,";"));
% frac_nan=sum(isnan(stay_lock_data),1)/(size(stay_lock_data,1)-1);
% stay_lock_data=stay_lock_data(:,frac_nan<=0.25);
% eval(strcat("leave_lock_data = leave_lock_res.",data_flag,";"));
% frac_nan=sum(isnan(leave_lock_data),1)/(size(leave_lock_data,1)-1);
% leave_lock_data=leave_lock_data(:,frac_nan<=0.25);
% 
% plot_data = [stay_lock_data leave_lock_data];
% eval(strcat("[",data_flag,"_err,",data_flag,"]=std(plot_data(1:end-1,:),1,1,'omitnan');"));
% eval(strcat("err_data=",data_flag,"_err;"));
% 
% fig=figure('Name',strcat(num2str(subid)," ",channel," ",stress_flag," ",tt_flag)); hold on;
% title(strcat(num2str(subid)," ",channel," ",stress_flag," ",tt_flag));
% 
% x = [1:size(stay_lock_data,2)+size(leave_lock_data,2)];
% xticks(1:size(stay_lock_data,2)+size(leave_lock_data,2));
% xticklabels([1:size(stay_lock_data,2) size(leave_lock_data,2)-1:-1:0]);
% xlabel("Trials (Stay1 + i | Leave - i)");
% ylabel(ylab);
% 
% xline(size(stay_lock_data,2)+0.5,'--k','LineWidth',2)
% displaynames = arrayfun(@(x)sprintf('Patch %d',x),1:size(plot_data,1)-1,'uni',0);
% displaynames = [displaynames {'Patch Average'}];
% p = plot(x,plot_data(1:end-1,:),'Marker','+','LineStyle',':','LineWidth',1);
% p = [p; plot(x,plot_data(end,:),'Marker','o','LineWidth',2.5,'Color','#7E2F8E')];hold on;
% errorbar(x,plot_data(end,:),err_data,'LineWidth',1,'Color','#7E2F8E');hold on;
% 
% legend(p,displaynames)
% hold off;

% data_flag = "slopes"; ylab = "Aperiodic exponent";
% data_flag = "delta_bp"; ylab = "Delta band power";
% data_flag = "theta_bp"; ylab = "Theta band power";
data_flag = "alpha_bp"; ylab = "Alpha band power";
% data_flag = "beta_bp"; ylab = "Beta band power";
% data_flag = "gamma_bp"; ylab = "Gamma band power";

eval(strcat("stay_lock_data = stay_lock_res.",data_flag,";"));
frac_nan=sum(isnan(stay_lock_data),1)/(size(stay_lock_data,1)-1);
stay_lock_data=stay_lock_data(:,frac_nan<=0.25);
stay_lock_data=stay_lock_data(:,1:int8(size(stay_lock_data,2)/2));

eval(strcat("leave_lock_data = leave_lock_res.",data_flag,";"));
frac_nan=sum(isnan(leave_lock_data),1)/(size(leave_lock_data,1)-1);
leave_lock_data=leave_lock_data(:,frac_nan<=0.25);
leave_lock_data=leave_lock_data(:,int8(size(leave_lock_data,2)/2):end);

plot_data = [stay_lock_data leave_lock_data];

fig=frg_plot(num_patches,data_flag,ylab,plot_data,stay_lock_data,leave_lock_data);
figname = strcat(num2str(subid)," ",channel," ",stress_flag," ",tt_flag);
fig.Name= figname;
title(figname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% poststim = eeg_fooof(preEEG_epoch,"channel",[1:EEG.nbchan],[0 1]*1000,100,freqs,fooof_settings);
% fooof_poststim = cell2mat(poststim.etc.FOOOF_results);
% postresp = eeg_fooof(preEEG_epoch,"channel",[1:EEG.nbchan],[1 2]*1000,100,freqs,fooof_settings);
% fooof_postresp = cell2mat(postresp.etc.FOOOF_results);
% fooof_results = [fooof_prestim fooof_poststim fooof_postresp];
% peak_data = {fooof_results(3,:).peak_params};
% %%
% figure('Name',EEG.chanlocs(chan_idx).labels);hold on;
% p1 = plot(fooof_results(chan_idx,1).freqs,10.^(fooof_results(chan_idx,1).fooofed_spectrum),'DisplayName','[-2 0] ms','LineWidth',2,'Color',"#0072BD");
% plot(fooof_results(chan_idx,1).freqs,10.^(fooof_results(chan_idx,1).ap_fit),'LineWidth',1,'LineStyle','--','Color',p1.Color);
% xline(peak_data{1}(:,1),'LineStyle',':','Color',p1.Color,'LineWidth',2);
% 
% p2 = plot(fooof_results(chan_idx,2).freqs,10.^(fooof_results(chan_idx,2).fooofed_spectrum),'DisplayName','[0 +1] ms','LineWidth',2,'Color',"#D95319");
% plot(fooof_results(chan_idx,2).freqs,10.^(fooof_results(chan_idx,2).ap_fit),'LineWidth',2,'LineStyle','--','Color',p2.Color);
% xline(peak_data{2}(:,1),'LineStyle',':','Color',p2.Color,'LineWidth',2);
% 
% p3 = plot(fooof_results(chan_idx,3).freqs,10.^(fooof_results(chan_idx,3).fooofed_spectrum),'DisplayName','[+1 +2] ms','LineWidth',2,'Color',"#EDB120");
% plot(fooof_results(chan_idx,3).freqs,10.^(fooof_results(chan_idx,3).ap_fit),'LineWidth',2,'LineStyle','--','Color',p3.Color);
% xline(peak_data{3}(:,1),'LineStyle',':','Color',p3.Color,'LineWidth',2);
% 
% ax=gca;ax.XScale='linear';ax.YScale='log';
% fmarks=cell2mat(peak_data(:));fmarks=sort(fmarks(:,1));fdiff=diff(fmarks);
% fmarks(fdiff<1)=[];xticks(fmarks);xticklabels(cellstr(num2str(fmarks(:),'%.1f')));
% yticklabels(cellstr(num2str(ax.YTick(:),'%.2f')));
% 
% xlabel('Frequency (Hz)');ylabel('PSD (\muV^2/Hz)')
% title(EEG.chanlocs(chan_idx).labels);legend([p1,p2,p3]);
% % end
%% Dataframes for EEGNet and HDDM:
% save(strcat('/home/decision_lab/work/eegnet/arl-eegmodels/',str(subid),'_pre.mat'),'-struct','preEEG_epoch')
% save(strcat('/home/decision_lab/work/eegnet/arl-eegmodels/',str(subid),'_post.mat'),'-struct','postEEG_epoch')
% Dataframe for hddm:
% preT = gendf(subid,preEEG_epoch,1,pre_short_end,pre_long_end);
% postT = gendf(subid,postEEG_epoch,0,post_short_end,post_long_end);
% subT = [preT;postT];
% title = strcat('/home/decision_lab/work/data_hddm/l0_eeg_hddm/raw/',string(subid),'_eeg_hddm.csv');
% writetable(subT,title);
% end
%% scratchpad
%     theta_peak_locs = find(fooof_prestim(chan_idx).peak_params(:,1)>2 & fooof_prestim(chan_idx).peak_params(:,1)<6);
%     if ~isempty(peak_locs)
%         stay_lock_patch_bandpower(1,idx) = mean(fooof_prestim(chan_idx).peak_params(peak_locs,2));
%     end
%     for patch_idx=1:num_patches
%     [s,m]=std(patch_bandpower,1,1,'omitnan');
%     p=errorbar(1:max(patch_trial_len),m,s, ...
%       'Marker','o','LineWidth',2.5,'DisplayName',"Patch Avg");hold on;
%     % for trial_idx=1:patch_trial_len(patch_idx)
%     %     p(patch_idx)=plot(1:max(patch_trial_len),patch_slopes(patch_idx,:), ...
%     %                     'LineWidth',2.5,'DisplayName',strcat("Patch ",num2str(patch_idx)));hold on;
%     % end
% end
% c=parula(min_stay_len+1);%c(1:end)=c(end:-1:1);
% patch_trials=single(~cellfun('isempty',patch_slopes));patch_trials(patch_trials==0)=nan;
% [[patch_slopes{patch_idx,:}], NaN(1,size(patch_slopes,2)-length([patch_slopes{patch_idx,:}]))];
% 
% leave_trial_idx(2:end+1)=leave_trial_idx(1:end);leave_trial_idx(1)=0;
% min_stay_len = min(min(diff(leave_trial_idx)-1),6);
% leave_trial_idx = get_leave_idx(epoch_data);
% patch_slopes(patch_idx,idx+max(patch_trial_len)-(2*idx-1)) %right allign, left pad
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%==Function Definitions==%%%%%%%%%%%%%%%%%%%%%%%%%
function leave_idx = get_leave_idx(ep_dat)
    leave_idx = []; 
    trigger_events = {'TRIGGER EVENT R'};
    for i = 1:ep_dat.trials
        trial_events = ep_dat.epoch(i).eventtype;
        [~,event_idx] = ismember(trial_events,trigger_events);
        if ismember(1,event_idx)
            leave_idx(end+1) = i;
        end
    end
end

function [stay_lock_res, leave_lock_res] = specparam(ep_dat,chan_idx,freqs)
    leave_trial_idxs = get_leave_idx(ep_dat);
    first_stay_idxs = [1 leave_trial_idxs(1:end-1)+1];
    patch_trial_idxs = arrayfun(@(f,g) (f:g),first_stay_idxs,leave_trial_idxs,'UniformOutput',false);
    num_patches = size(patch_trial_idxs,2);
    patch_trial_len = cell2mat(cellfun(@length,patch_trial_idxs,'UniformOutput',false));

    stay_lock_trial_idxs = reshape(cell2mat(cellfun(@(x) [x(1:end-1) nan(1,max(patch_trial_len)-size(x(1:end-1),2)-1)], ...
                                        patch_trial_idxs,'UniformOutput',false)), ...
                                        max(patch_trial_len)-1,num_patches)';
    leave_lock_trial_idxs = reshape(cell2mat(cellfun(@(x) [nan(1,max(patch_trial_len)-size(x,2)) x], ...
                                        patch_trial_idxs,'UniformOutput',false)), ...
                                        max(patch_trial_len),num_patches)';
    
    fooof_settings = struct('peak_width_limits',[2 10],'max_n_peaks',inf,'min_peak_height',0, ...
                            'peak_threshold',2,'aperiodic_mode','fixed','verbose',false);

    stay_lock_slopes = nan(num_patches+1,max(patch_trial_len)-1);
    stay_lock_delta_bp  = nan(num_patches+1,max(patch_trial_len)-1);
    stay_lock_theta_bp  = nan(num_patches+1,max(patch_trial_len)-1);
    stay_lock_alpha_bp  = nan(num_patches+1,max(patch_trial_len)-1);
    stay_lock_beta_bp   = nan(num_patches+1,max(patch_trial_len)-1);
    stay_lock_gamma_bp  = nan(num_patches+1,max(patch_trial_len)-1);
    
    leave_lock_slopes = nan(num_patches+1,max(patch_trial_len));
    leave_lock_delta_bp  = nan(num_patches+1,max(patch_trial_len));
    leave_lock_theta_bp  = nan(num_patches+1,max(patch_trial_len));
    leave_lock_alpha_bp  = nan(num_patches+1,max(patch_trial_len));
    leave_lock_beta_bp   = nan(num_patches+1,max(patch_trial_len));
    leave_lock_gamma_bp  = nan(num_patches+1,max(patch_trial_len));

    % patch trials alligned stay_lock and leave_lock 
    for idx=1:max(patch_trial_len)
        if idx ~= max(patch_trial_len)
            eval(strcat("stay_",num2str(idx),"= pop_select(ep_dat,'trial',stay_lock_trial_idxs(~isnan(stay_lock_trial_idxs(:,",num2str(idx),")),",num2str(idx),"));"));
            eval(strcat("stay_",num2str(idx),".xmin=-1;stay_",num2str(idx),".xmax=+1.9960;"))
            eval(strcat("stay_lock_fooof_prestim = eeg_fooof(stay_",num2str(idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"));
        %     eval(strcat("stay_lock_fooof_prestim = eeg_fooof(stay_",num2str(idx),",'channel',[1:ep_dat.nbchan],[-1 0]*1000,100,freqs,fooof_settings);"));
            fooof_prestim = cell2mat(stay_lock_fooof_prestim.etc.FOOOF_results(chan_idx));
        %     stay_lock_slopes(1,idx) = fooof_prestim(chan_idx).aperiodic_params(end);
            stay_lock_slopes(end,idx) = fooof_prestim.aperiodic_params(end);
            stay_lock_delta_bp(end,idx) = mean(fooof_prestim.bandpowers(1),'omitnan');
            stay_lock_theta_bp(end,idx) = mean(fooof_prestim.bandpowers(2),'omitnan');
            stay_lock_alpha_bp(end,idx) = mean(fooof_prestim.bandpowers(3),'omitnan');
            stay_lock_beta_bp( end,idx) = mean(fooof_prestim.bandpowers(4),'omitnan');
            stay_lock_gamma_bp(end,idx) = mean(fooof_prestim.bandpowers(5),'omitnan');
        end
        eval(strcat("leave_",num2str(max(patch_trial_len)-idx),"= pop_select(ep_dat,'trial',leave_lock_trial_idxs(~isnan(leave_lock_trial_idxs(:,", ...
                    num2str(idx),")),",num2str(idx),"));"));
        eval(strcat("leave_",num2str(max(patch_trial_len)-idx),".xmin=-1;leave_",num2str(max(patch_trial_len)-idx),".xmax=1.9960;"))
        eval(strcat("leave_lock_fooof_prestim = eeg_fooof(leave_",num2str(max(patch_trial_len)-idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"));
    %     eval(strcat("leave_lock_fooof_prestim = eeg_fooof(leave_",num2str(max(patch_trial_len)-idx),",'channel',[1:ep_dat.nbchan],[-1 0]*1000,100,freqs,fooof_settings);"));
        fooof_prestim = cell2mat(leave_lock_fooof_prestim.etc.FOOOF_results(chan_idx));
    %     leave_lock_slopes(1,idx) = fooof_prestim(chan_idx).aperiodic_params(end);
        leave_lock_slopes(end,end-idx+1) = fooof_prestim.aperiodic_params(end);
        leave_lock_delta_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(1),'omitnan');
        leave_lock_theta_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(2),'omitnan');
        leave_lock_alpha_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(3),'omitnan');
        leave_lock_beta_bp( end,end-idx+1) = mean(fooof_prestim.bandpowers(4),'omitnan');
        leave_lock_gamma_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(5),'omitnan');
    end
    % single trials
    for patch_idx=1:num_patches
        first_stay_idx=1;leave_trial_idx=patch_trial_len(patch_idx);
        patch = pop_select(ep_dat,'trial',patch_trial_idxs{patch_idx});
        patch.xmin=-1;patch.xmax=1.996;
    
        for idx = first_stay_idx:leave_trial_idx
            if idx~=leave_trial_idx
                eval(strcat("stay_",num2str(patch_idx),"_",num2str(idx),"= pop_select(patch,'trial',idx);"));
                eval(strcat("stay_",num2str(patch_idx),"_",num2str(idx),".xmin=-1;stay_",num2str(patch_idx),"_",num2str(idx),".xmax=1.9660;"));
                eval(strcat("stay_fooof_prestim = eeg_fooof(stay_",num2str(patch_idx),"_",num2str(idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"))
                fooof_prestim = cell2mat(stay_fooof_prestim.etc.FOOOF_results(chan_idx));
                stay_lock_slopes(patch_idx,idx) = fooof_prestim.aperiodic_params(end);
                stay_lock_delta_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(1),'omitnan');
                stay_lock_theta_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(2),'omitnan');
                stay_lock_alpha_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(3),'omitnan');
                stay_lock_beta_bp(  patch_idx,idx) = mean(fooof_prestim.bandpowers(4),'omitnan');
                stay_lock_gamma_bp( patch_idx,idx) = mean(fooof_prestim.bandpowers(5),'omitnan');
            end
            eval(strcat("leave_",num2str(patch_idx),"_",num2str(idx-1),"= pop_select(patch,'trial',leave_trial_idx-idx+1);"));
            eval(strcat("leave_",num2str(patch_idx),"_",num2str(idx-1),".xmin=-1;leave_",num2str(patch_idx),"_",num2str(idx-1),".xmax=1.9660;"));
            eval(strcat("leave_fooof_prestim = eeg_fooof(leave_",num2str(patch_idx),"_",num2str(idx-1),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"))
            fooof_prestim = cell2mat(leave_fooof_prestim.etc.FOOOF_results(chan_idx));
            leave_lock_slopes(patch_idx,end-idx+1) = fooof_prestim.aperiodic_params(end);
            leave_lock_delta_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(1),'omitnan');
            leave_lock_theta_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(2),'omitnan');
            leave_lock_alpha_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(3),'omitnan');
            leave_lock_beta_bp(  patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(4),'omitnan');
            leave_lock_gamma_bp( patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(5),'omitnan');        
        end
    end

    stay_lock_res = struct();
    stay_lock_res.slopes  =  stay_lock_slopes;
    stay_lock_res.delta_bp = stay_lock_delta_bp;
    stay_lock_res.theta_bp = stay_lock_theta_bp;
    stay_lock_res.alpha_bp = stay_lock_alpha_bp;
    stay_lock_res.beta_bp  = stay_lock_beta_bp;
    stay_lock_res.gamma_bp = stay_lock_gamma_bp;

    leave_lock_res = struct();
    leave_lock_res.slopes  =  leave_lock_slopes;
    leave_lock_res.delta_bp = leave_lock_delta_bp;
    leave_lock_res.theta_bp = leave_lock_theta_bp;
    leave_lock_res.alpha_bp = leave_lock_alpha_bp;
    leave_lock_res.beta_bp  = leave_lock_beta_bp;
    leave_lock_res.gamma_bp = leave_lock_gamma_bp;

    fprintf('SpecParam Done!\n')
end    

function fig = frg_plot(num_patches,data_flag,ylab,plot_data,stay_lock_data,leave_lock_data)
    eval(strcat("[",data_flag,"_err,",data_flag,"]=std(plot_data(1:end-1,:),1,1,'omitnan');"));
    eval(strcat("err_data=",data_flag,"_err;"));
    fig=figure();hold on;
    x = [1:size(stay_lock_data,2)+size(leave_lock_data,2)];
    xticks(1:size(stay_lock_data,2)+size(leave_lock_data,2));
    xticklabels([1:size(stay_lock_data,2) size(leave_lock_data,2)-1:-1:0]);
    xlabel("Trials (Stay1 + i | Leave - i)");
    ylabel(ylab);
    xline(size(stay_lock_data,2)+0.5,'--k','LineWidth',2)
    displaynames = arrayfun(@(x)sprintf('Patch %d',x),1:size(plot_data,1)-1,'uni',0);
    displaynames = [displaynames {'Patch Average'}];
    p = [];
    for idx=1:num_patches
        p = [p plot(x,plot_data(idx,:),'Marker','+','LineStyle',':','LineWidth',1)];
    end
    p = [p plot(x,plot_data(end,:),'Marker','o','LineWidth',2.5,'Color','#77AC30');];
    errorbar(x,plot_data(end,:),err_data,'LineWidth',1,'Color','#A2142F');
    legend(p,displaynames);hold off;
end

% function [x,y,y_fit,ph,lh,pl,ll,trange] = fit_erp(epoch_data,init_time,fin_time)
%     init_idx = find(epoch_data.times==init_time);
%     fin_idx = find(epoch_data.times==fin_time);
%     trange = sprintf('_%+4d_%+4d',init_time,fin_time);
%     x = transpose(epoch_data.times(init_idx:fin_idx));
%     erp = transpose(mean(mean(epoch_data.data([2 3 14],:,:),3),1));
%     y = erp(init_idx:fin_idx);
%     [f,~] = fit(x,y,'fourier4');
%     y_fit = f(x);
%     [ph,lh,wh,sh] = findpeaks(y_fit,x);
%     [pl,ll,wl,sl] = findpeaks(-y_fit,x);
% end

% function T = gendf(subid,epoch_data,preflag,short_end,long_end)
%     init_idx = find(epoch_data.times==300); fin_idx = find(epoch_data.times==500);
%     channel = 'CZ';cz_idx = find(strcmp({epoch_data.chanlocs.labels},{channel}));
%     channel = 'C3';c3_idx = find(strcmp({epoch_data.chanlocs.labels},{channel}));
%     channel = 'C4';c4_idx = find(strcmp({epoch_data.chanlocs.labels},{channel}));
%     
%     subj_idx = repmat(subid,epoch_data.trials,1);
%     
%     if preflag
%         condition = repmat(1,epoch_data.trials,1);
%     else
%         condition = repmat(2,epoch_data.trials,1);
%     end
%     
%     travel_time = zeros(epoch_data.trials,1);
% %     response_stay_default = ones(epoch_data.trials,1);
%     response = zeros(epoch_data.trials,1);
% 
% %     rt = zeros(epoch_data.trials,1);
%     max_cz = zeros(epoch_data.trials,1);
%     max_c3 = zeros(epoch_data.trials,1);
%     max_c4 = zeros(epoch_data.trials,1);
% 
% %     erp = mean(mean(epoch_data.data([cz_idx,c3_idx,c4_idx],:,:),3),1);
% %     x = transpose(epoch_data.times(init_idx:fin_idx));
% %     y = transpose(erp(init_idx:fin_idx));
% %     [f,~] = fit(x,y,'fourier4');
% %     y_fit = f(x);
% %     [ph,lh,wh,sh] = findpeaks(y_fit,x);
% %     [pl,ll,wl,sl] = findpeaks(-y_fit,x);
% 
%     for i = 1:epoch_data.trials
%         N_epochevent_idx = find(strcmp(epoch_data.epoch(i).eventtype,'TRIGGER EVENT N'),1);
%         N_event_idx = epoch_data.epoch(i).event(N_epochevent_idx);
%         N_urevent_idx = cell2mat(epoch_data.epoch(i).eventurevent(N_epochevent_idx));
% 
%         if short_end<long_end
%             if N_event_idx>0 && N_event_idx<short_end
%                 travel_time(i,1) = 5;
%             else
%                 travel_time(i,1) = 20;
%             end
%         else
%             if N_event_idx>0 && N_event_idx<long_end
%                 travel_time(i,1) = 20;
%             else
%                 travel_time(i,1) = 5;
%             end
%         end
% 
%         if ismember(epoch_data.urevent(N_urevent_idx+1).type,{'TRIGGER EVENT O','TRIGGER EVENT H','condition 22','condition 23','condition 24'})
%             response(i,1) = 1;
% %         elseif strcmp(epoch_data.urevent(N_urevent_idx+1).type,'TRIGGER EVENT R')
% %             response_stay_default(i,1) = 0;
% %         else
% %             response_leave_default(i,1) = nan;
% %             response_stay_default(i,1) = nan;
%         end
% 
% %         if sum(strcmp(epoch_data.epoch(i).eventtype,'TRIGGER EVENT O'))
% %             response_leave_default(i,1) = 1;
% %         elseif sum(strcmp(epoch_data.epoch(i).eventtype,'TRIGGER EVENT R'))
% %             response_stay_default(i,1) = 0;
% %             R_epochevent_idx = find(strcmp(epoch_data.epoch(i).eventtype,'TRIGGER EVENT R'));
% %             R_urevent_idx = cell2mat(epoch_data.epoch(i).eventurevent(R_epochevent_idx));
% %             travel_time(i,1) = (epoch_data.urevent(R_urevent_idx+1).latency - epoch_data.urevent(R_urevent_idx).latency)/epoch_data.srate;
% %         end
% %         rt(i,1) = (epoch_data.urevent(N_urevent_idx+1).latency - epoch_data.urevent(N_urevent_idx).latency)/epoch_data.srate;
% %         rt(i,1) = cell2mat(epoch_data.epoch(i).eventlatency(idx_N+1))/epoch_data.srate;
% 
% %         [max_val,max_idx] = max(epoch_data.data(cz_idx,N0_idx:N500_idx,i)); max_idx = max_idx+N0_idx;
% %         max_cz(i) = mean(epoch_data.data(cz_idx,max_idx-30:max_idx+30,i));
%         max_cz(i,1) = mean(epoch_data.data(cz_idx,init_idx:fin_idx,i));
%     
% %         [max_val,max_idx] = max(epoch_data.data(c3_idx,N0_idx:N500_idx,i)); max_idx = max_idx+N0_idx;
% %         max_c3(i) = mean(epoch_data.data(c3_idx,max_idx-30:max_idx+30,i));
%         max_c3(i,1) = mean(epoch_data.data(c3_idx,init_idx:fin_idx,i));
% 
% %         [max_val,max_idx] = max(epoch_data.data(c4_idx,N0_idx:N500_idx,i)); max_idx = max_idx+N0_idx;
% %         max_c4(i) = mean(epoch_data.data(c4_idx,max_idx-30:max_idx+30,i));
%         max_c4(i,1) = mean(epoch_data.data(c4_idx,init_idx:fin_idx,i));
%     end
%     
% %     T = table(subj_idx,condition,travel_time,response_stay_default,response_leave_default,max_cz,max_c3,max_c4);
%     T = table(subj_idx,condition,travel_time,response,max_cz,max_c3,max_c4);
% 
% end