%% Foraging EEG analysis:
% cd 'C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1';
clear;clc;eeglab;close;
% cd 'C:\Users\promitmoitra\Documents\GitHub\frg\'

% data_path = '/home/decision_lab/work/frg/foraging/Neuroflow/';
% badids = [37532 38058 39862 42125 43543 45194 46037 47678 47744 47801 48238 48278];
% %%% badids are [37532 38058 39862 42125 43543 '45194' 46037 '47678' 47744 47801 48238 48278];
% subids = readmatrix(fullfile(data_path,'subids.txt')); subids = setdiff(subids,badids);
% %%% sub 47801 has no poststress long travel time markers TRIGGER EVENT U
% 
global subid;
subid = 7873;%subids(subids==7873);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean line noise and bandpass filter:
preEEG = clean_data_with_zapline_plus_eeglab_wrapper(preEEG,struct('noisefreqs',60,'chunkLength',0,'adaptiveNremove',true, ...
                                                    'fixedNremove',1,'plotResults',0));
postEEG = clean_data_with_zapline_plus_eeglab_wrapper(postEEG,struct('noisefreqs',60,'chunkLength',0,'adaptiveNremove',true, ...
                                                    'fixedNremove',1,'plotResults',0));

global freqs
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

% global preshort_epoch; global prelong_epoch; global postshort_epoch; global postlong_epoch;

preshort_epoch = pop_epoch(preshortEEG,lock_event,epoch_trange);
prelong_epoch = pop_epoch(prelongEEG,lock_event,epoch_trange);
postshort_epoch = pop_epoch(postshortEEG,lock_event,epoch_trange);
postlong_epoch = pop_epoch(postlongEEG,lock_event,epoch_trange);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SpecParam
global channel; global chan_idx;
% channel = 'CZ';cz_idx = find(strcmp({EEG.chanlocs.labels},{channel}));chan_idx = cz_idx;
channel = 'C3';c3_idx = find(strcmp({EEG.chanlocs.labels},{channel}));chan_idx = c3_idx;
% channel = 'C4';c4_idx = find(strcmp({EEG.chanlocs.labels},{channel}));chan_idx = c4_idx;

big_plot("slopes")
big_plot("beta_bp")

% data_flag = "slopes";
% big_fig_exp = figure('Name',strcat(num2str(subid)," ",channel));
% 
% stress_flag='pre'; tt_flag='short';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag);%,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax1_exp = copyobj(plot_ax,big_fig_exp);
% 
% stress_flag='pre'; tt_flag='long';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag);%,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax2_exp = copyobj(plot_ax,big_fig_exp);
% 
% stress_flag='post'; tt_flag='short';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag);%,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax3_exp = copyobj(plot_ax,big_fig_exp);
% 
% stress_flag='post'; tt_flag='long';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag);%,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax4_exp = copyobj(plot_ax,big_fig_exp);
% 
% linkaxes([ax1_exp ax2_exp ax3_exp ax4_exp])
% figure(big_fig_exp);
% subplot(2,2,1,ax1_exp);
% subplot(2,2,2,ax2_exp);
% subplot(2,2,3,ax3_exp);
% subplot(2,2,4,ax4_exp);
% 
% han=axes(big_fig_exp,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,"Aperiodic exponent");
% xlabel(han,"Trials (Stay1 + i | Leave - i)");
% title(han,strcat(num2str(subid)," ",channel),'Position',[0.5 1.05 0.5]);
% saveas(gcf,strcat("./fig/",num2str(subid),"_",channel,"_slopes.png"));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% data_flag = "beta_bp";
% big_fig_beta = figure('Name',strcat(num2str(subid)," ",channel),'Visible','off');
% 
% stress_flag='pre'; tt_flag='short';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax1_beta = copyobj(plot_ax,big_fig_beta);
% 
% stress_flag='pre'; tt_flag='long';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax2_beta = copyobj(plot_ax,big_fig_beta);
% 
% stress_flag='post'; tt_flag='short';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax3_beta = copyobj(plot_ax,big_fig_beta);
% 
% stress_flag='post'; tt_flag='long';
% eval(strcat("epoch_data = ",stress_flag,tt_flag,"_epoch;"));
% [stay_lock_res,leave_lock_res,fig] = run_plot(epoch_data,data_flag,chan_idx,freqs);
% title(strcat(stress_flag," ",tt_flag));
% plot_ax = fig.get('CurrentAxes');
% ax4_beta = copyobj(plot_ax,big_fig_beta);
% 
% %save stay_lock_res,leave_lock_res
% 
% linkaxes([ax1_beta ax2_beta ax3_beta ax4_beta])
% figure(big_fig_beta);
% subplot(2,2,1,ax1_beta);
% subplot(2,2,2,ax2_beta);
% subplot(2,2,3,ax3_beta);
% subplot(2,2,4,ax4_beta);
% 
% han=axes(big_fig_beta,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,"Beta band power (\muV^2)",'Position',[-0.1 0.5 0]);
% xlabel(han,"Trials (Stay1 + i | Leave - i)");
% title(han,strcat(num2str(subid)," ",channel),'Position',[0.5 1.05 0.5]);
% saveas(gcf,strcat("./fig/",num2str(subid),"_",channel,"_beta.png"));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [stay_lock_res, leave_lock_res] = specparam(ep_dat)
    global chan_idx; global freqs;
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
            stay_lock_slopes(end,idx)   = fooof_prestim.aperiodic_params(end);
            stay_lock_delta_bp(end,idx) = mean(fooof_prestim.bandpowers(1),'omitnan');
            stay_lock_theta_bp(end,idx) = mean(fooof_prestim.bandpowers(2),'omitnan');
            stay_lock_alpha_bp(end,idx) = mean(fooof_prestim.bandpowers(3),'omitnan');
            stay_lock_beta_bp(end,idx)  = mean(fooof_prestim.bandpowers(4),'omitnan');
            stay_lock_gamma_bp(end,idx) = mean(fooof_prestim.bandpowers(5),'omitnan');
        end
        eval(strcat("leave_",num2str(max(patch_trial_len)-idx),"= pop_select(ep_dat,'trial',leave_lock_trial_idxs(~isnan(leave_lock_trial_idxs(:,", ...
                    num2str(idx),")),",num2str(idx),"));"));
        eval(strcat("leave_",num2str(max(patch_trial_len)-idx),".xmin=-1;leave_",num2str(max(patch_trial_len)-idx),".xmax=1.9960;"))
        eval(strcat("leave_lock_fooof_prestim = eeg_fooof(leave_",num2str(max(patch_trial_len)-idx),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"));
    %     eval(strcat("leave_lock_fooof_prestim = eeg_fooof(leave_",num2str(max(patch_trial_len)-idx),",'channel',[1:ep_dat.nbchan],[-1 0]*1000,100,freqs,fooof_settings);"));
        fooof_prestim = cell2mat(leave_lock_fooof_prestim.etc.FOOOF_results(chan_idx));
    %     leave_lock_slopes(1,idx) = fooof_prestim(chan_idx).aperiodic_params(end);
        leave_lock_slopes(end,end-idx+1)   = fooof_prestim.aperiodic_params(end);
        leave_lock_delta_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(1),'omitnan');
        leave_lock_theta_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(2),'omitnan');
        leave_lock_alpha_bp(end,end-idx+1) = mean(fooof_prestim.bandpowers(3),'omitnan');
        leave_lock_beta_bp(end,end-idx+1)  = mean(fooof_prestim.bandpowers(4),'omitnan');
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
                stay_lock_slopes(patch_idx,idx)   = fooof_prestim.aperiodic_params(end);
                stay_lock_delta_bp(patch_idx,idx) = mean(fooof_prestim.bandpowers(1),'omitnan');
                stay_lock_theta_bp(patch_idx,idx) = mean(fooof_prestim.bandpowers(2),'omitnan');
                stay_lock_alpha_bp(patch_idx,idx) = mean(fooof_prestim.bandpowers(3),'omitnan');
                stay_lock_beta_bp(patch_idx,idx)  = mean(fooof_prestim.bandpowers(4),'omitnan');
                stay_lock_gamma_bp(patch_idx,idx) = mean(fooof_prestim.bandpowers(5),'omitnan');
            end
            eval(strcat("leave_",num2str(patch_idx),"_",num2str(idx-1),"= pop_select(patch,'trial',leave_trial_idx-idx+1);"));
            eval(strcat("leave_",num2str(patch_idx),"_",num2str(idx-1),".xmin=-1;leave_",num2str(patch_idx),"_",num2str(idx-1),".xmax=1.9660;"));
            eval(strcat("leave_fooof_prestim = eeg_fooof(leave_",num2str(patch_idx),"_",num2str(idx-1),",'channel',",num2str(chan_idx),",[-1 0]*1000,100,freqs,fooof_settings);"))
            fooof_prestim = cell2mat(leave_fooof_prestim.etc.FOOOF_results(chan_idx));
            leave_lock_slopes(patch_idx,end-idx+1)   = fooof_prestim.aperiodic_params(end);
            leave_lock_delta_bp(patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(1),'omitnan');
            leave_lock_theta_bp(patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(2),'omitnan');
            leave_lock_alpha_bp(patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(3),'omitnan');
            leave_lock_beta_bp(patch_idx,end-idx+1)  = mean(fooof_prestim.bandpowers(4),'omitnan');
            leave_lock_gamma_bp(patch_idx,end-idx+1) = mean(fooof_prestim.bandpowers(5),'omitnan');        
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

function fig = frg_plot(num_patches,data_flag,plot_data,stay_lock_data,leave_lock_data)
    eval(strcat("[",data_flag,"_err,",data_flag,"]=std(plot_data(1:end-1,:),1,1,'omitnan');"));
    eval(strcat("err_data=",data_flag,"_err;"));
    fig=figure('Visible','off'); hold on; ax=gca;
    x = [1:size(stay_lock_data,2)+size(leave_lock_data,2)];
    xticks(ax,x);
    xticklabels(ax,[1:size(stay_lock_data,2) size(leave_lock_data,2)-1:-1:0]);
    xline(ax,size(stay_lock_data,2)+0.5,'--k','LineWidth',2)
%     displaynames = arrayfun(@(x)sprintf('Patch %d',x),1:size(plot_data,1)-1,'uni',0);
%     displaynames = [displaynames {'Patch Average'}];
    p = [];
    for idx=1:num_patches
        p = [p plot(ax,x,plot_data(idx,:),'Marker','+','LineStyle',':','LineWidth',1)];
    end
    p = [p plot(ax,x,plot_data(end,:),'Marker','o','LineWidth',2.5,'Color','#77AC30');];
    errorbar(ax,x,plot_data(end,:),err_data,'LineWidth',1,'Color','#A2142F');
%     legend(p,displaynames);
    hold off;
end

function [stay_lock_res,leave_lock_res,fig] = run_plot(ep_dat,data_flag)
    global chan_idx; global freqs;
    leave_trial_idxs = get_leave_idx(ep_dat);
    first_stay_idxs = [1 leave_trial_idxs(1:end-1)+1];
    patch_trial_idxs = arrayfun(@(f,g) (f:g),first_stay_idxs,leave_trial_idxs,'UniformOutput',false);
    num_patches = size(patch_trial_idxs,2);
    
    [stay_lock_res, leave_lock_res] = specparam(ep_dat);%,chan_idx,freqs);
    
    eval(strcat("stay_lock_data = stay_lock_res.",data_flag,";"));
    % frac_nan=sum(isnan(stay_lock_data),1)/(size(stay_lock_data,1)-1);
    % stay_lock_data=stay_lock_data(:,frac_nan<=0.25);
    stay_lock_data=stay_lock_data(:,1:int8(size(stay_lock_data,2)/2));
    
    eval(strcat("leave_lock_data = leave_lock_res.",data_flag,";"));
    % frac_nan=sum(isnan(leave_lock_data),1)/(size(leave_lock_data,1)-1);
    % leave_lock_data=leave_lock_data(:,frac_nan<=0.25);
    leave_lock_data=leave_lock_data(:,int8(size(leave_lock_data,2)-3):end);
    
    plot_data = [stay_lock_data leave_lock_data];
    fig = frg_plot(num_patches,data_flag,plot_data,stay_lock_data,leave_lock_data);
end    

function big_plot(data_flag)
    global subid; global channel;
    disp(subid)
    global preshort_epoch; global prelong_epoch; global postshort_epoch; global postlong_epoch;
    big_fig = figure('Name',strcat(num2str(subid)," ",channel));
    
    stress_flag='pre'; tt_flag='short';
    eval(strcat("ep_dat = ",stress_flag,tt_flag,"_epoch;"));
%     ep_dat = preshort_epoch;
    [stay_lock_res,leave_lock_res,fig] = run_plot(ep_dat,data_flag);
    title(strcat(stress_flag," ",tt_flag));
    plot_ax = fig.get('CurrentAxes');
    ax1 = copyobj(plot_ax,big_fig);
    
    stress_flag='pre'; tt_flag='long';
    eval(strcat("ep_dat = ",stress_flag,tt_flag,"_epoch;"));
%     ep_dat = prelong_epoch;
    [stay_lock_res,leave_lock_res,fig] = run_plot(ep_dat,data_flag);
    title(strcat(stress_flag," ",tt_flag));
    plot_ax = fig.get('CurrentAxes');
    ax2 = copyobj(plot_ax,big_fig);
    
    stress_flag='post'; tt_flag='short';
    eval(strcat("ep_dat = ",stress_flag,tt_flag,"_epoch;"));
%     ep_dat = postshort_epoch;
    [stay_lock_res,leave_lock_res,fig] = run_plot(ep_dat,data_flag);
    title(strcat(stress_flag," ",tt_flag));
    plot_ax = fig.get('CurrentAxes');
    ax3 = copyobj(plot_ax,big_fig);
    
    stress_flag='post'; tt_flag='long';
    eval(strcat("ep_dat = ",stress_flag,tt_flag,"_epoch;"));
%     ep_dat = postlong_epoch;
    [stay_lock_res,leave_lock_res,fig] = run_plot(ep_dat,data_flag);
    title(strcat(stress_flag," ",tt_flag));
    plot_ax = fig.get('CurrentAxes');
    ax4 = copyobj(plot_ax,big_fig);
    
    linkaxes([ax1 ax2 ax3 ax4])
    figure(big_fig);
    subplot(2,2,1,ax1);
    subplot(2,2,2,ax2);
    subplot(2,2,3,ax3);
    subplot(2,2,4,ax4);
    
    han=axes(big_fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    title(han,strcat(num2str(subid)," ",channel),'Position',[0.5 1.05 0.5]);
    xlabel(han,"Trials (Stay1 + i | Leave - i)");    
    if data_flag == "slopes"
        ylabel(han,"Aperiodic exponent",'Position',[-0.1 0.5 0]);
        fname = strcat("./fig/",num2str(subid),"_",channel,"_slopes.png");
    elseif data_flag == "beta_bp"
        ylabel(han,"Beta band power (\muV^2)",'Position',[-0.1 0.5 0]);
        fname = strcat("./fig/",num2str(subid),"_",channel,"_beta_bp.png");        
    end
    saveas(gcf,fname);
end