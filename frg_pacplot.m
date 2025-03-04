clear;clc;

dir_sep = '/';
eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
wrk_dir = '/home/decision_lab/work/github/frg/';

% dir_sep = '\';
% % eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1\";
% % eeglab_dir='E:\wrk\MATLAB Add-Ons\Collections\EEGLAB';
% wrk_dir = "C:\Users\promitmoitra\Documents\GitHub\frg\";
% % pac_dir = 'E:\wrk\Neuroflow\epoch\pac';

% cd(eeglab_dir); eeglab nogui;

data_path = [char(wrk_dir),'data/'];
epoch_dir = [data_path 'epoch/'];

cd(wrk_dir)
stai = readtable('frg_stai.csv');
anx_score = table2array(stai(:,["STAI_Trait"]));
trait_anx_cat = discretize(anx_score,[min(anx_score):12:max(anx_score)+12]);
stai.("trait_anx_cat") = trait_anx_cat;

anx_cats = {'low','mid','high'};
stress_conds = {'pre','post'};
tt_envs = {'short','long'};
trial_categories = {'early','mid','late','leave'};

badids = [37532 38058 39862 43543 45528 47801 48278];

% chan1 = 'FZ'; chan2 = {'CZ','T5'}; band1 = 'theta'; band2 = 'gamma';
% chan_pair = [chan1,'_',chan2];

%% TODO: FIND MODULAT(OR/ED) FREQS USING COMODT, PEAK FREQS USING FOOOF AND COMPARE
%% TODO: CHANNEL PAIR SELECTION?
%% SANDBOX:
% subids = [31730,43000];
% num_sub = 5;
% stai = sortrows(stai,"STAI_Trait");
% high_anx_subids = stai{end-(num_sub+1):end,["ParticipantID"]};
% high_anx_subids = setdiff(high_anx_subids,badids)
% 
% low_anx_subids = stai{1:num_sub,["ParticipantID"]};
% low_anx_subids = setdiff(low_anx_subids,badids)
%
% sub_idx = 1; global subid; subid = subids(sub_idx);
% sc_idx = 1;
% cat_idx = 1;
% % chan_pair = [chan1 '-' chan1];
% chan_pair = [chan1 '-' chan2];
% % chan_pair = [chan2 '-' chan2];
% 
% fname = strjoin({num2str(subids(sub_idx)),state_conditions{sc_idx},...
%                 trial_categories{cat_idx},chan1,chan2,band1,band2},'_');
% 
% eeg_dat = pop_loadset('filename',[fname '.set'],'filepath',pac_path);
% chan_pair_idx = find(contains([eeg_dat.etc.eegpac.labels],chan_pair));
% phase_f = eeg_dat.etc.eegpac(chan_pair_idx).params.freqs_phase;
% amp_f = eeg_dat.etc.eegpac(chan_pair_idx).params.freqs_amp;
% pacdat = eeg_dat.etc.eegpac(chan_pair_idx).ermipac.pacval;
% time_points = eeg_dat.etc.eegpac(chan_pair_idx).ermipac.times;
% % comodulogram(phase_f,amp_f,pacdat); hold on;cbar;
% % comodulogramt(phase_f,amp_f,time_points(time_points<0),pacdat,'npoints',11,'scale','log'); hold on;cbar;
% %%%Yes, it's a nested function call :')
% % pop_plotpac('plot_trialpac',...
% %             eeg_plotcomodt(eeg_dat,'plotindx',2,'timerange',[-500,0]));%,...
% %                             'freqval1',5.5,'freqval2',30))
% [outdata,erp] = frg_plottrialpac(eeg_dat,'plotindx',chan_pair_idx,...
%                             'timerange',[-500 0],'freqval1',5.5,'freqval2',30);
%% Extract ERMI - epoch windows: precue [-500 0], postcue [0 500] | freqs: 5.5 x 35
% cd(epoch_dir)
% set_fnames = dir('./*.set');
% fname_split = split({set_fnames(:).name},'_');
% subids = unique(str2num(char(fname_split(:,:,1))));
%
% mean_ermi = struct();
% for anx_cat_idx = 1:length(anx_cats)
%     anx_cat = anx_cats{anx_cat_idx};
%     subids = table2array(stai(stai.trait_anx_cat==anx_cat_idx,["ParticipantID"]));
%     subids = setdiff(subids,badids);
%     mean_ermi.(anx_cat) = struct();
%     for sc_idx = 1:length(stress_conds)
%         sc = stress_conds{sc_idx};
%         mean_ermi.(anx_cat).(sc) = struct();
%         for env_idx = 1:length(tt_envs)
%             env = tt_envs{env_idx};
%             mean_ermi.(anx_cat).(sc).(env) = struct();
%             mean_ermi.(anx_cat).(sc).(env).tpts = [];
%             sctt = [sc,env];
%             for cat_idx = 1:length(trial_categories)
%                 tcat = trial_categories{cat_idx};
%                 mean_ermi.(anx_cat).(sc).(env).(tcat) = struct();
% 
%                 for sub_idx = 1:length(subids)
%                     subid = subids(sub_idx);
%                     fname = strjoin({num2str(subid),sctt,tcat},'_');
% %                     fname = strjoin({num2str(subid),sctt,tcat,chan1,chan2,band1,band2},'_');
%                     try
%                         eeg_dat = pop_loadset('filename',[fname '.set'],'filepath',epoch_dir);
%                         for chan_pair_idx = 1:length([eeg_dat.etc.eegpac.labels])
%                             chan_pair = char(eeg_dat.etc.eegpac(chan_pair_idx).labels);
%                             chan_pair(strfind(chan_pair,'-')) = '_';
%                             if ~sum(contains(fieldnames(mean_ermi.(anx_cat).(sc).(env).(tcat)),chan_pair))
%                                 mean_ermi.(anx_cat).(sc).(env).(tcat).(chan_pair) = struct();
%                                 mean_ermi.(anx_cat).(sc).(env).(tcat).(chan_pair).ermi = [];
%                             end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             [out,ermi,tpts] = frg_plottrialpac(eeg_dat,...
%                                          'plotindx',chan_pair_idx,...
%                                          'freqval1',5.5,'freqval2',35,...
%                                          'timerange',[0 500]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              C1 = size(out,2);
%                              mean_ermi.(anx_cat).(sc).(env).(tcat).(chan_pair).ermi(:,end+1:end+C1) = out;
%                         end
%                     catch
%                         sprintf('%s not found',fname)
%                     end
%                 end
%             end
%         end
%         mean_ermi.(anx_cat).(sc).(env).tpts = tpts;
%     end
% end
% cd(wrk_dir)
% save('ermi_postcue.mat','mean_ermi')
%% Plotting:
% cd(wrk_dir);
% load('ermi_precue.mat')
% chan_pair = 'FZ_CZ';
% %%
% % anx_cats = {'low','mid'};
% % anx_cats = {'mid','high'};
% anx_cats = {'low','high'};
% stress_conds = {'pre','post'};
% tt_envs = {'short','long'};
% 
% % var_lvl1=anx_cats;
% % var_lvl2=stress_conds;
% % tt = 'short';
% % % tt = 'long';
% % var_fixed=tt;
% % var_order=[1,2,3];
% 
% var_lvl1=stress_conds;
% var_lvl2=tt_envs;
% anx_cat = 'low';
% % anx_cat = 'mid';
% % anx_cat = 'high';
% var_fixed=anx_cat;
% var_order=[3,1,2];
% 
% % var_lvl1=tt_envs;
% % var_lvl2=anx_cats;
% % sc = 'pre';
% % % sc = 'post';
% % var_fixed=sc;
% % var_order=[2,3,1];
% 
% figure('OuterPosition',[246 131 922 923]);hold on;
% lwid=2;msiz=8;
% hdl=[];lbl={};
% clrs=['b','r'];mrks=['o','s'];
% 
% for var_idx1 = 1:length(var_lvl1)
%     var_name1 = var_lvl1{var_idx1};
%     for var_idx2 = 1:length(var_lvl2)
%         var_name2 = var_lvl2{var_idx2};
%         var_n = {var_name1,var_name2,var_fixed};
%         lbl{end+1} = strjoin(var_n(var_order),'-');
%         for tc = 1:length(trial_categories)
%             tcat = trial_categories{tc};
%             y = eval(['mean_ermi.' strjoin(var_n(var_order),'.') '.(tcat).(chan_pair).ermi']);
%             y_t_mean = mean(y,1);
%             y_t_sub_mean = mean(y,'all');
%             y_err = std(y_t_mean)/sqrt(length(y_t_mean));
%             clr=clrs(var_idx1);mrk=mrks(var_idx2);
%             h = errorbar(tc,y_t_sub_mean,y_err,[clr mrk],'LineWidth',lwid,'MarkerSize',msiz);
%             xticks(1:4);xticklabels(trial_categories);xlim([0.5 4.5])
%         end
%         hdl(end+1) = h;
%     end
% end
% legend(hdl,lbl)
% %%
% % anx_cats = {'low','mid'};
% % anx_cats = {'mid','high'};
% % anx_cats = {'low','high'};
% anx_cats = {'low','mid','high'};
% stress_conds = {'pre','post'};
% tt_envs = {'short','long'};
% 
% % var_lvl1=anx_cats;
% % sc = 'pre';
% % % sc = 'post';
% % var_fixed1=sc;
% % tt = 'short';
% % % tt = 'long';
% % var_fixed2=tt;
% % var_order=[1,2,3];
% 
% % var_lvl1=stress_conds;
% % tt = 'short';
% % % tt = 'long';
% % var_fixed1=tt;
% % anx_cat = 'low';
% % % anx_cat = 'mid';
% % % anx_cat = 'high';
% % var_fixed2=anx_cat;
% % var_order=[3,1,2];
% 
% % var_lvl1=tt_envs;
% % anx_cat = 'low';
% % % anx_cat = 'mid';
% % % anx_cat = 'high';
% % var_fixed1=anx_cat;
% % sc = 'pre';
% % % sc = 'post';
% % var_fixed2=sc;
% % var_order=[2,3,1];
% 
% figure('OuterPosition',[246 131 922 923]);hold on;
% lwid=2;msiz=8;
% hdl=[];lbl={};
% clrs=['b','g','r'];mrks=['o','s','d'];
% 
% for var_idx1 = 1:length(var_lvl1)
%     var_name1 = var_lvl1{var_idx1};
%     var_n = {var_name1,var_fixed1,var_fixed2};
%     lbl{end+1} = strjoin(var_n(var_order),'-');
%     for tc = 1:length(trial_categories)
%         tcat = trial_categories{tc};
%         y = eval(['mean_ermi.' strjoin(var_n(var_order),'.') '.(tcat).(chan_pair).ermi']);
%         y_t_mean = mean(y,1);
%         y_t_sub_mean = mean(y,'all');
%         y_err = std(y_t_mean)/sqrt(length(y_t_mean));
%         clr=clrs(var_idx1);mrk=mrks(var_idx1);
%         h = errorbar(tc,y_t_sub_mean,y_err,[clr mrk],'LineWidth',lwid,'MarkerSize',msiz);
%         xticks(1:4);xticklabels(trial_categories);xlim([0.5 4.5])
%     end
%     hdl(end+1) = h;
% end
% legend(hdl,lbl)
% %%
% cd(wrk_dir);
% load('ermi_precue.mat')
% 
% var_n={};
% 
% anx_cat='low';sc='pre';tt='short';
% var_n{1} = {anx_cat,sc,tt};
% anx_cat='mid';sc='pre';tt='long';
% var_n{2} = {anx_cat,sc,tt};
% anx_cat='high';sc='pre';tt='long';
% var_n{3} = {anx_cat,sc,tt};
% 
% % anx_cat='low';sc='pre';tt='long';
% % var_n{1} = {anx_cat,sc,tt};
% % anx_cat='mid';sc='pre';tt='short';
% % var_n{2} = {anx_cat,sc,tt};
% % anx_cat='high';sc='post';tt='long';
% % var_n{3} = {anx_cat,sc,tt};
% 
% % anx_cat='mid';sc='pre';tt='long';
% % var_n{1} = {anx_cat,sc,tt};
% % anx_cat='mid';sc='post';tt='short';
% % var_n{2} = {anx_cat,sc,tt};
% % anx_cat='mid';sc='post';tt='long';
% % var_n{3} = {anx_cat,sc,tt};
% 
% num_plot=length(var_n);
% lwid=2;msiz=8;
% hdl=[];lbl={};
% clrs = ['b','g','r'];
% for idx = 1:num_plot
%     if idx==1
%         figure('OuterPosition',[246 131 922 923]);hold on;
%     end
%     lbl{end+1} = strjoin(var_n{idx},'-');
%     for tc = 1:length(trial_categories)
%         tcat = trial_categories{tc};
%         y = eval(['mean_ermi.' strjoin(var_n{idx},'.') '.(tcat).(chan_pair).ermi']);
%         y_t_mean = mean(y,1);
%         y_t_sub_mean = mean(y,'all');
%         y_err = std(y_t_mean)/sqrt(length(y_t_mean));
%         clr=clrs(idx);mrk='o';
%         h = errorbar(tc,y_t_sub_mean,y_err,[clr mrk],'LineWidth',lwid,'MarkerSize',msiz);
%         xticks(1:4);xticklabels(trial_categories);xlim([0.5 4.5])
%     end
%     hdl(end+1) = h;
% end
% legend(hdl,lbl)
%% CHECKING DIFFERENT CHANNELS:
cd(wrk_dir);
load('ermi_postcue.mat')

a = {'FZ','CZ','T5'};
[A,B] = meshgrid(a,a);
c = string(A)+'_'+string(B);
chan_pairs = char(c(:));
chan_pairs = chan_pairs([1:5 7 9],:);

var_n = chan_pairs;

% for cp_idx =1:length(chan_pairs)
%     cp = chan_pairs(cp_idx);
% end    

anx_cat='high';sc='pre';tt='short';

num_plot=length(var_n);
lwid=2;msiz=8;
hdl=[];lbl={};

for idx = 1:num_plot
    if idx==1
        figure('OuterPosition',[246 131 922 923]);hold on;
    end
    chan_pairs(idx,strfind(chan_pairs(idx,:),'_')) = '-';
    lbl{end+1} = chan_pairs(idx,:);%strjoin(var_n{idx},'-');
    y_mean = [];
    for tc = 1:length(trial_categories)
        tcat = trial_categories{tc};
        y = mean_ermi.(anx_cat).(sc).(tt).(tcat).(var_n(idx,:)).ermi;
        y_t_sub_mean = mean(y,'all');
        y_mean(end+1) = y_t_sub_mean;        
%         y_t_mean = mean(y,1);
%         y_err = std(y_t_mean)/sqrt(length(y_t_mean));
%         h = errorbar(tc,y_t_sub_mean,y_err,'LineWidth',lwid,'MarkerSize',msiz);
    end
    h = plot(1:4,y_mean,'-o','LineWidth',lwid,'MarkerSize',msiz);
    xticks(1:4);xticklabels(trial_categories);xlim([0.5 4.5])    
    hdl(end+1) = h;
end
legend(hdl,lbl)
%%
%%% Discuss what to plot: fluctuation around temporal mean (y_sub_err) 
%%% and across subject variability at each timepoint (y_t_err)
% %     title('Effect of anxiety and travel time - pre/post stress')
% %     labels = {'low anx short','high anx short','low anx long','high anx long'};
% %     title('Effect of anxiety and stress - short/long')
% %     labels = {'low anx pre','high anx pre','low anx post','high anx post'};
% %     title('Effect of stress and travel time - low/high')
% %     labels = {'pre short','pre long','post short','post long'};