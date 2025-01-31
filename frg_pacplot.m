clear;clc;
% eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
% wrk_dir = '/home/decision_lab/work/github/frg/';
% dir_sep = '/';

dir_sep = '\';
% eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1\";
eeglab_dir='E:\wrk\MATLAB Add-Ons\Collections\EEGLAB';
cd(eeglab_dir)
eeglab nogui;

wrk_dir = "C:\Users\promitmoitra\Documents\GitHub\frg\";
data_path = [char(wrk_dir),'data/'];
cd(data_path)
stai = readtable('frg_stai.csv');

% pac_dir = [data_path 'epoch/pac_data/'];
pac_dir = 'E:\wrk\Neuroflow\epoch\pac';
cd(pac_dir)

%%
% set_fnames = dir('./*.set');
% fname_split = split({set_fnames(:).name},'_');
% subids = unique(str2num(char(fname_split(:,:,1))));

anx_cats = {'high_anx','low_anx'};
state_conditions = {'preshort','prelong', 'postshort', 'postlong'};
trial_categories = {'early','mid','late','leave'};
chan1 = 'FZ'; chan2 = 'CZ'; band1 = 'theta'; band2 = 'gamma';

badids = [37532 38058 39862 43543 45528 47801 48278];
% subids = [31730,43000];
stai = sortrows(stai,"STAI_Trait");
high_anx_subids = stai{end-10:end,["ParticipantID"]};
high_anx_subids = setdiff(high_anx_subids,badids)

low_anx_subids = stai{1:10,["ParticipantID"]};
low_anx_subids = setdiff(low_anx_subids,badids)

%%
% sub_idx = 1; global subid; subid = subids(sub_idx);
% 
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

%%
mean_ermi = struct();

for sc_idx = 1:length(state_conditions)
    sc = state_conditions{sc_idx};
    mean_ermi.(sc) = struct();

    for cat_idx = 1:length(trial_categories)
        tcat = trial_categories{cat_idx};
        mean_ermi.(sc).(tcat) = struct();

        for anx_cat_idx = 1:length(anx_cats)
            anx_cat = char(anx_cats(anx_cat_idx));
            subids = eval([anx_cat '_subids;']);
            mean_ermi.(sc).(tcat).(anx_cat) = struct();
            mean_ermi.(sc).(tcat).(anx_cat).tpts = [];
            for sub_idx = 1:length(subids)
                subid = subids(sub_idx);
                fname = strjoin({num2str(subid),sc,tcat,chan1,chan2,band1,band2},'_');
                try
                    eeg_dat = pop_loadset('filename',[fname '.set'],'filepath',pac_dir);
                    for chan_pair_idx = 1:length([eeg_dat.etc.eegpac.labels])
                        chan_pair = char(eeg_dat.etc.eegpac(chan_pair_idx).labels);
                        chan_pair(strfind(chan_pair,'-')) = '_';
                        if ~sum(contains(fieldnames(mean_ermi.(sc).(tcat).(anx_cat)),chan_pair))
                            mean_ermi.(sc).(tcat).(anx_cat).(chan_pair) = struct();
                            mean_ermi.(sc).(tcat).(anx_cat).(chan_pair).ermi = [];
                        end
                        [out,ermi,tpts] = frg_plottrialpac(eeg_dat,...
                                     'plotindx',chan_pair_idx,...
                                     'freqval1',5.5,'freqval2',35,...
                                     'timerange',[-500 0]);
                         C1 = size(out,2);
                         mean_ermi.(sc).(tcat).(anx_cat).(chan_pair).ermi(:,end+1:end+C1) = out;
                    end
                catch
                    sprintf('%s not found',fname)
                    % pause(3)
                end

            end
            mean_ermi.(sc).(tcat).(anx_cat).tpts = tpts;
        end

    end

end

%%
% fig_upper = figure('OuterPosition',[-0.0158    0.3962    1.5536    0.4688]*1e3); hold on;
figure; hold on;
chan_pair = 'FZ_CZ';

% sc = 'preshort';
% sc = 'prelong';
sc = 'postshort';
% sc = 'postlong';

% anx_cat = 'high_anx';
anx_cat = 'low_anx';

% title([anx_cat,' ',sc])

x = mean_ermi.(sc).early.high_anx.tpts;

y1 = mean(mean_ermi.(sc).early.(anx_cat).(chan_pair).ermi,2);
y2 = mean(mean_ermi.(sc).mid.(anx_cat).(chan_pair).ermi  ,2);
y3 = mean(mean_ermi.(sc).late.(anx_cat).(chan_pair).ermi ,2);
y4 = mean(mean_ermi.(sc).leave.(anx_cat).(chan_pair).ermi,2);
% plot(x,y1,'DisplayName','early');
% plot(x,y2,'DisplayName','mid');
% plot(x,y3,'DisplayName','late');
% plot(x,y4,'DisplayName','leave');

y1_mean = mean(mean_ermi.(sc).early.(anx_cat).(chan_pair).ermi,'all');
y2_mean = mean(mean_ermi.(sc).mid.(anx_cat).(chan_pair).ermi  ,'all');
y3_mean = mean(mean_ermi.(sc).late.(anx_cat).(chan_pair).ermi ,'all');
y4_mean = mean(mean_ermi.(sc).leave.(anx_cat).(chan_pair).ermi,'all');
% plot(x,y1_mean*ones(size(x)),'DisplayName','early');
% plot(x,y2_mean*ones(size(x)),'DisplayName','mid');
% plot(x,y3_mean*ones(size(x)),'DisplayName','late');
% plot(x,y4_mean*ones(size(x)),'DisplayName','leave');

y1_sem = std(y1)/sqrt(length(y1));
y2_sem = std(y2)/sqrt(length(y2));
y3_sem = std(y3)/sqrt(length(y3));
y4_sem = std(y4)/sqrt(length(y4));
% errorbar(x,y1,y1_sem,'--k','HandleVisibility','off');
% errorbar(x,y2,y2_sem,'--k','HandleVisibility','off');
% errorbar(x,y3,y3_sem,'--k','HandleVisibility','off');
% errorbar(x,y4,y4_sem,'--k','HandleVisibility','off');
% legend;

anx_cat(strfind(anx_cat,'_'))='-';
errorbar(1:4,[y1_mean,y2_mean,y3_mean,y4_mean],[y1_sem,y2_sem,y3_sem,y4_sem],'DisplayName',[anx_cat,' ',sc])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anx_cat = 'high_anx';
% anx_cat = 'low_anx';

% fig_lower = figure('OuterPosition',[-0.0078    0.0362    1.5504    0.4480]*1e3); hold on;

% sc = 'preshort'; 
% sc = 'prelong'; 
% sc = 'postshort';
% sc = 'postlong'; 

y1 = mean(mean_ermi.(sc).early.(anx_cat).(chan_pair).ermi,2);
y2 = mean(mean_ermi.(sc).mid.(anx_cat).(chan_pair).ermi  ,2);
y3 = mean(mean_ermi.(sc).late.(anx_cat).(chan_pair).ermi ,2);
y4 = mean(mean_ermi.(sc).leave.(anx_cat).(chan_pair).ermi,2);
% plot(x,y1,'DisplayName','early');
% plot(x,y2,'DisplayName','mid');
% plot(x,y3,'DisplayName','late');
% plot(x,y4,'DisplayName','leave')

y1_mean = mean(mean_ermi.(sc).early.(anx_cat).(chan_pair).ermi,'all');
y2_mean = mean(mean_ermi.(sc).mid.(anx_cat).(chan_pair).ermi  ,'all');
y3_mean = mean(mean_ermi.(sc).late.(anx_cat).(chan_pair).ermi ,'all');
y4_mean = mean(mean_ermi.(sc).leave.(anx_cat).(chan_pair).ermi,'all');
% plot(x,y1_mean*ones(size(x)),'DisplayName','early');
% plot(x,y2_mean*ones(size(x)),'DisplayName','mid');
% plot(x,y3_mean*ones(size(x)),'DisplayName','late');
% plot(x,y4_mean*ones(size(x)),'DisplayName','leave');

y1_sem = std(y1)/sqrt(length(y1));
y2_sem = std(y2)/sqrt(length(y2));
y3_sem = std(y3)/sqrt(length(y3));
y4_sem = std(y4)/sqrt(length(y4));
% errorbar(x,y1,y1_sem,'--k','HandleVisibility','off');
% errorbar(x,y2,y2_sem,'--k','HandleVisibility','off');
% errorbar(x,y3,y3_sem,'--k','HandleVisibility','off');
% errorbar(x,y4,y4_sem,'--k','HandleVisibility','off');

anx_cat(strfind(anx_cat,'_'))='-';
errorbar(1:4,[y1_mean,y2_mean,y3_mean,y4_mean],[y1_sem,y2_sem,y3_sem,y4_sem],'DisplayName',[anx_cat,' ',sc])
legend;