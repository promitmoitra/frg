clear;clc;
eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
wrk_dir = '/home/decision_lab/work/github/frg/';
dir_sep = '/';

% eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1\";
% wrk_dir = "C:\Users\promitmoitra\Documents\GitHub\frg\";
% dir_sep = '\';

global data_path; data_path = [char(wrk_dir),'data/epoch/'];

pac_dir = [data_path '/pac/'];

cd(eeglab_dir)
eeglab nogui;
cd(data_path)

set_fnames = dir('./*.set');
fname_split = split({set_fnames(:).name},'_');
subids = unique(str2num(char(fname_split(:,:,1))));
% subids = [15089];
state_conditions = {'preshort','prelong', 'postshort', 'postlong'};
trial_categories = {'early','mid','late','leave'};
chan1 = 'FZ'; chan2 = 'T5';
theta = [4 8]; gamma = [30 90];
%%
for idx = 1:length(subids)
    global subid; subid = subids(idx);
    for sc_idx = 1:length(state_conditions)
        for cat_idx = 1:length(trial_categories)
            fname = [num2str(subid)  '_' state_conditions{sc_idx},'_' ...
                                trial_categories{cat_idx}];
            try
                fprintf([repmat('=',1,80),'\n','Loading epoched set:',fname,'\n',...
                        repmat('=',1,80),'\n'])
                ep_dat = pop_loadset('filename',[fname '.set'],...
                                    'filepath',data_path,'verbose','off');
                try
                    fprintf('Checking for pac set...\n')
                    pac_dat = pop_loadset('filename',[fname '_' chan1 '_'...
                                          chan2 '_' 'theta_gamma.set'],...
                                          'filepath',pac_dir,'verbose','off');
                    fprintf([fname ' processed. Loading next...\n'])
                    continue
                catch
                    fprintf([repmat('=',1,80),'\n','Initializing PAC calclulation for ',...
                            fname,'\n',repmat('=',1,80),'\n'])
                    chan1_idx = find(ismember({ep_dat.chanlocs.labels},{chan1}));
                    chan2_idx = find(ismember({ep_dat.chanlocs.labels},{chan2}));
                    ep_dat = pop_pac(ep_dat,'Channels',theta,gamma,...
                                    repelem([chan1_idx;chan2_idx],2,1)',repmat([chan1_idx;chan2_idx],2,1)',...
                                    'method','ermipac','nboot',200,'alpha',[],...
                                    'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
                    ep_dat = pop_saveset(ep_dat,'filename',[fname '_' chan1 '_'...
                                         chan2 '_' 'theta_gamma.set'],...
                                         'filepath',pac_dir,'savemode','onefile');
                end
            catch
                fprintf([repmat('=',1,80),'\n','No trials found for ',fname,'\n',...
                        repmat('=',1,80),'\n'])
                continue
            end
        end
    end
end
