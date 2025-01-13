clear;clc;
eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
wrk_dir = '/home/decision_lab/work/github/frg/';
dir_sep = '/';

% eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1\";
% wrk_dir = "C:\Users\promitmoitra\Documents\GitHub\frg\";
% dir_sep = '\';

global data_path; data_path = [char(wrk_dir),'data/epoch_data/'];

cd(eeglab_dir)
eeglab nogui;
cd(data_path)

set_fnames = dir('**/*.set');
fname_split = split({set_fnames(:).name},'_');
subids = unique(str2num(char(fname_split(:,:,1))));
state_conditions = {'preshort','prelong', 'postshort', 'postlong'};

%%
for idx = 1:length(subids)
    global subid; subid = subids(idx);
    for sc_idx = 1:length(state_conditions)
        early = pop_loadset('filename', [num2str(subid) '_' state_conditions{sc_idx} '_early.set'],...
                'filepath',data_path);
        try
            mid = pop_loadset('filename', [num2str(subid)  '_' state_conditions{sc_idx} '_mid.set'],...
                    'filepath',data_path);
        catch
            mid = struct();
            sprintf('no mid trials for %d',subid)
        end
        late = pop_loadset('filename', [num2str(subid) '_' state_conditions{sc_idx} '_late.set'],...
                'filepath',data_path);
        leave = pop_loadset('filename', [num2str(subid) '_' state_conditions{sc_idx} '_leave.set'],...
                'filepath',data_path);
%%
        pac_dir = [data_path '/pac_data/'];
        
        fz = 'FZ'; cz = 'CZ';
        fz_idx = find(ismember({early.chanlocs.labels},{fz}));
        cz_idx = find(ismember({early.chanlocs.labels},{cz}));
    
        tic
        if ~isempty(fieldnames(early))
        early = pop_pac(early,'Channels',[4 8],[30 90],[fz_idx fz_idx cz_idx],[fz_idx cz_idx cz_idx],...
                        'method','ermipac','nboot',200,'alpha',[],...
                        'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
        early = pop_saveset(early,'filename',[num2str(subid) '_' state_conditions{sc_idx} '_early_' fz '_' cz '_theta-gamma.set'],...
                            'filepath',pac_dir,'savemode','onefile');
        end

        if ~isempty(fieldnames(mid))
        mid = pop_pac(mid,'Channels',[4 8],[30 90],[fz_idx fz_idx cz_idx],[fz_idx cz_idx cz_idx],...
                        'method','ermipac','nboot',200,'alpha',[],...
                        'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
        mid = pop_saveset(mid,'filename',[num2str(subid) '_' state_conditions{sc_idx} '_mid_' fz '_' cz '_theta-gamma.set'],...
                            'filepath',pac_dir,'savemode','onefile');
        
        end

        if ~isempty(fieldnames(late))
        late = pop_pac(late,'Channels',[4 8],[30 90],[fz_idx fz_idx cz_idx],[fz_idx cz_idx cz_idx],...
                        'method','ermipac','nboot',200,'alpha',[],...
                        'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
        late = pop_saveset(late,'filename',[num2str(subid) '_' state_conditions{sc_idx} '_late_' fz '_' cz '_theta-gamma.set'],...
                            'filepath',pac_dir,'savemode','onefile');
        
        end

        if ~isempty(fieldnames(leave))
        leave = pop_pac(leave,'Channels',[4 8],[30 90],[fz_idx fz_idx cz_idx],[fz_idx cz_idx cz_idx],...
                        'method','ermipac','nboot',200,'alpha',[],...
                        'nfreqs1',10,'nfreqs2',20,'freqscale','log','bonfcorr',0);
        leave = pop_saveset(leave,'filename',[num2str(subid) '_' state_conditions{sc_idx} '_leave_' fz '_' cz '_theta-gamma.set'],...
                            'filepath',pac_dir,'savemode','onefile');
        end

        toc
    end
end    