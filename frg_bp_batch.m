function frg_bp_batch(chan)
    clc;
    rmpath_pat('NBTpublic')
    eeglab_dir = '/home/decision_lab/MATLAB Add-Ons/Collections/EEGLAB/';
    % eeglab_dir = "C:\Users\promitmoitra\Documents\MATLAB\eeglab2023.1";
    wrk_dir = '/home/decision_lab/work/github/frg';
    
    cd(eeglab_dir)
    eeglab nogui;
    cd(wrk_dir)
    
    data_path = [wrk_dir,'/data/'];
    badids = [37532 38058 39862 43543 45528 47801 48278];

%     chans = {'F8','FP1','FP2'};%{'FZ','FP1','FP2','F3','F4','F7','F8','FCz'};    
    
%     for chan=chans
        chan=char(chan);
        all_subids = readmatrix(fullfile(data_path,'subids.txt'));
        processed_ids = dir([data_path,'spec_data/',chan,'/*.mat']);
        proc_fnames = char({processed_ids(:).name});        
        for idx = 1:size(proc_fnames,1)
            sep_loc = strfind(proc_fnames(idx,:),'_');
            proc_id = proc_fnames(idx,sep_loc(1)+1:sep_loc(2)-1);
%             proc_id = char(regexp(proc_fnames(idx,:),'[0-9]','match'))'
            badids(end+1) = str2num(proc_id);
        end
        all_subids = setdiff(all_subids,badids);
%         batches = 1:4:length(all_subids)+1;
        subids = all_subids;%(batches(batch_idx):batches(batch_idx+1)-1);
        for sub_idx = 1:length(subids)
            clearvars -except sub_idx subids all_subids chan data_path
            subid = subids(sub_idx);
%             fprintf([repmat('=',1,80),'\n','Batch:',char(string(batch_idx)),'\n',char(string(subid)),'\n',repmat('=',1,80)])
            data_file = dir(fullfile(data_path,string(subid),'*.edf'));
            data_file_path = fullfile(data_file(end).folder,data_file(end).name);
            origchanlocs = readlocs([data_path,'Statnet_F3F4FCz.ced']);
            
            EEG = pop_biosig(data_file_path,'channels',1:19);
%             EEG = pop_select(EEG,'channel',1:19);
            EEG = pop_select(EEG,'rmchannel',{'EKG2'});

            EEG.data(end+1,:) = 0;EEG.nbchan = size(EEG.data,1);
            EEG = pop_chanedit(EEG,'load','./data/Statnet_F3F4FCz.ced');
            EEG = pop_reref(EEG,{'A1','A2'});
            EEG = pop_interp(EEG,origchanlocs);
            EEG = pop_reref(EEG,[]);
    
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
    
            band = [1 100];
    %       chans = {'F8','FP1','FP2'};%{'FZ','FP1','FP2','F3','F4','F7','F8','FCz'};
    %       chan = cell2mat(chan); 
            chan_idx = find(ismember({EEG.chanlocs.labels},chan));
            fprintf([repmat('=',1,80),'\n',EEG.chanlocs(chan_idx).labels,'\n',repmat('=',1,80),'\n'])
    
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
    
            state_conditions = {'pre','pre_short','pre_long','post_short','post_long'};
            spec_vars = {'times','freqs','ersp','exponents','offsets'};
            freq_bands = {'delta','theta','alpha','lbeta','hbeta','lgamma','hgamma'};
            
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
                end
            end
        
            spec_chan_dir = [data_path,'spec_data/',EEG.chanlocs(chan_idx).labels];
            if ~exist(spec_chan_dir,'dir')
                mkdir(spec_chan_dir)
            end
            save(sprintf('%s/spec_data/%s/spec_%d_%s',data_path,chan,subid,chan),'spec_data')
        end
%     end
cd ~
end