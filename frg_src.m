%% Load EEG and preproc...
EEG_rsbl = pop_rsbl(EEG,false,false,'mean','bsbl',1);

cng_idxs = find(~cellfun(@isempty,regexp(EEG_rsbl.etc.src.roi,"\w*anteriorcingulate\w*")));