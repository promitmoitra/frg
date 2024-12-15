function [times,frequencies,ersp,exponents,offsets,bp,ap_bp] = ap_exp_ts(eeg,chan_idx,band)
    [ersp,~,~,times,frequencies,~,~,~] = ...
    pop_newtimef(eeg,1,chan_idx,[eeg.xmin eeg.xmax]*1000,[1 0.6],...
                'freqs',[1 120],'baseline',NaN,'basenorm','off',...
                'plotersp','off','plotitc','off','scale','abs'); %'scale','log'

    ersp_abs = ersp;%10.^(ersp./10);
    lo_f = band(1); hi_f = band(2);
    half_win = 0.75; %%Sliding window!

    exponents = []; offsets = [];
    bp = []; ap_bp = [];
    for tpoint=times./1000
        x = frequencies(frequencies>lo_f & frequencies<hi_f);
        y = mean(ersp_abs((frequencies>lo_f & frequencies<hi_f),...
                        times>(tpoint-half_win)*1000 & times<(tpoint+half_win)*1000),2);
        fooof_msg = fprintf('FOOOFing... Timepoint:%.4f\n',tpoint);
        fooof_res = fooof(x,y,[x(1),x(end)],struct(),0);
        offsets(end+1) = fooof_res.aperiodic_params(1);
        exponents(end+1) = fooof_res.aperiodic_params(2);
        bp(:,end+1) = fooof_res.bandpowers;
        ap_bp(:,end+1) = fooof_res.ap_bandpowers;
        fprintf(repmat('\b', 1, (fooof_msg)));
    end
    sprintf('Done!')
end