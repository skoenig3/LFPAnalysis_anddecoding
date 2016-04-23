function [freq] = lfp_powerspectrum(data,lowhigh,whichtrials)
%calculate power spectrum for lfp data for 'low' or 'high' frequency

if ~strcmpi(whichtrials,'all')
    data.sampleinfo = data.sampleinfo(whichtrials,:);
    data.time = data.time(whichtrials);
    data.trial = data.trial(whichtrials);
    data.cfg.trl = data.cfg.trl(whichtrials,:);    
end

cfgfrq = [];

if strcmpi(lowhigh,'low')
    cfgfrq.output      = 'pow';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.taper       = 'hanning';
    cfgfrq.foi         = 2:1:30;
    cfgfrq.tapsmofrq   = 4;%4 for < 30 and 10 for > 30
    cfgfrq.t_ftimwin   = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = -0.75:.01:1.75;
    
    freq = ft_freqanalysis(cfgfrq,data);
    
elseif strcmpi(lowhigh,'high')
    cfgfrq.output      = 'pow';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.taper       = 'hanning';
    cfgfrq.foi         = 30:2:180;
    cfgfrq.tapsmofrq   = 10;%4 for < 30 and 10 for > 30
    cfgfrq.t_ftimwin   = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = -0.75:.01:1.75;
    
    freq = ft_freqanalysis(cfgfrq,data);
    
end

