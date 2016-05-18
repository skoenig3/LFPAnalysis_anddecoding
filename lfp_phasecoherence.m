function [stat] = lfp_phasecoherence(data,lowhigh,whichdata,time_window)
%written by Seth Konig 11/8/15
%calculates phase coherence of LFP locked to a defined event

if strcmpi(lowhigh,'low')
    
    cfgfrq = [];
    cfgfrq.output      = 'fourier';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.foi         = 4:1:30;
    cfgfrq.tapsmofrq   = 4;%4 for < 30 and 10 for > 30
    cfgfrq.taper       = 'hanning';
    cfgfrq.pad         = 'maxperlen';
    cfgfrq.keeptrials  = 'yes';
    cfgfrq.keeptapers  = 'yes';
    cfgfrq.complex     = 'complex';
    cfgfrq.t_ftimwin   = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = time_window;

    freq = ft_freqanalysis(cfgfrq, data);
    
    sizfrq=size(freq.label,1);
    freq.label{sizfrq+1} = 'fake';
    siz = size(freq.fourierspctrm);
    freq.fourierspctrm(:,sizfrq+1,:,:) = complex(ones(siz(1),1,siz(3),siz(4)), ...
        zeros(siz(1),1,siz(3),siz(4)));
    
    % create channelcmb
    lfpind=strmatch('AD',freq.label);
    cfgfrq.channelcmb=num2cell(zeros(length(lfpind),2)+nan);
    for k=1:length(lfpind)
        cfgfrq.channelcmb{k,1} = freq.label{lfpind(k)};
        cfgfrq.channelcmb{k,2} = 'fake';
    end
    
    % determine the descriptive spectral statistics
    cfgfrq.method      = 'coh';
    stat = ft_connectivityanalysis(cfgfrq, freq);
    
elseif strcmpi(lowhigh,'high')
    
    cfgfrq = [];
    cfgfrq.output      = 'fourier';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.foi         = 30:2:120;
    cfgfrq.taper       = 'hanning';
    cfgfrq.pad         = 'maxperlen';
    cfgfrq.keeptrials  = 'yes';
    cfgfrq.keeptapers  = 'yes';
    cfgfrq.complex     = 'complex';
    cfgfrq.t_ftimwin   = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = time_window;
    
    if ~strcmpi(whichdata,'all')
        data.sampleinfo = data.sampleinfo(whichdata,:);
        data.time = data.time(whichdata);
        data.trial = data.trial(whichdata);
        data.cfg.trl = data.cfg.trl(whichdata,:);
    end
    
    freq = ft_freqanalysis(cfgfrq, data);
    
    sizfrq=size(freq.label,1);
    freq.label{sizfrq+1} = 'fake';
    siz = size(freq.fourierspctrm);
    freq.fourierspctrm(:,sizfrq+1,:,:) = complex(ones(siz(1),1,siz(3),siz(4)), ...
        zeros(siz(1),1,siz(3),siz(4)));
    
    % create channelcmb
    lfpind=strmatch('AD',freq.label);
    cfgfrq.channelcmb=num2cell(zeros(length(lfpind),2)+nan);
    for k=1:length(lfpind)
        cfgfrq.channelcmb{k,1} = freq.label{lfpind(k)};
        cfgfrq.channelcmb{k,2} = 'fake';
    end
    
    % determine the descriptive spectral statistics
    cfgfrq.method      = 'coh';
    stat = ft_connectivityanalysis(cfgfrq, freq);
end