function [stat] = lfp_phasecoherence(data,whichdata)
%written by Seth Konig 11/8/15
%calculates phase coherence of LFP locked to a defined event

cfgfrq = [];
cfgfrq.output      = 'fourier';
cfgfrq.method      = 'mtmconvol';
cfgfrq.foi         = 4:2:30;
numfoi             = length(cfgfrq.foi);
cfgfrq.taper       = 'hanning';
cfgfrq.pad         = 'maxperlen';
cfgfrq.keeptrials  = 'yes';
cfgfrq.keeptapers  = 'yes';
cfgfrq.complex     = 'complex';
cfgfrq.t_ftimwin   = 0.2 * ones(1,numfoi);
cfgfrq.toi         = -0.75:.01:1.6;

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