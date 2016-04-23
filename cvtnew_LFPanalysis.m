function  cvtnew_LFPanalysis(data_dir,preprocessed_data_file,figure_dir)
%function to look at LFP data during cvtnew task

%%%---import & reformat data so that spikes are locked to dot position---%%%
load([data_dir preprocessed_data_file],'data','cfg','multiunit','num_units','meta','hdr');

[LFPchans] = find_desired_channels(hdr,data,'LFP');


%set/get some general important info
num_trials = length(cfg.trl);
Fs = data(1).fsample; %should be 1000
ITIstart_code = 15;
dot_on_code = 25;
dot_clrchng_code = 27;
bar_code_response = 4; %monkey made a move
reward_code = 3;
imageX = 800; %horizontal size of the screen
imageY = 600; %horizontal size of the screen


Fline       = 60; %removes harmonics so 120, 180, etc.
LFPchannels = find_desired_channels(hdr,data,'LFP');
Fs = hdr.Fs;

all_LFP = NaN(length(cfg.trl),3000);

for t = 1:length(cfg.trl);
    if any(cfg.trl(t).allval == reward_code); %only take correct trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
        pathon =  cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %meta data starts at 1st 666 which appears before dot actually turns on in event array
        dot_clrchng = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start;
        responded = cfg.trl(t).alltim(cfg.trl(t).allval == bar_code_response)-trial_start;
        
        
        LFPdata = [];
        for channels = 1:4
            LFPdata = [LFPdata; data(LFPchans(channels)).values{t}(pathon:dot_clrchng)];
        end
        
        % taken directly from fieldtrip ft_preproc_dftfilter.m
        % determine the size of the data
        [Nchans, Nsamples] = size(LFPdata);
        
        % determine the largest integer number of line-noise cycles that fits in the data
        sel = 1:round(floor(Nsamples * Fline/Fs) * Fs/Fline);
        
        % temporarily remove mean to avoid leakage
        mdat = mean(LFPdata(:,sel),2);
        zeroedLFPdat  = LFPdata - mdat(:,ones(1,Nsamples));
        
        % fit a sin and cos to the signal and subtract them
        time  = (0:Nsamples-1)/Fs;
        tmp  = exp(j*2*pi*Fline*time);                    % complex sin and cos
        % ampl = 2*dat*tmp'/Nsamples;                  % estimated amplitude of complex sin and cos
        ampl = 2*zeroedLFPdat(:,sel)*tmp(sel)'/length(sel);     % estimated amplitude of complex sin and cos on integer number of cycles
        est  = ampl*tmp;                               % estimated signal at this frequency
        %filt = dat - est;                              % subtract estimated signal
        LFPdata = LFPdata - est; %do want to remove the mean
        LFPdata = real(LFPdata);
        LFPdata = LFPdata./(std(LFPdata')'*ones(1,length(LFPdata))); %normalize since impedance changes absolute power
        
        all_LFP(t,1:length(LFPdata)) = LFPdata(1,:);
    end
end

figure
plot(nanmean(all_LFP(:,1:1500)));
xlabel('time')

%%


%%