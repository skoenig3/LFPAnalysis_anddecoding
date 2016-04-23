clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures\'; %where to save figures

cvtnew_files = {'PW140729_4-sorted.nex','PW140730_4-sorted.nex',...
    'PW140806_4-sorted.nex','PW140829_4-sorted.nex',...
    'PW140903_3-sorted.nex','PW140909_3-sorted.nex','PW140916_3-sorted.nex',...
    'PW141007_4-sorted.nex','PW141015_4-sorted.nex','PW141009_4-sorted.nex',...
    'PW141008_4-sorted.nex','PW141024_4-sorted.nex','PW141028_4-sorted.nex',...
    'PW141105_4-sorted.nex','PW150130_3-sorted.nex',...
   'PW150205_4-sorted.nex'};


%set/get some general important info
trial_start_code = 15;
dot_on_code = 25;
dot_clrchng_code = 27;
fixation_on_cross = 8;
bar_code_response = 4; %monkey made a move
reward_code = 3;
imageX = 800; %horizontal size of the screen
imageY = 600; %horizontal size of the screen

%time is + 300 ms
short = [700 1133]; %short duration trials
mediumm = [1134 1567];
long = [1568 2500];%cap should be 2000 but cortex can have some lag

all_short = cell(1,length(cvtnew_files));
all_medium = cell(1,length(cvtnew_files));
all_long = cell(1,length(cvtnew_files));

Fs = 1000;

Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity


%row 1 locked to fixation cross hair, row 2 locked to image onset, col by file
%all unfiltered/sorted trials
avg_lfp = cell(2,length(cvtnew_files));
all_high_powerspctrm = cell(2,3,length(cvtnew_files));
all_low_powerspctrm = cell(2,3,length(cvtnew_files));
all_coherence = cell(2,3,length(cvtnew_files));
all_coherence2High = cell(2,3,length(cvtnew_files));

for file = 1:length(cvtnew_files)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import & Preprocess Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load([data_dir cvtnew_files{file}(1:end-11) '-preprocessed'],'cfg','meta');
    
    %preallocate space and parallel structure of cfg
    successful_trials = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if sum(cfg.trl(t).allval == 3) > 0;
            successful_trials(t) = t;
        end
    end
    successful_trials = laundry(successful_trials);
    cfg.trl = cfg.trl(successful_trials);
    cfg.trlold = cfg.trl; %store since may need if and we're going to write over
    
    
    num_trials = length(successful_trials);
    all_event_times = NaN(length(successful_trials),2);
    for t = 1:num_trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        all_event_times(t,1) = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_on_cross)-trial_start;
        all_event_times(t,2) = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    buffer1 = 2048;
    buffer2 = 2048;
    
    %keep track of sequences and item number for later analysis and store
    %by event/eye movement
    dot_aligned = NaN(num_trials,3); %align to dot turning on
    fixation_aligned = NaN(num_trials,3); %aligned to fixation on cross
    index = 1;
    for trl = 1:size(all_event_times,1);
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        fixation_aligned(trl,:) = [all_event_times(trl,1)-1-buffer1+startind ...
            all_event_times(trl,1)-1+buffer2+startind -buffer1];
        dot_aligned(trl,:) = [all_event_times(trl,2)-1-buffer1+startind ...
            all_event_times(trl,2)-1+buffer2+startind -buffer1];
    end
    
    
    %clean up/remove any NaNs (occur when saccades couldn't be detected)
    dot_aligned = laundry(dot_aligned);
    fixation_aligned = laundry(fixation_aligned);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Find Short-Long Trials---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    short_trials = [];
    medium_trials = [];
    long_trials = [];
    
    durs = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        dot_duration = cfg.trl(t).alltim(cfg.trl(t).allval == 27)-cfg.trl(t).alltim(cfg.trl(t).allval == 25);
        durs(t) = dot_duration;
        if dot_duration < short(1)
            disp('Trial duration too short?')
        elseif dot_duration <= short(2)
            short_trials =  [short_trials t];
        elseif dot_duration <= mediumm(2)
            medium_trials =  [medium_trials t];
        elseif dot_duration <= long(2)
            long_trials =  [long_trials t];
        else
            disp('Trial duration too long?')
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import LFP data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %also filters out line noise
    
    %code to find spike/LFP channels
    LFPchannels = find_desired_channels(cfg,'LFP');
    
    % read the data from file and preprocess them
    cfg.channel       = cfg.channel(LFPchannels)';
    cfg.dftfilter     = 'yes';
    cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
    cfg.padding       = 1;
    
    %saccade aligned LFP data
    cfg.trl = dot_aligned;
    dot_aligneddata = ft_preprocessing(cfg);
    
    %fixation aligned LFP data
    
    cfg.trl = fixation_aligned;
    fixation_aligneddata = ft_preprocessing(cfg);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Time Locked Analysis---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear cfg
    cfg.channel       = 'all';
    cfg.covariance    = 'no';
    cfg.keeptrials    = 'yes';
    cfg.removemean    = 'yes';
    cfg.vartrllength  =  2;
    
    %fixation aligned LFP data
    timelock = ft_timelockanalysis(cfg,fixation_aligneddata);
    figure
    plot(timelock.time,timelock.avg')
    xlim([-0.75 1.0])
    line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
    xlabel('Time (sec)')
    ylabel('Voltage (uV)')
    title('fixation on cross Aligned')
    avg_lfp{1,file} = timelock;
    save_and_close_fig(figure_dir,['CVTNEW_fixationcross_LFP' cvtnew_files{file}(1:end-11)])

    
    %saccade aligned LFP data
    timelock = ft_timelockanalysis(cfg,dot_aligneddata);
    figure
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.7])
    line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
    xlabel('Time (sec)')
    ylabel('Voltage (uV)')
    title('Dot Aligned')
    avg_lfp{2,file} = timelock;
    save_and_close_fig(figure_dir,['CVTNEW_DotonLFP' cvtnew_files{file}(1:end-11)])

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Power Spectrum for Short, Medium and Long Trials---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %for low frequencies
    
    c = NaN(6,2);
    freq = cell(2,3);
    triallength = {short_trials,medium_trials,long_trials};
    min_length = [short(1) mediumm(1) long(1)]/1000;
    
    plotnum = [1 3 5; 2 4 6];
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    
    figure
    for event = 1:2;
        for d = 1:3
            if event == 1;
                eventdata = fixation_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = dot_aligneddata;
                title2 = 'Aligned to Dot on';
            end
            
            freq = lfp_powerspectrum(eventdata,'low',triallength{d},'none');
            
            all_low_powerspctrm{event,d,file} = freq;
            
            if d == 1
                title1 = 'Short: ';
            elseif d == 2
                title1 = 'Medium: ';
            else
                title1 = 'Long: ';
            end
            
            subplot(3,2,plotnum(event,d))
            ft_singleplotTFR(cfgfrq, freq);
            xlim([-0.3 min_length(d)])
            title([title1 title2]);
            c(plotnum(event,d),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for d = 1:3
            subplot(3,2,plotnum(event,d))
            caxis([minc maxc])
        end
    end
    save_and_close_fig(figure_dir,['CVTNEW_LowFrequencyPower' cvtnew_files{file}(1:end-11)])
    
    
    %for high frequencies
    
    c = NaN(6,2);
    freq = cell(2,3);
    triallength = {short_trials,medium_trials,long_trials};
    min_length = [short(1) mediumm(1) long(1)]/1000;
    
    plotnum = [1 3 5; 2 4 6];
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    figure
    for event = 1:2;
        for d = 1:3
            if event == 1;
                eventdata = fixation_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = dot_aligneddata;
                title2 = 'Aligned to Dot on';
            end
            
            freq{event,d} = lfp_powerspectrum(eventdata,'high',triallength{d},'none');
            all_high_powerspctrm{event,d,file} = freq{event,d};
            
            if d == 1
                title1 = 'Short: ';
            elseif d == 2
                title1 = 'Medium: ';
            else
                title1 = 'Long: ';
            end
            
            subplot(3,2,plotnum(event,d))
            ft_singleplotTFR(cfgfrq, freq{event,d});
            xlim([-0.75 min_length(d)])
            title([title1 title2]);
            c(plotnum(event,d),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for d = 1:3
            subplot(3,2,plotnum(event,d))
            caxis([minc maxc])
        end
    end
    save_and_close_fig(figure_dir,['CVTNEW_HighFrequencyPower' cvtnew_files{file}(1:end-11)])

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Phase Coherence---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %for low frequencies
    
    c = NaN(6,2);
    freq = cell(2,3);
    triallength = {short_trials,medium_trials,long_trials};
    min_length = [short(1) mediumm(1) long(1)]/1000;
    
    plotnum = [1 3 5; 2 4 6];
    
    figure
    for event = 1:2;
        for d = 1:3
            if event == 1;
                eventdata = fixation_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = dot_aligneddata;
                title2 = 'Aligned to Dot on';
            end
            
            stat = lfp_phasecoherence(eventdata,'all');
            
            all_coherence{event,d,file} = stat;
            
            if d == 1
                title1 = 'Short: ';
            elseif d == 2
                title1 = 'Medium: ';
            else
                title1 = 'Long: ';
            end
            
            subplot(3,2,plotnum(event,d))
            imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
            axis xy
            xlim([-0.3 min_length(d)])
            line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
            colorbar
            colormap jet
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)')
            c(plotnum(event,d),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for d = 1:3
            subplot(3,2,plotnum(event,d))
            caxis([minc maxc])
        end
    end
    save_and_close_fig(figure_dir,['CVTNEW_HighFrequencyCoherence' cvtnew_files{file}(1:end-11)])

    %for high frequencies
    
    c = NaN(6,2);
    freq = cell(2,3);
    triallength = {short_trials,medium_trials,long_trials};
    min_length = [short(1) mediumm(1) long(1)]/1000;
    
    plotnum = [1 3 5; 2 4 6];
    
    figure
    for event = 1:2;
        for d = 1:3
            if event == 1;
                eventdata = fixation_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = dot_aligneddata;
                title2 = 'Aligned to Dot on';
            end
            
            stat = lfp_phasecoherence2High(eventdata,'all');
            
            all_coherence2High{event,d,file} = stat;
            
            if d == 1
                title1 = 'Short: ';
            elseif d == 2
                title1 = 'Medium: ';
            else
                title1 = 'Long: ';
            end
            
            subplot(3,2,plotnum(event,d))
            imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
            axis xy
            xlim([-0.3 min_length(d)])
            line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
            colorbar
            colormap jet
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)')
            c(plotnum(event,d),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for d = 1:3
            subplot(3,2,plotnum(event,d))
            caxis([minc maxc])
        end
    end
    save_and_close_fig(figure_dir,['CVTNEW_HighFrequencyCoherence' cvtnew_files{file}(1:end-11)])

    
end
save('CVTnewLFPdata')