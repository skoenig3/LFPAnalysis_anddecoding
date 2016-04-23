clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Figures\';

listsq_files = {'PW140725_3-sorted.nex','PW140729_3-sorted.nex',...
    'PW140730_3-sorted.nex','PW140801_3-sorted.nex','PW140805_3-sorted.nex','PW140806_3-sorted.nex',...
    'PW140825_3-sorted.nex','PW140826_3-sorted.nex','PW140829_3-sorted.nex','PW140908_3-sorted.nex',...
    'PW140910_3-sorted.nex','PW140915_3-sorted.nex','PW140101_3-sorted.nex',...
    'PW141006_3-sorted.nex','PW141007_3-sorted.nex','PW141013_3-sorted.nex','PW141014_3-sorted.nex',...
    'PW141009_3-sorted.nex','PW141015_3-sorted.nex','PW141008_3-sorted.nex','PW141010_3-sorted.nex',...
    'PW141024_3-sorted.nex','PW141027_3-sorted.nex','PW141028_3-sorted.nex','PW141029_3-sorted.nex',...
    'PW141023_3-sorted.nex','PW141031_3-sorted.nex','PW141103_3-sorted.nex','PW141105_3-sorted.nex',...
    'PW141106_3-sorted.nex','PW141110_3-sorted.nex','PW150121_3-sorted.nex','PW150122_3-sorted.nex',...
    'PW150123_3-sorted.nex','PW150126_3-sorted.nex','PW150127_3-sorted.nex','PW150204_3-sorted.nex',...
    'PW150205_3-sorted.nex'};


% Important parametrs/values
twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
imgon_code = 23;
imgoff_code = 24;
fixspoton = 35;
fixation_oncross = 8;
trial_start_code = 15;
Fs = 1000;
fixrange = [4 15]; %first fixation to last fixation to look at for behavioral changes

Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity

%row 1 locked to fixation cross hair, row 2 locked to image onset, col by file
%all unfiltered/sorted trials
avg_lfp = cell(2,length(listsq_files));
all_high_powerspctrm = cell(2,length(listsq_files));
all_low_powerspctrm = cell(2,length(listsq_files));
all_fixdurs = cell(4,length(listsq_files)); %novel low novel high, rep low rep high
all_percent_change = cell(1,length(listsq_files));
all_coherence = cell(2,length(listsq_files));

%same thing but divided by "memory" strenght
all_low_recog_coherence = cell(2,length(listsq_files));
all_high_recog_coherence = cell(2,length(listsq_files));
all_recogdiff_coherence = cell(2,length(listsq_files));
all_low_recog_coherence2High = cell(2,length(listsq_files));
all_high_recog_coherence2High = cell(2,length(listsq_files));
all_recogdiff_coherence2High = cell(2,length(listsq_files));
all_lowrecog_low_power = cell(2,length(listsq_files));
all_highrecog_low_power = cell(2,length(listsq_files));
all_recogdiff_low_power = cell(2,length(listsq_files));
all_lowrecog_high_power = cell(2,length(listsq_files));
all_highrecog_high_power = cell(2,length(listsq_files));
all_recogdiff_high_power = cell(2,length(listsq_files));

for file = 1:length(listsq_files)
    disp(['File #' num2str(file) ': ' listsq_files{file}(1:end-11) ])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import & Preprocess Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([data_dir,listsq_files{file}(1:end-11)  '-preprocessed'],'cfg','fixationstats','item_set');
    
    %get important task specific information
    [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
    [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
    
    
    %preallocate space and parallel structure of cfg
    successful_trials = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if any(cfg.trl(t).allval == imgon_code) %in which image was displayed or 1st item in sequence was displayed
            if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)% sequence trials and only want rewarded ones
                successful_trials(t) = t;
            end
        end
    end
    successful_trials = laundry(successful_trials);
    fixationstats = fixationstats(successful_trials);
    cfg.trl = cfg.trl(successful_trials);
    
    num_trials = length(successful_trials);
    all_event_times = NaN(num_trials,2);
    event_codes = [fixation_oncross imgon_code];
    for t = 1:num_trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        for event = 1:length(event_codes);
            event_time = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(event));
            event_time = event_time(1);%fixation code can be found multiple times just want the first one
            all_event_times(t,event) = event_time-trial_start;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---process eye data and fixation durations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %only want trials in which we have a novel and repeat pair
    trials_with_img_pairs = NaN(1,num_trials);
    for t = 1:num_trials
        if sum(which_img == which_img(t)) == 2;
            trials_with_img_pairs(t) = t;
        end
    end
    trials_with_img_pairs = laundry(trials_with_img_pairs);
    
    if length(trials_with_img_pairs)/2 < 48
        continue %need more data than this
    end
    
    fixationstats = fixationstats(trials_with_img_pairs);
    cfg.trl = cfg.trl(trials_with_img_pairs);
    all_event_times = all_event_times(trials_with_img_pairs,:);
    which_img = which_img(trials_with_img_pairs);
    novel_vs_repeat = novel_vs_repeat(trials_with_img_pairs);
    num_trials = length(which_img);
    
    novel_fixdurs = NaN(num_trials/2,75);
    repeat_fixdurs = NaN(num_trials/2,75);
    for t = 1:length(fixationstats);
        if any(cfg.trl(t).allval == imgon_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code);%to ignore 1st fixation
            fixationtimes = fixationstats{t}.fixationtimes;
            if ~isempty(fixationtimes)
                tooearly = find(fixationtimes(1,:) <= imgon-trial_start);
                fixationtimes(:,tooearly) = [];
                fixdurs = diff(fixationtimes)+1;
                if novel_vs_repeat(t) == 1 %then novel
                    novel_fixdurs(which_img(t),1:length(fixdurs))=fixdurs;
                else
                    repeat_fixdurs(which_img(t),1:length(fixdurs))=fixdurs;
                end
            end
        end
    end
    novel_fixdurs = laundry(novel_fixdurs);
    repeat_fixdurs = laundry(repeat_fixdurs);
    
    all_fixdurs{1,file} = nanmean(novel_fixdurs(:,1:18));
    all_fixdurs{2,file} = nanmean(repeat_fixdurs(:,1:18));
    
    
    percent_change = NaN(1,num_trials/2);
    for trial = 1:num_trials/2
        percent_change(trial) = 100*(nanmean(repeat_fixdurs(trial,fixrange(1):fixrange(2)))-...
            nanmean(novel_fixdurs(trial,fixrange(1):fixrange(2))))/...
            nanmean(novel_fixdurs(trial,fixrange(1):fixrange(2)));
    end
    all_percent_change{file} = percent_change;
    
    low = find(percent_change <= prctile(percent_change,33));
    high = find(percent_change >= prctile(percent_change,67));
    
    all_fixdurs{1,file} = nanmean(novel_fixdurs(low,1:18));
    all_fixdurs{2,file} = nanmean(novel_fixdurs(high,1:18));
    all_fixdurs{3,file} = nanmean(repeat_fixdurs(low,1:18));
    all_fixdurs{4,file} = nanmean(repeat_fixdurs(high,1:18));
    
    figure
    subplot(1,3,1)
    hold on
    plot(nanmean(novel_fixdurs(:,1:18)))
    plot(nanmean(repeat_fixdurs(:,1:18)),'r')
    hold off
    legend('Novel','Repeat')
    xlabel('Fixation #')
    ylabel('Fixation duration (ms)')
    title('Average Fixation Duration')
    
    subplot(1,3,2)
    hist(percent_change,25)
    xlabel('% Change in Fixation Durations')
    
    subplot(1,3,3)
    hold on
    plot(nanmean(novel_fixdurs(low,1:18)),'b')
    plot(nanmean(repeat_fixdurs(low,1:18)),'r')
    plot(nanmean(novel_fixdurs(high,1:18)),'g')
    plot(nanmean(repeat_fixdurs(high,1:18)),'m')
    hold off
    legend('Novel (low)','Repeat (low)','Novel (high)','Repeat (high)')
    xlabel('Fixation #')
    ylabel('Fixation duration (ms)')
    title('Average Fixation Duration')
    save_and_close_fig(figure_dir,['ListSQ_List_FixationDurations_' listsq_files{file}(1:end-11)])
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    buffer1 = 2048;
    buffer2 = 2048;
    
    
    fixon_aligned = NaN(num_trials/2,3);
    img_aligned = NaN(num_trials/2,3);
    index = 1; 
    for t = 1:length(cfg.trl);
        if novel_vs_repeat(t) == 1
            startind = cfg.trl(t).begsmpind;%trial start index in all data
            
            fixon_aligned(index,:) = [all_event_times(t,1)-1-buffer1+startind ...
                all_event_times(t,2)-1+buffer2+startind -buffer1];
            img_aligned(index,:) = [all_event_times(t,2)-1-buffer1+startind ...
                all_event_times(t,2)-1+buffer2+startind -buffer1];
            
            index = index+1;
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
    
    %LFP data aligned to fixation on crosshair
    cfg.trl = fixon_aligned;
    fixon_aligneddata = ft_preprocessing(cfg);
    
    %LFP data aligned to image onset
    cfg.trl = img_aligned;
    img_aligneddata = ft_preprocessing(cfg);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Time Locked Analysis---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear cfg
    cfg.channel       = 'all';
    cfg.covariance    = 'no';
    cfg.keeptrials    = 'yes';
    cfg.removemean    = 'yes';
    cfg.vartrllength  =  2;
    
    %LFP data aligned to fixation on crosshair
    timelock = timelocklfp(cfg,fixon_aligneddata,figure_dir,['AverageLFP_fixationcross_' listsq_files{file}(1:end-11)]);
    avg_lfp{1,file} = timelock.avg';
    timelockfix = timelock.time;
    
    %LFP data aligned to image onset
    timelock = timelocklfp(cfg,img_aligneddata,figure_dir,['AverageLFP_ImageOnset_' listsq_files{file}(1:end-11)]);
    avg_lfp{2,file} = timelock.avg';
    timelockimage = timelock.time;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Low Frequency Power Spectrum---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = NaN(2,2);
    
    figure
    
    freq = lfp_powerspectrum(fixon_aligneddata,'low','all',[-1.25 -1.0]);
    all_low_powerspctrm{1,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,1)
    ft_singleplotTFR(cfgfrq, freq);
    title('Fix on Cross Aligned')
    c(1,:) = caxis;
    
    
    freq = lfp_powerspectrum(img_aligneddata,'low','all',[-1.25 -1.0]);
    all_low_powerspctrm{2,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,2)
    ft_singleplotTFR(cfgfrq, freq);
    title('Image Aligned')
    c(2,:) = caxis;
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    subtitle('Low Frequency Power Spectrum')
    save_and_close_fig(figure_dir,['ListSQ_List_LowFreqPowerSpectrum_' listsq_files{file}(1:end-11)])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate High Frequency Power Spectrum---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(2,2);
    
    figure
    
    freq = lfp_powerspectrum(fixon_aligneddata,'high','all',[-0.5 -0.25]);
    all_high_powerspctrm{1,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,1)
    ft_singleplotTFR(cfgfrq, freq);
    title('Fix on Cross Aligned')
    c(1,:) = caxis;
    
    
    freq = lfp_powerspectrum(img_aligneddata,'high','all',[-1.25 -1.0]);
    all_high_powerspctrm{2,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,2)
    ft_singleplotTFR(cfgfrq, freq);
    title('Image Aligned')
    c(2,:) = caxis;
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    subtitle('High Frequency Power Spectrum')
    save_and_close_fig(figure_dir,['ListSQ_List_HighFreqPowerSpectrum_' listsq_files{file}(1:end-11)])
    clear freq;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Phase Coherence---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(2,2);
    figure
    
    stat = lfp_phasecoherence(fixon_aligneddata,'all');
    all_coherence{1,file} = stat;
    statlow = stat;
    
    subplot(1,2,1)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    c(1,:) = caxis;
    
    stat = lfp_phasecoherence(img_aligneddata,'all');
    all_coherence{2,file} = stat;
    
    subplot(1,2,2)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    c(2,:) = caxis;
    
    subtitle('Phase Coherence')
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    save_and_close_fig(figure_dir,['ListSQ_List_PhaseCoherence_' listsq_files{file}(1:end-11)])
    clear stat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---High vs low Recognition Memory Low Freq Power Spectrum---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(6,2);
    freq = cell(2,2);
    lowhightrials = {low,high};
    plotnum = [1 2; 4 5];
    baseline = {[-0.5 -0.25],[-1.25 -1.0]};
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    
    figure
    for event = 1:2;
        for lh = 1:2
            if event == 1;
                eventdata = fixon_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = img_aligneddata;
                title2 = 'Aligned to Image Onset';
            end
            
            freq{event,lh} = lfp_powerspectrum(eventdata,'low',lowhightrials{lh},baseline{event});
            
            if lh == 1
                title1 = 'Low Recognition: ';
                all_lowrecog_low_power{event,file} = freq{event,lh};
            else
                title1 = 'High Recognition: ';
                all_highrecog_low_power{event,file} = freq{event,lh};
            end
            
            subplot(2,3,plotnum(event,lh))
            ft_singleplotTFR(cfgfrq, freq{event,lh});
            title([title1 title2]);
            c(plotnum(event,lh),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for lh = 1:2
            subplot(2,3,plotnum(event,lh))
            caxis([minc maxc])
        end
    end
    
    %calculate and plot difference btwn low and high recognition trials
    plotnum = [3 6];
    for event = 1:2
        freqd = freq{event,1};
        freqd.powspctrm = freq{event,2}.powspctrm-freq{event,1}.powspctrm;
        
        all_recogdiff_low_power{event,file} = freqd;
        
        subplot(2,3,plotnum(event))
        ft_singleplotTFR(cfgfrq, freqd);
        if event == 1;
            title('Aligned to Fix on Cross');
        else
            title('Aligned to Image Onset');
        end
    end
    subtitle('Low Frequency Power Spectrum High vs Low Recognition Memory')
    save_and_close_fig(figure_dir,['ListSQ_List_Recogntion_LowFrqPower_' listsq_files{file}(1:end-11)])
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---High vs low Recognition Memory High Freq Power Spectrum---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    c = NaN(6,2);
    freq = cell(2,2);
    lowhightrials = {low,high};
    plotnum = [1 2; 4 5];
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    figure
    for event = 1:2;
        for lh = 1:2
            if event == 1;
                eventdata = fixon_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = img_aligneddata;
                title2 = 'Aligned to Image Onset';
            end
            
            freq{event,lh} = lfp_powerspectrum(eventdata,'high',lowhightrials{lh},baseline{event});
            
            if lh == 1
                title1 = 'Low Recognition: ';
                all_lowrecog_high_power{event,file} = freq{event,lh};
            else
                title1 = 'High Recognition: ';
                all_highrecog_high_power{event,file} = freq{event,lh};
            end
            
            subplot(2,3,plotnum(event,lh))
            ft_singleplotTFR(cfgfrq, freq{event,lh});
            title([title1 title2]);
            c(plotnum(event,lh),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for lh = 1:2
            subplot(2,3,plotnum(event,lh))
            caxis([minc maxc])
        end
    end
    
    %calculate and plot difference btwn low and high recognition trials
    plotnum = [3 6];
    for event = 1:2
        freqd = freq{event,1};
        freqd.powspctrm = freq{event,2}.powspctrm-freq{event,1}.powspctrm;
        
        all_recogdiff_high_power{event,file} = freqd;
        
        subplot(2,3,plotnum(event))
        ft_singleplotTFR(cfgfrq, freqd);
        if event == 1;
            title('Aligned to Fix on Cross');
        else
            title('Aligned to Image Onset');
        end
    end
    subtitle('High Frequency Power Spectrum High vs Low Recognition Memory')
    save_and_close_fig(figure_dir,['ListSQ_List_Recogntion_HighFrqPower_' listsq_files{file}(1:end-11)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Phase Coherence for High vs Low Recognition---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(6,2);
    stat = cell(2,2);
    lowhightrials = {low,high};
    plotnum = [1 2; 4 5];
    
    figure
    for event = 1:2;
        for lh = 1:2
            if event == 1;
                eventdata = fixon_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = img_aligneddata;
                title2 = 'Aligned to Image Onset';
            end
            
            stat{event,lh} = lfp_phasecoherence(eventdata,lowhightrials{lh});
            
            if lh == 1
                title1 = 'Low Recognition: ';
                all_low_recog_coherence{event,file} = stat{event,lh};
            else
                title1 = 'High Recognition: ';
                all_high_recog_coherence{event,file} =  stat{event,lh};
            end
            
            subplot(2,3,plotnum(event,lh))
            imagesc(stat{event,lh}.time,stat{event,lh}.freq,squeeze(abs(stat{event,lh}.cohspctrm(1,:,:))))
            axis xy
            xlim([-0.3 0.5])
            line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
            colorbar
            colormap jet
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)')
            title([title1 title2]);
            c(plotnum(event,lh),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for lh = 1:2
            subplot(2,3,plotnum(event,lh))
            caxis([minc maxc])
        end
    end
    
    %calculate and plot difference btwn low and high recognition trials
    plotnum = [3 6];
    for event = 1:2
        statd = stat{event,1};
        statd.cohspctrm = stat{event,2}.cohspctrm-stat{event,1}.cohspctrm;
        
        all_recogdiff_coherence{event,file} = statd;
        
        subplot(2,3,plotnum(event))
        imagesc(statd.time,statd.freq,squeeze(abs(statd.cohspctrm(1,:,:))))
        axis xy
        xlim([-0.3 0.5])
        line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
        colorbar
        colormap jet
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)')
        if event == 1;
            title('Aligned to Fix on Cross');
        else
            title('Aligned to Image Onset');
        end
    end
    subtitle('Phase Coherence High vs Low Recognition Memory')
    save_and_close_fig(figure_dir,['ListSQ_List_Recogntion_PhaseCoh_' listsq_files{file}(1:end-11)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate High Frequency Phase Coherence for High vs Low Recognition---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(6,2);
    stat = cell(2,2);
    lowhightrials = {low,high};
    plotnum = [1 2; 4 5];
    
    figure
    for event = 1:2;
        for lh = 1:2
            if event == 1;
                eventdata = fixon_aligneddata;
                title2 = 'Aligned to Fix on Cross';
            else
                eventdata = img_aligneddata;
                title2 = 'Aligned to Image Onset';
            end
            
            stat{event,lh} = lfp_phasecoherence2High(eventdata,lowhightrials{lh});
            stathigh = stat{event,lh};
            
            if lh == 1
                title1 = 'Low Recognition: ';
                all_low_recog_coherence2High{event,file} = stat{event,lh};
            else
                title1 = 'High Recognition: ';
                all_high_recog_coherence2High{event,file} =  stat{event,lh};
            end
            
            subplot(2,3,plotnum(event,lh))
            imagesc(stat{event,lh}.time,stat{event,lh}.freq,squeeze(abs(stat{event,lh}.cohspctrm(1,:,:))))
            axis xy
            xlim([-0.3 0.5])
            line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
            colorbar
            colormap jet
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)')
            title([title1 title2]);
            c(plotnum(event,lh),:) = caxis;
        end
    end
    clear eventdata
    
    %rescale plots
    minc = min(c(:,1));
    maxc = max(c(:,2));
    for event = 1:2
        for lh = 1:2
            subplot(2,3,plotnum(event,lh))
            caxis([minc maxc])
        end
    end
    
    %calculate and plot difference btwn low and high recognition trials
    plotnum = [3 6];
    for event = 1:2
        statd = stat{event,1};
        statd.cohspctrm = stat{event,2}.cohspctrm-stat{event,1}.cohspctrm;
        
        all_recogdiff_coherence2High{event,file} = statd;
        
        subplot(2,3,plotnum(event))
        imagesc(statd.time,statd.freq,squeeze(abs(statd.cohspctrm(1,:,:))))
        axis xy
        xlim([-0.3 0.5])
        line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
        colorbar
        colormap jet
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)')
        if event == 1;
            title('Aligned to Fix on Cross');
        else
            title('Aligned to Image Onset');
        end
    end
    subtitle('Phase Coherence High vs Low Recognition Memory')
    save_and_close_fig(figure_dir,['ListSQ_List_Recogntion_PhaseCoh_' listsq_files{file}(1:end-11)])
end

save('All_ListSQ_List_LFPdata')
emailme('LFP analysis for List of ListSQ is done')


%% Plot Summary Data across all sessions


%for all data
lfps = cell(1,2);
high_powerspctrm = cell(1,2);
low_powerspctrm = cell(1,2);
percent_change = [];
coherence = cell(1,2);

%for fixation durations
novlow = [];
novhigh = [];
replow = [];
rephigh = [];

%same thing but divided by "memory" strength
low_recog_coherence = cell(1,2);
high_recog_coherence = cell(1,2);
recogdiff_coherence = cell(1,2);
lowrecog_low_power =cell(1,2);
highrecog_low_power = cell(1,2);
recogdiff_low_power = cell(1,2);
lowrecog_high_power = cell(1,2);
highrecog_high_power = cell(1,2);
recogdiff_high_power = cell(1,2);



col = 2;
eventsize = [4858 4097];
for event = 1:2
    
    lfps{event} = zeros(eventsize(event),4);
    high_powerspctrm{event} = zeros(size(all_high_powerspctrm{event,col}.powspctrm));
    low_powerspctrm{event} = zeros(size(all_low_powerspctrm{event,col}.powspctrm));
    coherence{event} = zeros(size(all_coherence{event,col}.cohspctrm));
    
    %same thing but divided by "memory" strength
    low_recog_coherence{event}  = zeros(size(all_low_recog_coherence{event,col}.cohspctrm));
    high_recog_coherence{event}  =zeros(size(all_high_recog_coherence{event,col}.cohspctrm));
    recogdiff_coherence{event}  = zeros(size(all_recogdiff_coherence{event,col}.cohspctrm));
    
    low_recog_coherence2High{event}  = zeros(size(all_low_recog_coherence2High{event,col}.cohspctrm));
    high_recog_coherence2High{event}  =zeros(size(all_high_recog_coherence2High{event,col}.cohspctrm));
    recogdiff_coherence2High{event}  = zeros(size(all_recogdiff_coherence2High{event,col}.cohspctrm));
    
    lowrecog_low_power{event}  = zeros(size(all_lowrecog_low_power{event,col}.powspctrm));
    highrecog_low_power{event}  = zeros(size(all_highrecog_low_power{event,col}.powspctrm));
    recogdiff_low_power{event}  = zeros(size(all_recogdiff_low_power{event,col}.powspctrm));
    lowrecog_high_power{event}  = zeros(size(all_lowrecog_high_power{event,col}.powspctrm));
    highrecog_high_power{event}  = zeros(size(all_highrecog_high_power{event,col}.powspctrm));
    recogdiff_high_power{event}  = zeros(size(all_recogdiff_high_power{event,col}.powspctrm));
       
end

total_sessions = 0;
for file = 1:length(listsq_files)
    if ~isempty(avg_lfp{1,file})
        total_sessions = total_sessions+1;
        for event = 1:2
            
            lfpsize = size(avg_lfp{event,file},1);
            if lfpsize == eventsize(event)
                lfps{event} = lfps{event} + avg_lfp{event,file};
            elseif lfpsize < eventsize(event) %should be larger
                disp('data size too small')
            else
                extra = lfpsize-eventsize(event);
                lfps{event} = lfps{event} + avg_lfp{event,file}(1:end-extra,:);
            end
            high_powerspctrm{event} = high_powerspctrm{event} +  all_high_powerspctrm{event,file}.powspctrm;
            low_powerspctrm{event} = low_powerspctrm{event} +  all_low_powerspctrm{event,file}.powspctrm;
            
            percent_change = [percent_change all_percent_change{file}];
            
            coherence{event} = coherence{event}+all_coherence{event,file}.cohspctrm;
            
            low_recog_coherence{event} = low_recog_coherence{event} + all_low_recog_coherence{event,file}.cohspctrm;
            high_recog_coherence{event} = high_recog_coherence{event} + all_high_recog_coherence{event,file}.cohspctrm;
            recogdiff_coherence{event} = recogdiff_coherence{event} + all_recogdiff_coherence{event,file}.cohspctrm;
            
            low_recog_coherence2High{event} = low_recog_coherence2High{event} + all_low_recog_coherence2High{event,file}.cohspctrm;
            high_recog_coherence2High{event} = high_recog_coherence2High{event} + all_high_recog_coherence2High{event,file}.cohspctrm;
            recogdiff_coherence2High{event} = recogdiff_coherence2High{event} + all_recogdiff_coherence2High{event,file}.cohspctrm;
            
            lowrecog_low_power{event} = lowrecog_low_power{event}+ all_lowrecog_low_power{event,file}.powspctrm;
            highrecog_low_power{event} = highrecog_low_power{event}+ all_highrecog_low_power{event,file}.powspctrm;
            recogdiff_low_power{event} = recogdiff_low_power{event}+ all_recogdiff_low_power{event,file}.powspctrm;
            
            lowrecog_high_power{event} = lowrecog_high_power{event}+ all_lowrecog_high_power{event,file}.powspctrm;
            highrecog_high_power{event} = highrecog_high_power{event}+ all_highrecog_high_power{event,file}.powspctrm;
            recogdiff_high_power{event} = recogdiff_high_power{event}+ all_recogdiff_high_power{event,file}.powspctrm;
            
        end
        %for fixation durations
        novlow = [novlow; all_fixdurs{1,file}];
        novhigh = [novhigh; all_fixdurs{2,file}];
        replow = [replow; all_fixdurs{3,file}];
        rephigh = [rephigh; all_fixdurs{4,file}];
        
    end
end

%---Plot fixation durations by low vs high memory---%
figure
hold on
plot(nanmean(novlow),'b')
plot(nanmean(novhigh),'g')
plot(nanmean(replow),'r')
plot(nanmean(rephigh),'m')
errorb(nanmean(novlow),nanstd(novlow)./sqrt(total_sessions));
errorb(nanmean(novhigh),nanstd(novhigh)./sqrt(total_sessions));
errorb(nanmean(replow),nanstd(replow)./sqrt(total_sessions));
errorb(nanmean(rephigh),nanstd(rephigh)./sqrt(total_sessions));
hold off
xlabel('Fixation #')
ylabel('Fixation Duration (ms)')
legend('Novel low','Novel high','Rep low','Rep high')
title('Fixation durations for high vs low recognition memory')
save_and_close_fig(figure_dir,'Fixation durations for high vs low recognition memory');

%---Plot Average LFP waveform across sessions---%
figure
subplot(1,2,1)
plot(timelockfix(1:eventsize(1)),lfps{1}/total_sessions)
xlim([-0.75 1.0])
line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
xlabel('Time (sec)')
ylabel('Voltage (uV)')
title('Aligned to Fix on Cross Across Sessions')
yl1 = ylim;

subplot(1,2,2)
plot(timelockimage(1:eventsize(2)),lfps{2}/total_sessions)
xlim([-0.75 1.0])
line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
xlabel('Time (sec)')
ylabel('Voltage (uV)')
title('Aligned to Image Onset Across Sessions')
yl2 = ylim;

%resacale plots to same height
ylmin = min([yl1 yl2]);
ylmax = max([yl1 yl2]);
subplot(1,2,1)
ylim([ylmin ylmax]);
subplot(1,2,2)
ylim([ylmin ylmax]);
save_and_close_fig(figure_dir,['SessionAverageLFP_PW'])


%---Plot Power spectrum of LFP---%
c1 = NaN(2,2);
c2 = NaN(2,2);

cfgfrq = [];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

figure
freqs = all_high_powerspctrm{1,2};
freqs.powspctrm = high_powerspctrm{1}/total_sessions;
subplot(2,2,1)
ft_singleplotTFR(cfgfrq, freqs);
title('Aligned to Fix on Cross Across Sessions')
c1(1,:) = caxis;

freqs = all_low_powerspctrm{1,2};
freqs.powspctrm = low_powerspctrm{1}/total_sessions;
subplot(2,2,3)
ft_singleplotTFR(cfgfrq, freqs);
c2(2,:) = caxis;

freqs = all_high_powerspctrm{2,2};
freqs.powspctrm = high_powerspctrm{2}/total_sessions;
subplot(2,2,2)
ft_singleplotTFR(cfgfrq, freqs);
title('Aligned to Image Onset Across Sessions')
c1(1,:) = caxis;

freqs = all_low_powerspctrm{2,2};
freqs.powspctrm = low_powerspctrm{2}/total_sessions;
subplot(2,2,4)
ft_singleplotTFR(cfgfrq, freqs);
c2(2,:) = caxis;

%rescale axis
minc1 = min(c1(:,1));
maxc1 = max(c1(:,2));
minc2 = min(c2(:,1));
maxc2 = max(c2(:,2));

subplot(2,2,1)
caxis([minc1 maxc1])
subplot(2,2,2)
caxis([minc1 maxc1])
subplot(2,2,3)
caxis([minc2 maxc2])
subplot(2,2,4)
caxis([minc2 maxc2])
save_and_close_fig(figure_dir,['SessionAveragePowerSpectrum_PW'])

%---Plot Histogram of % Change in Fixation Duration---%
figure
hist(percent_change,100);
xlabel('% Change in Fixation Durations from novel to repeat')
ylabel('Count')
xlim([-100 400])
save_and_close_fig(figure_dir,['SessionPercentChange_PW'])


%---Plot Phase Coherence---%

c1 = NaN(2,2);
figure

subplot(1,2,1)
imagesc(statlow.time,statlow.freq,squeeze(abs(coherence{1}(1,:,:)))/total_sessions)
axis xy
xlim([-0.3 0.5])
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
colorbar
colormap jet
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('Aligned to Fixation on Cross')
c1(1,:) = caxis;

subplot(1,2,2)
imagesc(statlow.time,statlow.freq,squeeze(abs(coherence{2}(1,:,:)))/total_sessions)
axis xy
xlim([-0.3 0.5])
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
colorbar
colormap jet
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('Aligned to Image Onset')
c1(2,:) = caxis;

minc1 = min(c1(:,1));
maxc1 = max(c1(:,2));

subplot(1,2,1)
caxis([minc1 maxc1])
subplot(1,2,2)
caxis([minc1 maxc1])

save_and_close_fig(figure_dir,['SessionPhaseCoherence_PW'])

%---Plot Phase Coherence for low vs high  recogntion---%
plotnum = [1 3 5; 2 4 6];
c1 = cell(1,2);

figure
for event = 1:2
    
    subplot(3,2,plotnum(event,1))
    imagesc(statlow.time,statlow.freq,squeeze(abs(low_recog_coherence{event}(1,:,:)))/total_sessions)
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Low Recog: Aligned to Fixation on Cross')
    c1{event} = [c1{event}; caxis];
    
    subplot(3,2,plotnum(event,2))
    imagesc(statlow.time,statlow.freq,squeeze(abs(high_recog_coherence{event}(1,:,:)))/total_sessions)
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('High Recog: Aligned to Fixation on Cross')
    c1{event} = [c1{event}; caxis];
    
    subplot(3,2,plotnum(event,3))
    imagesc(statlow.time,statlow.freq,squeeze(abs(recogdiff_coherence{event}(1,:,:)))/total_sessions)
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('High - Low: Aligned to Fixation on Cross')
    
end

for event = 1:2
    
    minc1 = min(c1{event}(:,1));
    maxc1 = max(c1{event}(:,2));
    
    subplot(3,2,plotnum(event,1))
    caxis([minc1 maxc1])
    
    subplot(3,2,plotnum(event,2))
    caxis([minc1 maxc1])
    
end

save_and_close_fig(figure_dir,['SessionPhaseCoherence_recognitionmemory_PW'])


%---Plot Phase coherence2High for low vs high  recogntion---%
plotnum = [1 3 5; 2 4 6];
c1 = cell(1,2);

figure
for event = 1:2
    
    subplot(3,2,plotnum(event,1))
    imagesc(stathigh.time,stathigh.freq,squeeze(abs(low_recog_coherence2High{event}(1,:,:)))/total_sessions)
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Low Recog: Aligned to Fixation on Cross')
    c1{event} = [c1{event}; caxis];
    
    subplot(3,2,plotnum(event,2))
    imagesc(stathigh.time,stathigh.freq,squeeze(abs(high_recog_coherence2High{event}(1,:,:)))/total_sessions)
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('High Recog: Aligned to Fixation on Cross')
    c1{event} = [c1{event}; caxis];
    
    subplot(3,2,plotnum(event,3))
    imagesc(stathigh.time,stathigh.freq,squeeze(abs(recogdiff_coherence2High{event}(1,:,:)))/total_sessions)
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('High - Low: Aligned to Fixation on Cross')
    
end

for event = 1:2
    
    minc1 = min(c1{event}(:,1));
    maxc1 = max(c1{event}(:,2));
    
    subplot(3,2,plotnum(event,1))
    caxis([minc1 maxc1])
    
    subplot(3,2,plotnum(event,2))
    caxis([minc1 maxc1])
    
end

save_and_close_fig(figure_dir,['SessionPhasecoherence2High2High_recognitionmemory_PW'])

%---Plot Power spectrum of LFP for High vs Low Recognition---%
% for low frequencies
c1 = NaN(2,2);
c2 = NaN(2,2);

cfgfrq = [];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

plotnum = [1 3 5; 2 4 6];
c1 = cell(1,2);


figure
for event = 1:2
    
    freqs = all_low_powerspctrm{1,2};
    freqs.powspctrm = lowrecog_low_power{event}/total_sessions;
    subplot(3,2,plotnum(event,1))
    ft_singleplotTFR(cfgfrq, freqs);
    c1{event} = [c1{event}; caxis];
    if event == 1
        title('Low Recog: Aligned to Fix on Cross')
    else
        title('Low Recog: Aligned to Image Onset Across Sessions')
    end
    
    freqs = all_low_powerspctrm{1,2};
    freqs.powspctrm = highrecog_low_power{event}/total_sessions;
    subplot(3,2,plotnum(event,2))
    ft_singleplotTFR(cfgfrq, freqs);
    c1{event} = [c1{event}; caxis];
    if event == 1
        title('High Recog: Aligned to Fix on Cross')
    else
        title('High Recog: Aligned to Image Onset Across Sessions')
    end
    
    freqs = all_low_powerspctrm{1,2};
    freqs.powspctrm = recogdiff_low_power{event}/total_sessions;
    subplot(3,2,plotnum(event,3))
    ft_singleplotTFR(cfgfrq, freqs);
    c1{event} = [c1{event}; caxis];
    if event == 1
        title('Memory Effect: Aligned to Fix on Cross')
    else
        title('Memory Effect: Aligned to Image Onset Across Sessions')
    end
    
    
end

for event = 1:2
    
    minc1 = min(c1{event}(:,1));
    maxc1 = max(c1{event}(:,2));
    
    subplot(3,2,plotnum(event,1))
    caxis([minc1 maxc1])
    
    subplot(3,2,plotnum(event,2))
    caxis([minc1 maxc1])
end
save_and_close_fig(figure_dir,['SessionLowFreqMemory_recognitionmemory_PW'])

% for high freq now
c1 = NaN(2,2);
c2 = NaN(2,2);

cfgfrq = [];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

plotnum = [1 3 5; 2 4 6];
c1 = cell(1,2);


figure
for event = 1:2
    
    freqs = all_high_powerspctrm{1,2};
    freqs.powspctrm = lowrecog_high_power{event}/total_sessions;
    subplot(3,2,plotnum(event,1))
    ft_singleplotTFR(cfgfrq, freqs);
    c1{event} = [c1{event}; caxis];
    if event == 1
        title('Low Recog: Aligned to Fix on Cross')
    else
        title('Low Recog: Aligned to Image Onset Across Sessions')
    end
    
    freqs = all_high_powerspctrm{1,2};
    freqs.powspctrm = highrecog_high_power{event}/total_sessions;
    subplot(3,2,plotnum(event,2))
    ft_singleplotTFR(cfgfrq, freqs);
    c1{event} = [c1{event}; caxis];
    if event == 1
        title('High Recog: Aligned to Fix on Cross')
    else
        title('High Recog: Aligned to Image Onset Across Sessions')
    end
    
    freqs = all_high_powerspctrm{1,2};
    freqs.powspctrm = recogdiff_high_power{event}/total_sessions;
    subplot(3,2,plotnum(event,3))
    ft_singleplotTFR(cfgfrq, freqs);
    c1{event} = [c1{event}; caxis];
    if event == 1
        title('Memory Effect: Aligned to Fix on Cross')
    else
        title('Memory Effect: Aligned to Image Onset Across Sessions')
    end
    
    
end

for event = 1:2
    
    minc1 = min(c1{event}(:,1));
    maxc1 = max(c1{event}(:,2));
    
    subplot(3,2,plotnum(event,1))
    caxis([minc1 maxc1])
    
    subplot(3,2,plotnum(event,2))
    caxis([minc1 maxc1])
end
save_and_close_fig(figure_dir,['SessionHighFreqMemory_recognitionmemory_PW'])