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
minfixdur = 100;%how long of fixation is required after a saccade

Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity

%row 1 locked to fixation cross hair, row 2 locked to image onset, col by file
%all unfiltered/sorted trials
% avg_lfp = cell(2,length(listsq_files));
% all_high_powerspctrm = cell(2,length(listsq_files));
% all_low_powerspctrm = cell(2,length(listsq_files));
% all_coherence = cell(2,length(listsq_files));
% all_coherence2High = cell(2,length(listsq_files));


for file = 3:length(listsq_files)
    disp(['File #' num2str(file) ': ' listsq_files{file}(1:end-11) ])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import & Preprocess Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([data_dir,listsq_files{file}(1:end-11)  '-preprocessed'],'cfg','fixationstats','item_set');
    
    %get important task specific information
    [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
    [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
    
    if str2double(item_set(7:8)) <= 6
        maximgdur = 7000*1.25;
        %images were up for 7 seconds cumulative
    else
        maximgdur = 5000*1.25;
        %images were up for 5 seconds cumulative
    end
    
    
    %preallocate space and parallel structure of cfg
    successful_trials = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if any(cfg.trl(t).allval == imgon_code) %in which image was displayed or 1st item in sequence was displayed
            if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)% sequence trials and only want rewarded ones
                imgdur = cfg.trl(t).alltim(cfg.trl(t).allval == imgoff_code)-cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code);
                if imgdur < maximgdur %so didn't look away too too much, trying to only take the best data
                    successful_trials(t) = t;
                end               
            end
        end
    end
    successful_trials = laundry(successful_trials);
    fixationstats = fixationstats(successful_trials);
    cfg.trl = cfg.trl(successful_trials);
    num_trials = length(successful_trials);
    
    
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
    which_img = which_img(trials_with_img_pairs);
    novel_vs_repeat = novel_vs_repeat(trials_with_img_pairs);
    num_trials = length(which_img);
    
    saccade_times = NaN(num_trials,75);
    for t = 1:length(fixationstats);
        if any(cfg.trl(t).allval == imgon_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code);%to ignore 1st fixation
            saccadetimes = fixationstats{t}.saccadetimes;
            fixationtimes = fixationstats{t}.fixationtimes;
            if ~isempty(fixationtimes)
                %ignore saccades before the first 250 ms of image on to
                %remove image onset confounds
                
                if saccadetimes(1,1) < fixationtimes(1,1) %saccade is the first index
                    sacindex = 0;
                else
                    sacindex = +1;
                end
                timeoutside = cumsum(isnan(fixationstats{t}.XY(1,:)));
                timeoutside = find(timeoutside >= 1000); %looked away a lot, again trying to cut out more data but with reason
                if ~isempty(timeoutside)
                    timeoutside = timeoutside(1);
                    toolate = find(saccadetimes(1,:) > timeoutside);
                    saccadetimes(:,toolate) = NaN;
                    toolate = find(fixationtimes(1,:) > timeoutside);
                    fixationtimes(:,toolate) = NaN;
                end
                
                tooearly = find(saccadetimes(1,:) <= imgon-trial_start+250);
                saccadetimes(:,tooearly) = NaN;
                tooearly = find(fixationtimes(1,:) <= imgon-trial_start+250);
                fixationtimes(:,tooearly) = NaN;
                
                %time outside image is ignored so not 1:1 saccaddes to
                %fixations
                fixdur =  NaN(1,size(saccadetimes,2));
                for sac = 1:size(saccadetimes,2);
                    sacend = saccadetimes(2,sac);
                    nextfix = find(fixationtimes(1,:) > sacend);
                    if ~isempty(nextfix)
                        nextfix = nextfix(1);
                        if (fixationtimes(1,nextfix) - sacend) == 1 %should be == 1
                            fixdur(sac) = fixationtimes(2,nextfix)-fixationtimes(1,nextfix)+1;
                        end
                    end
                end
                
                goodsacs = find(~isnan(saccadetimes(1,:)) & fixdur > minfixdur);
                saccade_times(t,1:length(goodsacs)) = saccadetimes(1,goodsacs);
            end
        end
    end
    sumnan = nansum(saccade_times);
    badcols = find(sumnan == 0);
    saccade_times(:,badcols)  = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    buffer1 = 2048;
    buffer2 = 2048;
    
    
    nov_saccade_aligned = NaN(ceil(num_trials/2),3);
    rep_saccade_aligned = NaN(ceil(num_trials/2),3);
    nov_index = 1;
    rep_index = 1;
    for trl = 1:length(cfg.trl);
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        
        sacs = saccade_times(trl,:);
        maxind = find(~isnan(sacs));
        for sac = 1:max(maxind)
            if ~isnan(saccade_times(sac))
                if novel_vs_repeat(trl) == 1 %novel image
                    nov_saccade_aligned(nov_index,:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                        saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                    nov_index = nov_index+1;
                else %repeat image
                    rep_saccade_aligned(rep_index,:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                        saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                    rep_index = rep_index+1;
                end
            end
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
    cfg.trl = nov_saccade_aligned;
    novsac_aligneddata = ft_preprocessing(cfg);
    
    %LFP data aligned to image onset
    cfg.trl = rep_saccade_aligned;
    repsac_aligneddata = ft_preprocessing(cfg);
    
    
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
    figurename = [listsq_files{file}(1:end-11) '_novel_saccade_lockedLFP'];
    timelock = timelocklfp(cfg,novsac_aligneddata,figure_dir,figurename);
    avg_lfp{1,file} = timelock.avg';
    timelockfix = timelock.time;
    
    %LFP data aligned to image onset
    figurename = [listsq_files{file}(1:end-11) '_repeat_saccade_lockedLFP'];
    timelock = timelocklfp(cfg,repsac_aligneddata,figure_dir,figurename);
    avg_lfp{2,file} = timelock.avg';
    timelockimage = timelock.time;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Low Frequency Power Spectrum---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = NaN(2,2);
    
    figure
    
    freq = lfp_powerspectrum(novsac_aligneddata,'low','all',[-0.5 -0.25]);
    all_low_powerspctrm{1,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,1)
    ft_singleplotTFR(cfgfrq, freq);
    title('Novel Saccades')
    c(1,:) = caxis;
    
    
    freq = lfp_powerspectrum(repsac_aligneddata,'low','all',[-1.25 -1.0]);
    all_low_powerspctrm{2,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,2)
    ft_singleplotTFR(cfgfrq, freq);
    title('Repeat Saccades')
    c(2,:) = caxis;
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    subtitle('Low Frequency Power Spectrum aligned to Saccades')
    save_and_close_fig(figure_dir,['ListSQ_List_Saccade_LowFreqPowerSpectrum_' listsq_files{file}(1:end-11)])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate High Frequency Power Spectrum---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(2,2);
    
    figure
    
    freq = lfp_powerspectrum(novsac_aligneddata,'high','all',[-0.5 -0.25]);
    all_high_powerspctrm{1,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,1)
    ft_singleplotTFR(cfgfrq, freq);
    title('Novel Saccades')
    c(1,:) = caxis;
    
    
    freq = lfp_powerspectrum(repsac_aligneddata,'high','all',[-1.25 -1.0]);
    all_high_powerspctrm{2,file} = freq;
    
    cfgfrq = [];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(1,2,2)
    ft_singleplotTFR(cfgfrq, freq);
    title('Repeat Saccades')
    c(2,:) = caxis;
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    subtitle('High Frequency Power Spectrum aligned to Saccades')
    save_and_close_fig(figure_dir,['ListSQ_List_Saccade_HighFreqPowerSpectrum_' listsq_files{file}(1:end-11)])
    clear freq;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Low Frequency Phase Coherence---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(2,2);
    figure
    
    stat = lfp_phasecoherence(novsac_aligneddata,'all');
    all_coherence{1,file} = stat;
    
    subplot(1,2,1)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Novel Saccades')
    c(1,:) = caxis;
    
    stat = lfp_phasecoherence(repsac_aligneddata,'all');
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
    title('Repeat Saccades')
    c(2,:) = caxis;
    
    subtitle('Low Freqyency Phase Coherence Aligned to Saccades')
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    save_and_close_fig(figure_dir,['ListSQ_List_Saccade_PhaseCoherence_' listsq_files{file}(1:end-11)])
    clear stat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate High Frequency Phase coherence---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = NaN(2,2);
    figure
    
    stat = lfp_phasecoherence2High(novsac_aligneddata,'all');
    all_coherence2High{1,file} = stat;
    
    subplot(1,2,1)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Novel Saccades')
    c(1,:) = caxis;
    
    stat = lfp_phasecoherence2High(repsac_aligneddata,'all');
    all_coherence2High{2,file} = stat;
    
    subplot(1,2,2)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Repeat Saccades')
    c(2,:) = caxis;
    
    subtitle('High Freqency Phase coherence2High Aligned to Saccades')
    
    %rescale so everything is the same
    minc = min(c(:,1));
    maxc = max(c(:,2));
    subplot(1,2,1)
    caxis([minc maxc])
    subplot(1,2,2)
    caxis([minc maxc]);
    save_and_close_fig(figure_dir,['ListSQ_List_Saccade_PhaseCoherence2High_' listsq_files{file}(1:end-11)])
    clear stat;
    
end