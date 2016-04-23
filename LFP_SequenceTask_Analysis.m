clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures\';

% multiunits = {[0 0 1],... %for plotting and other purpose 1 mulitunit 0 single unit
%     [0 0 0 1 1 1 0 0 0],[zeros(1,9) 1],[zeros(1,5) 1],[0 1 zeros(1,5) 1 zeros(1,5) 1],...
%     [0 0 0 1 0 0],[0 0 0],[0 0 0 0],[ 0 1 0 0 0 0],[zeros(1,7)],[0 0 0],[0 0 0 0],...
%     [zeros(1,4) 1],0,[1 0 0], [0 1 0],zeros(1,7),zeros(1,6),[zeros(1,8) 1 0],...
%     zeros(1,11),zeros(1,10),zeros(1,3),[zeros(1,4) ones(1,2) zeros(1,7)],...
%     zeros(1,5),zeros(1,7),zeros(1,8),[1 0 0 1],zeros(1,9),[0 0 0 1 0 0 0],...
%     zeros(1,8),zeros(1,4),[zeros(1,5) 1 0],zeros(1,3),[1 0],[0 0],...
%     zeros(1,9),zeros(1,8),zeros(1,3),[0 1 1 1],zeros(1,20)};

% cch25_files = {'PW140725_2.nex','PW140728_2_2-sorted.nex','PW140729_2-sorted.nex',...
%     'PW140730_2-sorted.nex','PW140801_2.nex','PW140805_2.nex','PW140806_2.nex',...
%     'PW140825_1.nex','PW140826_1.nex','PW140829_1.nex','PW140908_1.nex','PW140910_1.nex',...
%     'PW140915_1.nex','PW140917_1.nex','PW140101_1.nex','PW141006_1.nex',...
%     'PW141007_1.nex','PW141013_1.nex','PW141014_1.nex','PW141009_1.nex','PW141015_1.nex',...
%     'PW141008_1.nex','PW141010_1.nex','PW141024_1.nex','PW141027_1.nex','PW141028_1.nex',...
%     'PW141029_1.nex','PW141023_1.nex','PW141031_1.nex','PW141103_1.nex','PW141105_1.nex',...
%     'PW141106_1.nex','PW141110_1.nex','PW150121_1.nex','PW150122_1.nex','PW150123_1.nex',...
%     'PW150126_1.nex','PW150127_1.nex','PW150204_1.nex','PW150205_1.nex'};
% listsq_files = {'PW140725_3-sorted.nex','PW140728_2_3-sorted.nex','PW140729_3-sorted.nex',...
%     'PW140730_3-sorted.nex','PW140801_3-sorted.nex','PW140805_3-sorted.nex','PW140806_3-sorted.nex',...
%     'PW140825_3-sorted.nex','PW140826_3-sorted.nex','PW140829_3-sorted.nex','PW140908_3-sorted.nex',...
%     'PW140910_3-sorted.nex','PW140915_3-sorted.nex','PW140917_3-sorted.nex','PW140101_3-sorted.nex',...
%     'PW141006_3-sorted.nex','PW141007_3-sorted.nex','PW141013_3-sorted.nex','PW141014_3-sorted.nex',...
%     'PW141009_3-sorted.nex','PW141015_3-sorted.nex','PW141008_3-sorted.nex','PW141010_3-sorted.nex',...
%     'PW141024_3-sorted.nex','PW141027_3-sorted.nex','PW141028_3-sorted.nex','PW141029_3-sorted.nex',...
%     'PW141023_3-sorted.nex','PW141031_3-sorted.nex','PW141103_3-sorted.nex','PW141105_3-sorted.nex',...
%     'PW141106_3-sorted.nex','PW141110_3-sorted.nex','PW150121_3-sorted.nex','PW150122_3-sorted.nex',...
%     'PW150123_3-sorted.nex','PW150126_3-sorted.nex','PW150127_3-sorted.nex','PW150204_3-sorted.nex',...
%     'PW150205_3-sorted.nex'};
% item_sets =  {'ListSQ02.itm','ListSQ03.itm','ListSQ04.itm','ListSQ05.itm','ListSQ06.itm',...
%     'ListSQ07.itm','ListSQ08.itm','ListSQ09.itm','ListSQ10.itm','ListSQ12.itm','ListSQ13.itm',...
%     'ListSQ14.itm','ListSQ15.itm','ListSQ16.itm','ListSQ17.itm','ListSQ19.itm',...
%     'ListSQ20.itm','ListSQ23.itm','ListSQ24.itm','ListSQ21.itm','ListSQ25.itm',...
%     'ListSQ21.itm','ListSQ22.itm','ListSQ29.itm','ListSQ30.itm','ListSQ31.itm',...
%     'ListSQ32.itm','ListSQ28.itm','listSQ33.itm','listsq34.itm','Listsq35.itm',...
%     'ListSQ37.itm','ListSQ39.itm','ListSQVA.itm','ListSQVB.itm','ListSQVC.itm',...
%     'ListSQVD.itm','ListSQVE.itm','ListSQVH.itm','ListSQVI.itm'};

listsq_files = {'TO160217_3-sorted.nex'};
item_sets = {'ListSQ34.itm'};

% Important parametrs/values
twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
predicted_thresh = 10;% percent of saccades that must be < thresh ms to item to constitue as a learned-predicted item
% for Vivian 5% was the false positive rate for 150 ms so wanted at least double the False positive rate
predicted_rt = 138;%maximum "reaction time" for what constitutes as predictive, everything else is reactive
Fs = 1000;

Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity

for file = 1:length(listsq_files)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import & Preprocess Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([data_dir,listsq_files{file}(1:end-11)  '-preprocessed'],'cfg','fixationstats','item_set');
    
    %get important task specific information
    [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
    overlap = find((sequence_locations{1}(1,:) == sequence_locations{2}(1,:)) &...
        (sequence_locations{1}(2,:) == sequence_locations{2}(2,:)));
    
    %preallocate space and parallel structure of cfg
    successful_sequence_trials = NaN(1,length(cfg.trl));
    which_sequence = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
            which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            successful_sequence_trials(t) = t;
        end
    end
    successful_sequence_trials = laundry(successful_sequence_trials);
    which_sequence = laundry(which_sequence);
    
    num_trials = length(successful_sequence_trials);
    all_event_times = NaN(length(successful_sequence_trials),8);
    for t = 1:num_trials
        trial_start = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == trial_start_code);
        for event = 1:length(event_codes);
            all_event_times(t,event) = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == event_codes(event))-trial_start;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---process eye data locked to trial events---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixationstats = fixationstats(successful_sequence_trials);
    cfg.trl = cfg.trl(successful_sequence_trials);
    
    saccade_start_time = NaN(length(fixationstats),4);%when did saccade to item start
    fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
    reaction_time = NaN(length(fixationstats),4); %how long after a stimulus disappears do they saccade
    time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
    fixation_accuracy = NaN(length(fixationstats),4); %how far off
    fixation_duration = NaN(length(fixationstats),4); %fixation duration
    extrafixations = NaN(length(fixationstats),4); %how many noncross hair fixations do they mak
    for trial = 1:num_trials
        locs = sequence_locations{which_sequence(trial)};
        
        %convert to DVA for this analysis
        locs(1,:) = (locs(1,:)-400)/24;
        locs(2,:) = (locs(2,:)-300)/24;
        fixationstats{trial}.fixations(1,:) = (fixationstats{trial}.fixations(1,:)-400)/24;
        fixationstats{trial}.fixations(2,:) = (fixationstats{trial}.fixations(2,:)-300)/24;
        
        event_codes = cfg.trl(trial).allval;
        event_codes(event_codes == 100)= 0;
        event_codes(1) = 100;%eye data starts for recording right away
        event_times = cfg.trl(trial).alltim;
        
        trialdata = analyze_sequence_trial(fixationstats{trial},locs,fixwin,...
            event_codes,event_times);
        
        time_to_fixation(trial,:) = trialdata.t2f;
        fixation_duration(trial,:) = trialdata.fixation_duration;
        reaction_time(trial,:)  = trialdata.time_to_leave;
        extrafixations(trial,:) = trialdata.extrafixations;
        fixation_accuracy(trial,:) =  trialdata.accuracy;
        
        fixation_numbers = trialdata.fixationnums; %fixation number for each item
        fixationtimes = fixationstats{trial}.fixationtimes;
        saccadetimes = fixationstats{trial}.saccadetimes;
        
        for item = 1:4
            if ~isnan(fixation_numbers(item))
                
                fixation_start = fixationtimes(1,fixation_numbers(item));
                saccadeind = find(saccadetimes(2,:)+1 ==  fixation_start);
                
                if ~isempty(saccadeind)
                    saccade_start_time(trial,item) = saccadetimes(1,saccadeind);
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    buffer1 = 2048;
    buffer2 = 2048;
    
    %keep track of sequences and item number for later analysis and store
    %by event/eye movement
    itemnums = NaN(1,4*num_trials);%item number for event/eyemovement
    predicted = NaN(1,4*num_trials);%1 if less than RT thresh 0 if greater than
    context = NaN(1,4*num_trials); %sequence 1 or 2
    saccade_aligned = NaN(4*num_trials,3); %align to saccade to item
    fixation_aligned = NaN(4*num_trials,3); %item aligned to fixation on item
    item_aligned = NaN(4*num_trials,3); %aligned to item appearing
    index = 1;
    for trl = 1:size(saccade_start_time,1);
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        for item = 1:4
            if ~isnan(saccade_start_time(trl,item))
                itemnums(index) = item;
                context(index) = which_sequence(trl);
                if time_to_fixation(trl,item) < predicted_rt
                    predicted(index) = 1;
                else
                    predicted(index) = 0;
                end
                
                saccade_aligned(index,:) = [saccade_start_time(trl,item)-1-buffer1+startind ...
                    saccade_start_time(trl,item)-1+buffer2+startind -buffer1];
                fixation_aligned(index,:) = [fixation_start_time(trl,item)-1-buffer1+startind ...
                    fixation_start_time(trl,item)-1+buffer2+startind -buffer1];
                item_aligned(index,:) = [all_event_times(trl,2*item-1)-1-buffer1+startind ...
                    all_event_times(trl,2*item-1)-1+buffer2+startind -buffer1];
                
                index = index+1;
            end
        end
    end
    
    %clean up/remove any NaNs (occur when saccades couldn't be detected)
    itemnums = laundry(itemnums);
    context = laundry(context);
    predicted = laundry(predicted);
    saccade_aligned = laundry(saccade_aligned);
    fixation_aligned = laundry(fixation_aligned);
    item_aligned = laundry(item_aligned);
    
    
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
    cfg.trl = saccade_aligned;
    saccade_aligneddata = ft_preprocessing(cfg);
    
    %fixation aligned LFP data
    cfg.trl = fixation_aligned;
    fixation_aligneddata = ft_preprocessing(cfg);
    
    %item aligned LFP data
    cfg.trl = item_aligned;
    item_aligneddata = ft_preprocessing(cfg);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Time Locked Analysis---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear cfg
    cfg.channel       = 'all';
    cfg.covariance    = 'no';
    cfg.keeptrials    = 'yes';
    cfg.removemean    = 'yes';
    cfg.vartrllength  =  2;
    
    %saccade aligned LFP data
    timelock = ft_timelockanalysis(cfg,saccade_aligneddata);
    figure
    plot(timelock.time,timelock.avg')
    xlim([-0.75 1.0])
    line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
    xlabel('Time (sec)')
    ylabel('Voltage (uV)')
    title('Saccade Aligned')
    
    %fixation aligned LFP data
    timelock = ft_timelockanalysis(cfg,fixation_aligneddata);
    figure
    plot(timelock.time,timelock.avg')
    xlim([-0.75 1.0])
    line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
    xlabel('Time (sec)')
    ylabel('Voltage (uV)')
    title('fixation Aligned')
    
    %item aligned LFP data
    timelock = ft_timelockanalysis(cfg,item_aligneddata);
    figure
    plot(timelock.time,timelock.avg')
    xlim([-0.75 1.0])
    line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
    xlabel('Time (sec)')
    ylabel('Voltage (uV)')
    title('item Aligned')
    
    
end
%%
% Change in frquence/power over time by item
figure
freq = cell(1,4);
for item = 1:4
    
    clear cfgfrq
    cfgfrq.output      = 'pow';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.taper       = 'hanning';
    cfgfrq.foi         = 2:1:30;
    cfgfrq.tapsmofrq = 8;
    cfgfrq.t_ftimwin       = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = -0.75:.01:1.0;
    
    desired_trials = find(itemnums == item);
    data = saccade_aligneddata;
    data.sampleinfo = data.sampleinfo(desired_trials,:);
    data.time = data.time(desired_trials);
    data.trial = data.trial(desired_trials);
    data.cfg.trl = data.cfg.trl(desired_trials,:);
    freq{item} = ft_freqanalysis(cfgfrq, data);
    
    cfgfrq = [];
    % cfgfrq.baseline     = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    % cfgfrq.channel      = 'AD01';
    cfgfrq.zparam       = 'powspctrm';
    
    subplot(2,2,item)
    ft_singleplotTFR(cfgfrq, freq{item});
    title('item Aligned')
end

sb = 1:16;
sb = reshape(sb,[4,4])';
figure
for item = 1:4
    for it = 1:4
        if item < it
            frq = freq{1};
            frq.powspctrm = freq{item}.powspctrm-freq{it}.powspctrm;
            subplot(4,4,sb(item,it))
            ft_singleplotTFR(cfgfrq, frq);
            title('Saccade Aligned')
        end
    end
end
subtitle('Power Spectrum')

%% calculate inter-item coherence by item

freqc = cell(1,4);
stat_nov = cell(1,4);

figure
for item = 1:4
    
    clear cfgfrq
    cfgfrq.output      = 'fourier';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.foi         = 4:2:24;
    numfoi             = length(cfgfrq.foi);
    cfgfrq.taper       = 'hanning';
    cfgfrq.pad         = 'maxperlen';
    cfgfrq.keeptrials  = 'yes';
    cfgfrq.keeptapers  = 'yes';
    cfgfrq.complex     = 'complex';
    cfgfrq.t_ftimwin   = 0.5 * ones(1,numfoi);
    cfgfrq.toi         = -0.5:.01:0.75;
    
    desired_trials = find(itemnums == item);
    data = saccade_aligneddata;
    data.sampleinfo = data.sampleinfo(desired_trials,:);
    data.time = data.time(desired_trials);
    data.trial = data.trial(desired_trials);
    data.cfg.trl = data.cfg.trl(desired_trials,:);
    freqc{item} = ft_freqanalysis(cfgfrq, data);
    
    sizfrq=size(freqc{item}.label,1);
    freqc{item}.label{sizfrq+1} = 'fake';
    siz = size(freqc{item}.fourierspctrm);
    freqc{item}.fourierspctrm(:,sizfrq+1,:,:) = complex(ones(siz(1),1,siz(3),siz(4)), ...
        zeros(siz(1),1,siz(3),siz(4)));
    
    % create channelcmb
    lfpind=strmatch('AD',freqc{item}.label);
    cfgfrq.channelcmb=num2cell(zeros(length(lfpind),2)+nan);
    for k=1:length(lfpind)
        cfgfrq.channelcmb{k,1} = freqc{item}.label{lfpind(k)};
        cfgfrq.channelcmb{k,2} = 'fake';
    end
    
    % determine the descriptive spectral statistics
    cfgfrq.method      = 'coh';
    stat_nov{item} = ft_connectivityanalysis(cfgfrq, freqc{item});
    
    subplot(2,2,item)
    imagesc(stat_nov{item}.time,stat_nov{item}.freq,squeeze(abs(stat_nov{item}.cohspctrm(1,:,:))))
    axis xy
    xlim([-0.3 0.5])
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar
    colormap jet
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
end
subtitle('Power Spectrum by Item')

figure
sb = 1:16;
sb = reshape(sb,[4,4])';
for item = 1:4
    for it = 1:4
        if item < it
            coh = stat_nov{item}.cohspctrm-stat_nov{it}.cohspctrm;
            subplot(4,4,sb(item,it))
            imagesc(stat_nov{1}.time,stat_nov{1}.freq,squeeze(abs(coh(1,:,:))))
            axis xy
            xlim([-0.3 0.5])
            line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
            colorbar
            colormap jet
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)')
        end
    end
end
subtitle('Difference in Power Spectrum by Item')
%%
% Change in frquence/power over time by predicted or not
freq = cell(1,2);
stat_nov = cell(1,2);
for status = 0:1
    
    clear cfgfrq
    cfgfrq.output      = 'pow';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.taper       = 'hanning';
    cfgfrq.foi         = 2:2:180;
    cfgfrq.tapsmofrq = 8;
    cfgfrq.t_ftimwin       = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = -0.75:.01:1.0;
    
    desired_trials = find(predicted(1:end-1) == status & itemnums(1:end-1) == 3);
    data = saccade_aligneddata;
    data.sampleinfo = data.sampleinfo(desired_trials,:);
    data.time = data.time(desired_trials);
    data.trial = data.trial(desired_trials);
    data.cfg.trl = data.cfg.trl(desired_trials,:);
    freq{status+1} = ft_freqanalysis(cfgfrq, data);
    
    cfgfrq = [];
    % cfgfrq.baseline     = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    % cfgfrq.channel      = 'AD01';
    cfgfrq.zparam       = 'powspctrm';
    
    figure(1)
    subplot(1,2,status+1)
    ft_singleplotTFR(cfgfrq, freq{status+1});
    if status == 0
        title('Reactive')
    else
        title('Predictive')
    end
    

    %for coherence
    clear cfgfrq
    cfgfrq.output      = 'fourier';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.foi         = 1:2:30;
    numfoi             = length(cfgfrq.foi);
    cfgfrq.taper       = 'hanning';
    cfgfrq.pad         = 'maxperlen';
    cfgfrq.keeptrials  = 'yes';
    cfgfrq.keeptapers  = 'yes';
    cfgfrq.complex     = 'complex';
    cfgfrq.t_ftimwin   = 0.5 * ones(1,numfoi);
    cfgfrq.toi         = -0.5:.01:0.75;
    
    freq2 = ft_freqanalysis(cfgfrq, data);
    
    sizfrq=size(freq2.label,1);
    freq2.label{sizfrq+1} = 'fake';
    siz = size(freq2.fourierspctrm);
    freq2.fourierspctrm(:,sizfrq+1,:,:) = complex(ones(siz(1),1,siz(3),siz(4)), ...
        zeros(siz(1),1,siz(3),siz(4)));
    
    % create channelcmb
    lfpind=strmatch('AD',freq2.label);
    cfgfrq.channelcmb=num2cell(zeros(length(lfpind),2)+nan);
    for k=1:length(lfpind)
        cfgfrq.channelcmb{k,1} = freq2.label{lfpind(k)};
        cfgfrq.channelcmb{k,2} = 'fake';
    end
    
    % determine the descriptive spectral statistics
    cfgfrq.method      = 'coh';
    stat_nov{status+1} = ft_connectivityanalysis(cfgfrq, freq2);
    
%     figure(2)
%     subplot(1,2,status+1)
%     if status == 0
%         title('Reactive')
%     else
%         title('Predictive')
%     end
%     imagesc(stat_nov.time,stat_nov.freq,squeeze(abs(stat_nov.cohspctrm(1,:,:))))
%     axis xy
%     xlim([-0.3 0.5])
%     line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
%     colorbar
%     colormap jet
%     xlabel('Time (sec)')
%     ylabel('frequency (Hz)')
end

figure(1)
subtitle('Power Spectrum')
figure(2)
subtitle('Phase Coherence')

figure
frq = freq{1};
frq.powspctrm = freq{2}.powspctrm-freq{1}.powspctrm;
ft_singleplotTFR(cfgfrq, frq);
title('Predictive-Reactive: Saccade Aligned')


figure
coh = squeeze(abs(stat_nov{1}.cohspctrm(1,:,:)));
coh2 = squeeze(abs(stat_nov{2}.cohspctrm(1,:,:)));
imagesc(stat_nov{1}.time,stat_nov{1}.freq,coh2-coh);
axis xy
xlim([-0.5 0.5])
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
colorbar
colormap jet
xlabel('Time (sec)')
ylabel('frequency (Hz)')
title('Predictive-Reactive: Saccade Aligned')
%% high frequencey phase coherence

freq = cell(1,2);
stat = cell(1,2);
for status = 0:1
    
    clear cfgfrq
    cfgfrq.output      = 'pow';
    cfgfrq.method      = 'mtmconvol';
    cfgfrq.taper       = 'hanning';
    cfgfrq.foi         = 2:2:180;
    cfgfrq.tapsmofrq = 8;
    cfgfrq.t_ftimwin       = 7./cfgfrq.foi;  % 7 cycles per time window
    cfgfrq.toi         = -0.75:.01:1.0;
    
    desired_trials = find(predicted(1:end-1) == status & itemnums(1:end-1) ~= 1);
    data = saccade_aligneddata;
    data.sampleinfo = data.sampleinfo(desired_trials,:);
    data.time = data.time(desired_trials);
    data.trial = data.trial(desired_trials);
    data.cfg.trl = data.cfg.trl(desired_trials,:);

  %Saccades during novel images during ListSQ List trials
    stat{status+1} = lfp_phasecoherence2High(data,'all');
    subplot(1,2,status+1)
    imagesc(stat{status+1}.time,stat{status+1}.freq,squeeze(abs(stat{status+1}.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
end
%%
% % do the spectral analysis - time-frequency
% clear cfgfrq
% cfgfrq.output      = 'fourier';
% cfgfrq.method      = 'mtmconvol';
% cfgfrq.foi         = 0:2:20;
% numfoi             = length(cfgfrq.foi);
% cfgfrq.taper       = 'hanning';
% cfgfrq.pad         = 'maxperlen';
% cfgfrq.keeptrials  = 'yes';
% cfgfrq.keeptapers  = 'yes';
% cfgfrq.complex     = 'complex';
% cfgfrq.t_ftimwin   = 0.5 * ones(1,numfoi);
% cfgfrq.toi         = -0.5:.01:0.75;
%
% freq = ft_freqanalysis(cfgfrq, saccade_aligneddata);
%
% sizfrq=size(freq.label,1);
% freq.label{sizfrq+1} = 'fake';
% siz = size(freq.fourierspctrm);
% freq.fourierspctrm(:,sizfrq+1,:,:) = complex(ones(siz(1),1,siz(3),siz(4)), ...
%     zeros(siz(1),1,siz(3),siz(4)));
%
% % create channelcmb
% lfpind=strmatch('AD',freq.label);
% cfgfrq.channelcmb=num2cell(zeros(length(lfpind),2)+nan);
% for k=1:length(lfpind)
%     cfgfrq.channelcmb{k,1} = freq.label{lfpind(k)};
%     cfgfrq.channelcmb{k,2} = 'fake';
% end
%
% % determine the descriptive spectral statistics
% cfgfrq.method      = 'coh';
% stat_nov = ft_connectivityanalysis(cfgfrq, freq);
%
% figure
% imagesc(stat_nov.time,stat_nov.freq,squeeze(abs(stat_nov.cohspctrm(1,:,:))))
% axis xy
% xlim([-0.3 0.5])
% line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
% colorbar
% colormap jet
% xlabel('Time (sec)')
% ylabel('Frequency (Hz)')
%%

% clear cfgfrq
% cfgfrq.output      = 'pow';
% cfgfrq.method      = 'mtmconvol';
% cfgfrq.taper       = 'hanning';
% cfgfrq.foi         = 0:2:30;
% numfoi             = length(cfgfrq.foi);
% cfgfrq.t_ftimwin   = 0.5 * ones(1,numfoi);
% cfgfrq.toi         = -0.75:.01:1.0;
%
% freq = ft_freqanalysis(cfgfrq, fixation_aligneddata);
%
% cfgfrq = [];
% % cfgfrq.baseline     = [-0.75 -0.5];
% cfgfrq.baselinetype = 'absolute';
% cfgfrq.maskstyle    = 'saturation';
% % cfgfrq.channel      = 'AD01';
% cfgfrq.zparam       = 'powspctrm';
% figure
% ft_singleplotTFR(cfgfrq, freq);

%%

clear cfgfrq
cfgfrq.output      = 'pow';
cfgfrq.method      = 'mtmconvol';
cfgfrq.taper       = 'hanning';
cfgfrq.foi         = 4:2:180;
cfgfrq.tapsmofrq = 8;
cfgfrq.t_ftimwin       = 7./cfgfrq.foi;  % 7 cycles per time window
cfgfrq.toi         = -0.75:.01:1.0;

freq = ft_freqanalysis(cfgfrq, saccade_aligneddata);

cfgfrq = [];
% cfgfrq.baseline     = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
% cfgfrq.channel      = 'AD01';
cfgfrq.zparam       = 'powspctrm';
figure
ft_singleplotTFR(cfgfrq, freq);
title('Saccade Aligned Power Spectrum2')

%% calculate inter-fixation coherence
%
% % do the spectral analysis - time-frequency
% clear cfgfrq
% cfgfrq.output      = 'fourier';
% cfgfrq.method      = 'mtmconvol';
% cfgfrq.foi         = 2:2:20;
% numfoi             = length(cfgfrq.foi);
% cfgfrq.taper       = 'hanning';
% cfgfrq.pad         = 'maxperlen';
% cfgfrq.keeptrials  = 'yes';
% cfgfrq.keeptapers  = 'yes';
% cfgfrq.complex     = 'complex';
% cfgfrq.t_ftimwin   = 0.5 * ones(1,numfoi);
% cfgfrq.toi         = -0.5:.01:0.75;
%
% freq = ft_freqanalysis(cfgfrq, saccade_aligneddata);
%
% sizfrq=size(freq.label,1);
% freq.label{sizfrq+1} = 'fake';
% siz = size(freq.fourierspctrm);
% freq.fourierspctrm(:,sizfrq+1,:,:) = complex(ones(siz(1),1,siz(3),siz(4)), ...
%     zeros(siz(1),1,siz(3),siz(4)));
%
% % create channelcmb
% lfpind=strmatch('AD',freq.label);
% cfgfrq.channelcmb=num2cell(zeros(length(lfpind),2)+nan);
% for k=1:length(lfpind)
%     cfgfrq.channelcmb{k,1} = freq.label{lfpind(k)};
%     cfgfrq.channelcmb{k,2} = 'fake';
% end
%
% % determine the descriptive spectral statistics
% cfgfrq.method      = 'coh';
% stat_nov = ft_connectivityanalysis(cfgfrq, freq);
%
% figure
% imagesc(stat_nov.time,stat_nov.freq,squeeze(abs(stat_nov.cohspctrm(1,:,:))))
% axis xy
% xlim([-0.3 0.5])
% line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
% colorbar
% colormap jet
% xlabel('Time (sec)')
% ylabel('Frequency (Hz)')