cfgfrq = [];
cfgfrq.baseline = 'no';
cfgfrq.baselinetype = [];
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';


%written 11/10 && 11/11 2015 by Seth Konig using/combining/modify code
%previously written to anlayze task speartely. Code focuses on eye movments
%or lack there off!!!!

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Across Tasks Figures\'; %where to save figures

%sessions with recordings on the same day
listsq_files = {'PW140729_3-sorted.nex','PW140730_3-sorted.nex','PW140806_3-sorted.nex',...
    'PW140829_3-sorted.nex','PW141007_3-sorted.nex','PW141008_3-sorted.nex',...
    'PW141009_3-sorted.nex','PW141015_3-sorted.nex','PW141024_3-sorted.nex',...
    'PW141028_3-sorted.nex','PW150205_3-sorted.nex'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important parametrs/values for ListSQ List
ITI_start = 15;
ITI_end = 16;
imgon_code = 23;
imgoff_code = 24;
fixspoton = 35;
imageX = 800; %horizontal size of the screen
imageY = 600; %horizontal size of the screen
minfixdur = 100;%how long of fixation is required after a saccade
maxoutside = 1000;%how long cumulatively have they looked outside the image before the data is ignored
%1000 gives them a reasonable buffer


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Important parameters/values for ListSQ Sequence
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
predicted_rt = 156;%maximum "reaction time" for what constitutes as predictive, everything else is reactive
reward_code = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set/get some general important info for CVTNEW task
dot_on_code = 25;
dot_clrchng_code = 27;
bar_code_response = 4; %monkey made a move
%time is + 300 ms
short = [700 1133]; %short duration trials
mediumm = [1134 1567];
long = [1568 2500];%cap should be 2000 but cortex can have some lag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Important parameters/values across all tasks
trial_start_code = 15;
fixation_on_cross = 8;
reward_code = 3;
Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity
Fs = 1000;
buffer1 = 2048;
buffer2 = 4096;


low_freq_data = cell(1,length(listsq_files));
low_freq_coherence_data = cell(1,length(listsq_files));
high_freq_data = cell(1,length(listsq_files));
high_freq_coherence_data = cell(1,length(listsq_files));
for file = 1:length(listsq_files)
    
    
    disp(['Running File #' num2str(file) ': ' listsq_files{file}(1:end-13) ])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Sequence portion of ListSQ---%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import & Preprocess Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %reload because wrote over cfg
    load([data_dir,listsq_files{file}(1:end-11)  '-preprocessed'],'cfg','fixationstats','item_set');
    
    %get important task specific information
    [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
    overlap = find((sequence_locations{1}(1,:) == sequence_locations{2}(1,:)) &...
        (sequence_locations{1}(2,:) == sequence_locations{2}(2,:)));
    
    %preallocate space and parallel structure of cfg
    successful_sequence_trials = NaN(1,length(cfg.trl));
    which_sequence = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if sum(cfg.trl(t).allval == reward_code) >= 5; %in which sequence trials were rewarded
            which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            successful_sequence_trials(t) = t;
        end
    end
    successful_sequence_trials = laundry(successful_sequence_trials);
    which_sequence = laundry(which_sequence);
    
    event_codes = [15]; %odd indeces item on, even indeces item off
    num_trials = length(successful_sequence_trials);
    all_event_times = NaN(length(successful_sequence_trials),1);
    for t = 1:num_trials
        trial_start = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == trial_start_code);
        for event = 1:length(event_codes);
            all_event_times(t,event) = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == event_codes(event))-trial_start;
        end
    end
    cfg.trl = cfg.trl(successful_sequence_trials);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    
    %keep track of sequences and item number for later analysis and store
    %by event/eye movement
    ITI_aligned = NaN(num_trials,3); %aligned to item appearing
    index = 1;
    for trl = 1:length(cfg.trl)
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        ITI_aligned(index,:) = [all_event_times(trl)-1-buffer1+startind ...
            all_event_times(trl)-1+buffer2+startind -buffer1];
        index = index+1;
    end
    
    %clean up/remove any NaNs (occur when saccades couldn't be detected)
    ITI_aligned = laundry(ITI_aligned);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import LFP data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %also filters out line noise
    
    %code to find spike/LFP channels
    LFPchannels = find_desired_channels(cfg,'LFP');
    
    % read the data from file and preprocess them
    cfg.channel       = cfg.channel(LFPchannels)';
    cfg.dftfilter     = 'yes';
    cfg.dftfreq       = Fline;
    cfg.padding       = 1;
    
    %saccade aligned LFP data
    ITIcfg = cfg;
    ITIcfg.trl = ITI_aligned;
    ITI_aligneddata = ft_preprocessing(ITIcfg);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Plot time locked LFP---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     timelock = ft_timelockanalysis(ITIcfg,ITI_aligneddata);
    %     plot(timelock.time,timelock.avg')
    %     xlim([-0.75 1.75])
    %     xlabel('Time from Dot on (sec)')
    %     ylabel('Voltage (uV)')
    %     title('ITI LFP')
    
    
    %---Plot All Low Freq Power Spectrum---%
    
    
    cfgfrq = [];
    cfgfrq.baseline = 'no';
    cfgfrq.baselinetype = [];
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    freq = lfp_powerspectrum(ITI_aligneddata,'low','all');
    low_freq_data{file} = freq;
    
    %     figure
    %     hold on
    %     ft_singleplotTFR(cfgfrq, freq);
    %     line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    %     xlim([-0.75 1.75])
    %     xlabel('Time from ITI start (sec)')
    %     hold off
    
    
    %---Low Freq Coherence---%
    
    stat = lfp_phasecoherence(ITI_aligneddata,'all');
    low_freq_coherence_data{file} = stat;
    
    
    %     figure
    %     imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    %     line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    %     colorbar, axis xy
    %     colormap jet
    %     ylabel('Frequency (Hz)')
    %     xlim([-0.75 1.75])
    %     xlabel('Time from ITI start (sec)')
    
    
    %---Plot All high Freq Power Spectrum---%
    
    
    cfgfrq = [];
    cfgfrq.baseline = 'no';
    cfgfrq.baselinetype = [];
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    freq = lfp_powerspectrum(ITI_aligneddata,'high','all');
    high_freq_data{file} = freq;
    
    %     figure
    %     hold on
    %     ft_singleplotTFR(cfgfrq, freq);
    %     line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    %     xlim([-0.75 1.75])
    %     xlabel('Time from Saccade (sec)')
    %     hold off
    %
    
    %---high Freq Coherence---%
    
    stat = lfp_phasecoherence2High(ITI_aligneddata,'all');
    high_freq_coherence_data{file} = stat;
    
    
    %     figure
    %     imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    %     line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    %     colorbar, axis xy
    %     colormap jet
    %     ylabel('Frequency (Hz)')
    %     xlim([-0.75 1.75])
    %     xlabel('Time from Saccade (sec)')
    
    
end
save('SequenceOnly_ITI_across_tasks')
emailme('Done with Across task data analysis')
%%

all_low_freq_data = [];
all_low_freq_coherence_data = [];
all_high_freq_data = [];
all_high_freq_coherence_data = [];
total_sessions = 0;
for file = 1:length(listsq_files);
    if ~isempty(low_freq_data{file})
        if file == 1
            all_low_freq_data= low_freq_data{file}.powspctrm;
            all_low_freq_coherence_data = low_freq_coherence_data{file}.cohspctrm;
            all_high_freq_data= high_freq_data{file}.powspctrm;
            all_high_freq_coherence_data = high_freq_coherence_data{file}.cohspctrm;f
        else
            all_low_freq_data = all_low_freq_data+low_freq_data{file}.powspctrm;
            all_low_freq_coherence_data =  all_low_freq_coherence_data+low_freq_coherence_data{file}.cohspctrm;
            all_high_freq_data = all_high_freq_data+high_freq_data{file}.powspctrm;
            all_high_freq_coherence_data =  all_high_freq_coherence_data+high_freq_coherence_data{file}.cohspctrm;
        end
        total_sessions = total_sessions+1;
    end
end

cfgfrq = [];
cfgfrq.baseline = 'no';
cfgfrq.baselinetype = [];
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

%%
%---Across Trial low frequency Power---%
freq = low_freq_data{1};
freq.powspctrm =  all_low_freq_data/total_sessions;

figure
hold on
ft_singleplotTFR(cfgfrq, freq);
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
hold off
xlim([-0.75 1.75])
xlabel('Time from ITI start (sec)')
%%
%---Across Trial High frequency Power---%
freq = high_freq_data{1};
freq.powspctrm =  all_high_freq_data/total_sessions;

figure
hold on
ft_singleplotTFR(cfgfrq, freq);
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
hold off
xlim([-0.75 1.75])
xlabel('Time from ITI start (sec)')

%%


figure

stat = low_freq_coherence_data{1};
stat.cohspctrm = all_low_freq_coherence_data/total_sessions;
imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
colorbar, axis xy
colormap jet
ylabel('Frequency (Hz)')
xlim([-0.75 1.75])
xlabel('Time from ITI start(sec)')
%%
figure

stat = high_freq_coherence_data{1};
stat.cohspctrm = all_high_freq_coherence_data/total_sessions;
imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
colorbar, axis xy
colormap jet
ylabel('Frequency (Hz)')
xlim([-0.75 1.75])
xlabel('Time from ITI start(sec)')

%%


baselinetypes = {'absolute','relchange','relative'};


baselinetype = baselinetypes{2};
baseline = [0.25 0.75];% baseline time interval; cfgfrq.baseline

timeVec = freq.time;

tidx = find(timeVec >= baseline(1) & timeVec <= baseline(2));

TFdata = all_low_freq_data;

TFbl = NaN(size(TFdata,1),size(TFdata,2));
for k=1:size(TFdata,2) % loop frequencies
    for l=1:size(TFdata,1) % loop channels
        TFbl(l,k) = nanmean(TFdata(l,k,tidx),3); %compute average baseline power
        if TFbl(l,k) == 0,
            error('Average baseline power is zero');
        end
    end
end

if strcmpi(baselinetype,'relative')
    for k=1:size(TFdata,2) % loop frequencies
        for l=1:size(TFdata,1) % loop channels
            TFdata(l,k,:) = TFdata(l,k,:) / TFbl(l,k);     % compute relative change (i.e. ratio)
        end
    end
    
elseif strcmpi(baselinetype,'absolute')
    for k=1:size(TFdata,2) % loop frequencies
        for l=1:size(TFdata,1) % loop channels
            TFdata(l,k,:) = TFdata(l,k,:) - TFbl(l,k);        % subtract baseline power
        end
    end
    
elseif strcmpi(baselinetype,'relchange')
    for k=1:size(TFdata,2) % loop frequencies
        for l=1:size(TFdata,1) % loop channels
            TFdata(l,k,:) = ((TFdata(l,k,:) - TFbl(l,k)) / TFbl(l,k)); % compute relative change
        end
    end
    
end