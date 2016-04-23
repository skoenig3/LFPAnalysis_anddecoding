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
buffer2 = 2048;


l_trials = NaN(1,length(listsq_files));
c_trials = NaN(1,length(listsq_files));

low_freq_data = cell(4,length(listsq_files));
low_freq_coherence_data = cell(4,length(listsq_files));

item_low_freq_data = cell(4,length(listsq_files));
item_low_freq_coherence_data = cell(4,length(listsq_files));
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
    
    event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
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
    time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
    for trial = 1:num_trials/4%shouldnt be divided by 4
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
        
        fixation_numbers = trialdata.fixationnums; %fixation number for each item
        fixationtimes = fixationstats{trial}.fixationtimes;
        saccadetimes = fixationstats{trial}.saccadetimes;
        if saccadetimes(1,1) < fixationtimes(1,1); %then started with a saccade
            sacind = 0;
        else%started with a fixation
            sacind = 1;
        end
        for item = 1:4
            if ~isnan(fixation_numbers(item))
                if fixation_numbers(item)+sacind  >= 1
                    try
                        saccade_start_time(trial,item) = saccadetimes(1,fixation_numbers(item)+sacind);
                    catch
                        saccade_start_time(trial,item) = NaN;%should only occur if saccade was off screen or a blink occured otherwise should work
                        disp('Saccade Start Time not found')
                    end
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    
    %keep track of sequences and item number for later analysis and store
    %by event/eye movement
    itemnums = NaN(1,4*num_trials);%item number for event/eyemovement
    predicted = NaN(1,4*num_trials);%1 if less than RT thresh 0 if greater than
    context = NaN(1,4*num_trials); %sequence 1 or 2
    saccade_aligned = NaN(4*num_trials,3); %align to saccade to item
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
    cfg.dftfreq       = Fline;
    cfg.padding       = 1;
    
    %saccade aligned LFP data
    seqcfg = cfg;
    seqcfg.trl = saccade_aligned;
    seq_saccade_aligneddata = ft_preprocessing(seqcfg);
    
    %item aligned LFP data
    itmcfg = cfg;
    itmcfg.trl = item_aligned;
    item_aligneddata = ft_preprocessing(itmcfg);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%---Plot All Low Freq Power Spectrum---%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     cfgfrq = [];
%     cfgfrq.baseline = 'no';
%     cfgfrq.baselinetype = [];
%     cfgfrq.maskstyle    = 'saturation';
%     cfgfrq.zparam       = 'powspctrm';
%     
%     yl = NaN(4,2);
%     
%     figure
%     for it = 1:4
%         freq = lfp_powerspectrum(seq_saccade_aligneddata,'low',find(itemnums == it),[-0.65 -0.4]);
%         low_freq_data{it,file} = freq;
%         
%         subplot(2,2,it)
%         hold on
%         ft_singleplotTFR(cfgfrq, freq);
%         hold off
%         line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
%         xlim([-0.75 0.75])
%         xlabel('Time from Saccade (sec)')
%         yl(it,:) = caxis;
%         title(['Item #' num2str(it)])
%     end
%     
%     ymin = min(yl(:,1));
%     ymax = max(yl(:,2));
%     for it = 1:4
%         subplot(2,2,it)
%         caxis([ymin ymax])
%     end
%     save_and_close_fig(figure_dir,['SEquence_only_Low_FreqPow_' listsq_files{file}(1:end-13)])
%     
%     
%     %---Low Freq Coherence---%
%     yl = NaN(4,2);
%     figure
%     for it = 1:4
%         stat = lfp_phasecoherence(seq_saccade_aligneddata,find(itemnums == it));
%         low_freq_coherence_data{it,file} = stat;
%         subplot(2,2,it)
%         imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
%         line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
%         colorbar, axis xy
%         colormap jet
%         ylabel('Frequency (Hz)')
%         xlim([-0.75 1.0])
%         xlabel('Time from Saccade (sec)')
%         yl(it,:) = caxis;
%         title(['Item #' num2str(it)])
%     end
%     ymin = min(yl(:,1));
%     ymax = max(yl(:,2));
%     for it = 1:4
%         subplot(2,2,it)
%         caxis([ymin ymax])
%     end
%     
%     save_and_close_fig(figure_dir,['SEquence_only_Low_FreqCoh_' listsq_files{file}(1:end-13)])
%     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All Low Freq Power Spectrum---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfgfrq = [];
    cfgfrq.baseline = [-0.650 -0.4];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(4,2);
    
    figure
    for it = 1:4
        freq = lfp_powerspectrum(item_aligneddata,'low',find(itemnums == it),[-0.65 -0.4]);
        item_low_freq_data{it,file} = freq;
        
        subplot(2,2,it)
        hold on
        ft_singleplotTFR(cfgfrq, freq);
        hold off
        line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
        xlim([-0.75 0.75])
        xlabel('Time from Item On (sec)')
        yl(it,:) = caxis;
        title(['Item #' num2str(it)])
    end
    
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for it = 1:4
        subplot(2,2,it)
        caxis([ymin ymax])
    end
    save_and_close_fig(figure_dir,['SEquence_onlyItem_Low_FreqPow_' listsq_files{file}(1:end-13)])
    
    
    %---Low Freq Coherence---%
    yl = NaN(4,2);
    figure
    for it = 1:4
        stat = lfp_phasecoherence(item_aligneddata,find(itemnums == it));
        item_low_freq_coherence_data{it,file} = stat;
        subplot(2,2,it)
        imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
        line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
        colorbar, axis xy
        colormap jet
        ylabel('Frequency (Hz)')
        xlim([-0.75 1.0])
        xlabel('Time from Item On (sec)')
        yl(it,:) = caxis;
        title(['Item #' num2str(it)])
    end
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for it = 1:4
        subplot(2,2,it)
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['SEquence_onlyItem_Low_FreqCoh_' listsq_files{file}(1:end-13)])
    
    
    
end
save('SequenceOnly_across_tasks')
emailme('Done with Across task data analysis')
%%

all_low_freq_data = cell(1,4);
all_low_freq_coherence_data = cell(1,4);
total_sessions = 0;
for file = 1:length(listsq_files);
    if ~isempty(low_freq_data{file})
        for it = 1:4
            if file == 1
                all_low_freq_data{it}= low_freq_data{it,file}.powspctrm;
                all_low_freq_coherence_data{it} = low_freq_coherence_data{it,file}.cohspctrm;
            else
                all_low_freq_data{it}= all_low_freq_data{it}+low_freq_data{it,file}.powspctrm;
                all_low_freq_coherence_data{it} =  all_low_freq_coherence_data{it}+low_freq_coherence_data{it,file}.cohspctrm;
            end
        end
        total_sessions = total_sessions+1;
    end
end

%%
cfgfrq = [];
cfgfrq.baseline = [-0.650 -0.4];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

yl = NaN(4,2);

figure
for it = 1:4
    freq = low_freq_data{1};
    freq.powspctrm =  all_low_freq_data{it}/total_sessions;
    
    subplot(2,2,it)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    hold off
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(it,:) = caxis;
    title(['Item #' num2str(it)])
end

% ymin = min(yl(:,1));
% ymax = max(yl(:,2));
% for it = 1:4
%     subplot(2,2,it)
%     caxis([ymin ymax])
% end

%%
yl = NaN(4,2);
figure
for it = 1:4
    stat = low_freq_coherence_data{1};
    stat.cohspctrm = all_low_freq_coherence_data{it}/total_sessions;
    subplot(2,2,it)
    imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(it,:) = caxis;
    title(['Item #' num2str(it)])
end
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for it = 1:4
    subplot(2,2,it)
    caxis([ymin ymax])
end