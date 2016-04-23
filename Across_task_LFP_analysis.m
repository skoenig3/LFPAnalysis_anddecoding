%written 11/10 && 11/11 2015 by Seth Konig using/combining/modify code
%previously written to anlayze task speartely. Code focuses on eye movments
%or lack there off!!!!

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Across Tasks Figures\'; %where to save figures

%sessions with recordings on the same day
% listsq_files = {'PW140729_3-sorted.nex','PW140730_3-sorted.nex','PW140806_3-sorted.nex',...
%     'PW140829_3-sorted.nex','PW141007_3-sorted.nex','PW141008_3-sorted.nex',...
%     'PW141009_3-sorted.nex','PW141015_3-sorted.nex','PW141024_3-sorted.nex',...
%     'PW141028_3-sorted.nex','PW150205_3-sorted.nex'};
% 
% cvtnew_files = {'PW140729_4-sorted.nex','PW140730_4-sorted.nex','PW140806_4-sorted.nex',...
%     'PW140829_4-sorted.nex','PW141007_4-sorted.nex','PW141008_4-sorted.nex',...
%     'PW141009_4-sorted.nex','PW141015_4-sorted.nex','PW141024_4-sorted.nex',...
%     'PW141028_4-sorted.nex','PW150205_4-sorted.nex'};

listsq_files = {'TO160107_3-sorted.nex'};
cvtnew_files = {'TO160107_4-sorted.nex'};

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

event_time_lock_data = cell(6,length(listsq_files));
low_freq_data = cell(6,length(listsq_files));
high_freq_data =  cell(6,length(listsq_files));
low_freq_coherence_data = cell(6,length(listsq_files));
high_freq_coherence_data = cell(6,length(listsq_files));

for file = 1:length(listsq_files)
    
    %double check files are from the same day
    date1 = str2double(listsq_files{file}(3:8));
    date2 = str2double(cvtnew_files{file}(3:8));
    
    if date1 ~= date2
        disp(num2str(file))
    end
    
    disp(['Running File #' num2str(file) ': ' listsq_files{file}(1:end-13) ])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---List portion of ListSQ---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    l_trials(file) = length(trials_with_img_pairs)/2;
    if length(trials_with_img_pairs)/2 < 48
        disp(['Skipping File #' num2str(file) ': ' listsq_files{file}(1:end-11)])
        continue %need more data than this
    end
    
    fixationstats = fixationstats(trials_with_img_pairs);
    cfg.trl = cfg.trl(trials_with_img_pairs);
    which_img = which_img(trials_with_img_pairs);
    novel_vs_repeat = novel_vs_repeat(trials_with_img_pairs);
    num_trials = length(which_img);
    
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
    fixation_durations = NaN(num_trials,75); %keep track of fixation duration following saccade because we may want to sort things out
    for t = 1:length(fixationstats);
        if any(cfg.trl(t).allval == imgon_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code);%to ignore 1st fixation
            imgoff =  cfg.trl(t).alltim(cfg.trl(t).allval == imgoff_code);%to ignore last fixation if cutoff
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
                timeoutside = find(timeoutside >= maxoutside); %looked away a lot, again trying to cut out more data but with reason
                if ~isempty(timeoutside)
                    timeoutside = timeoutside(1);
                    toolate = find(saccadetimes(1,:) > timeoutside);
                    saccadetimes(:,toolate) = NaN;
                    toolate = find(fixationtimes(1,:) > timeoutside);
                    fixationtimes(:,toolate) = NaN;
                end
                
                toolate = find(saccadetimes(2,:)) > imgoff-trial_start;
                saccadetimes(:,toolate) = NaN;
                toolate = find(fixationtimes(2,:)) > imgoff-trial_start;
                fiationtimes(:,toolate) = NaN;
                
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
                fixation_durations(t,1:length(goodsacs)) = fixdur(goodsacs);
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
    
    nov_saccade_aligned = NaN(4*num_trials,3);
    rep_saccade_aligned = NaN(4*num_trials,3);
    novfixation_duration = NaN(4*num_trials);
    repfixation_duration = NaN(4*num_trials);
    
    nov_index = 1;
    rep_index = 1;
    for trl = 1:length(cfg.trl);
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        
        sacs = saccade_times(trl,:);
        maxind = find(~isnan(sacs));
        for sac = 1:max(maxind)
            if ~isnan(saccade_times(trl,sac))
                if novel_vs_repeat(trl) == 1 %novel image
                    nov_saccade_aligned(nov_index,:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                        saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                    novfixation_duration(nov_index) = fixation_durations(trl,sac);
                    nov_index = nov_index+1;
                else %repeat image
                    rep_saccade_aligned(rep_index,:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                        saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                    repfixation_duration(nov_index) = fixation_durations(trl,sac);
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
    cfg.dftfreq       = Fline;
    cfg.padding       = 1;
    
    %LFP data aligned saccades on novel images
    cfgnovsac = cfg;
    cfgnovsac.trl = nov_saccade_aligned;
    novsac_aligneddata = ft_preprocessing(cfgnovsac);
    
    %LFP data aligned saccades on repeat images
    cfgrepsac = cfg;
    cfgrepsac.trl = rep_saccade_aligned;
    repsac_aligneddata = ft_preprocessing(cfgrepsac);
    
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Covert Attention Task---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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
    
    %fixation aligned LFP data
    cvtfixcfg = cfg;
    cvtfixcfg.trl = fixation_aligned;
    cvtfix_aligneddata = ft_preprocessing(cvtfixcfg);
    
    %dot aligned LFP data
    dotcfg = cfg;
    dotcfg.trl = dot_aligned;
    dot_aligneddata = ft_preprocessing(dotcfg);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All LFP "Waveforms"---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yl = NaN(6,2);
    figure
    
    %Saccades during novel images during ListSQ List trials
    timelock = ft_timelockanalysis(cfgnovsac,novsac_aligneddata);
    event_time_lock_data{1,file} = timelock;
    subplot(2,3,1)
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    ylabel('Voltage (uV)')
    yl(1,:) = ylim;
    title('List Novel Images')
    
    %Saccades during repeat images during ListSQ List trials
    timelock = ft_timelockanalysis(cfgrepsac,repsac_aligneddata);
    event_time_lock_data{4,file} = timelock;
    subplot(2,3,4)
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.75])
    line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
    xlabel('Time from Saccade (sec)')
    ylabel('Voltage (uV)')
    yl(4,:) = ylim;
    title('List Repeat Images')
    
    %Saccades during ListSQ sequence trials
    timelock = ft_timelockanalysis(seqcfg,seq_saccade_aligneddata);
    event_time_lock_data{2,file} = timelock;
    subplot(2,3,2)
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    ylabel('Voltage (uV)')
    yl(2,:) = ylim;
    title('Sequence Trials')
    
    %Item Aligned during ListSQ sequence trials
    timelock = ft_timelockanalysis(itmcfg,item_aligneddata);
    event_time_lock_data{5,file} = timelock;
    subplot(2,3,5)
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.75])
    xlabel('Time from Item On (sec)')
    ylabel('Voltage (uV)')
    yl(5,:) = ylim;
    title('Sequence Trials')
    
    %CVT new Aligned to fixation on Cross hair
    timelock = ft_timelockanalysis(cvtfixcfg,cvtfix_aligneddata);
    event_time_lock_data{3,file} = timelock;
    subplot(2,3,3)
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.75])
    xlabel('Time from Fixation (sec)')
    ylabel('Voltage (uV)')
    yl(3,:) = ylim;
    title('CVT Trials')
    
    %CVT new Aligned to Dot on
    timelock = ft_timelockanalysis(dotcfg,dot_aligneddata);
    event_time_lock_data{6,file} = timelock;
    subplot(2,3,6)
    plot(timelock.time,timelock.avg')
    xlim([-0.75 0.75])
    xlabel('Time from Dot on (sec)')
    ylabel('Voltage (uV)')
    yl(6,:) = ylim;
    title('CVT Trials')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:6
        subplot(2,3,sb)
        hold on
        line([0 0],[ymin ymax],'Color','k','LineStyle','--')
        hold off
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        grid(gca,'minor')
        ylim([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['AcrossTaskLFPs-' listsq_files{file}(1:end-13)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All Low Freq Power Spectrum---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(6,2);
    
    figure
    
    %Saccades during novel images during ListSQ List trials
    freq = lfp_powerspectrum(novsac_aligneddata,'low','all');
    low_freq_data{1,file} = freq;
    subplot(2,3,1)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    %Saccades during repeat images during ListSQ List trials
    freq = lfp_powerspectrum(repsac_aligneddata,'low','all');
    low_freq_data{4,file} = freq;
    subplot(2,3,4)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(4,:) = caxis;
    title('List Repeat Images')
    
    %Saccades during ListSQ sequence trials
    freq = lfp_powerspectrum(seq_saccade_aligneddata,'low','all');
    low_freq_data{2,file} = freq;
    subplot(2,3,2)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    yl(4,:) = caxis;
    xlabel('Time from Saccade (sec)')
    title('Sequence Trials')
    
    %Item Aligned during ListSQ sequence trials
    freq = lfp_powerspectrum(item_aligneddata,'low','all');
    low_freq_data{5,file} = freq;
    subplot(2,3,5)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time Item On (sec)')
    yl(4,:) = caxis;
    title('Sequence Trials')
    
    %CVT new Aligned to fixation on Cross hair
    freq = lfp_powerspectrum(cvtfix_aligneddata,'low','all');
    low_freq_data{3,file} = freq;
    subplot(2,3,3)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Fixation (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %CVT new Aligned to Dot on
    freq = lfp_powerspectrum(dot_aligneddata,'low','all');
    low_freq_data{6,file} = freq;
    subplot(2,3,6)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Dot on (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:6
        subplot(2,3,sb)
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['AcrossTask-LowFreqPowerSpectrum-' listsq_files{file}(1:end-13)]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All High Freq Power Spectrum---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(6,2);
    
    figure
    
    %Saccades during novel images during ListSQ List trials
    freq = lfp_powerspectrum(novsac_aligneddata,'high','all');
    high_freq_data{1,file} = freq;
    subplot(2,3,1)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    %Saccades during repeat images during ListSQ List trials
    freq = lfp_powerspectrum(repsac_aligneddata,'high','all');
    high_freq_data{4,file} = freq;
    subplot(2,3,4)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(4,:) = caxis;
    title('List Repeat Images')
    
    %Saccades during ListSQ sequence trials
    freq = lfp_powerspectrum(seq_saccade_aligneddata,'high','all');
    high_freq_data{2,file} = freq;
    subplot(2,3,2)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    yl(4,:) = caxis;
    xlabel('Time from Saccade (sec)')
    title('Sequence Trials')
    
    %Item Aligned during ListSQ sequence trials
    freq = lfp_powerspectrum(item_aligneddata,'high','all');
    high_freq_data{5,file} = freq;
    subplot(2,3,5)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time Item On (sec)')
    yl(4,:) = caxis;
    title('Sequence Trials')
    
    %CVT new Aligned to fixation on Cross hair
    freq = lfp_powerspectrum(cvtfix_aligneddata,'high','all');
    high_freq_data{3,file} = freq;
    subplot(2,3,3)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Fixation (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %CVT new Aligned to Dot on
    freq = lfp_powerspectrum(dot_aligneddata,'high','all');
    high_freq_data{6,file} = freq;
    subplot(2,3,6)
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 0.75])
    xlabel('Time from Dot on (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:6
        subplot(2,3,sb)
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['AcrossTask-highFreqPowerSpectrum-' listsq_files{file}(1:end-13)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All Low Freq Coherence---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure
    
    %Saccades during novel images during ListSQ List trials
    stat = lfp_phasecoherence(novsac_aligneddata,'all');
    low_freq_coherence_data{1,file} = stat;
    subplot(2,3,1)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(4,:) = caxis;
    title('List Novel Images')
    
    %Saccades during repeat images during ListSQ List trials
    stat = lfp_phasecoherence(repsac_aligneddata,'all');
    low_freq_coherence_data{4,file} = stat;
    subplot(2,3,4)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(4,:) = caxis;
    title('List Repeat Images')
    
    %Saccades during ListSQ sequence trials
    stat = lfp_phasecoherence(seq_saccade_aligneddata,'all');
    low_freq_coherence_data{2,file} = stat;
    subplot(2,3,2)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    yl(4,:) = caxis;
    xlabel('Time from Saccade (sec)')
    title('Sequence Trials')
    
    %Item Aligned during ListSQ sequence trials
    stat = lfp_phasecoherence(item_aligneddata,'all');
    low_freq_coherence_data{5,file} = stat;
    subplot(2,3,5)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time Item On (sec)')
    yl(4,:) = caxis;
    title('Sequence Trials')
    
    %CVT new Aligned to fixation on Cross hair
    stat = lfp_phasecoherence(cvtfix_aligneddata,'all');
    low_freq_coherence_data{3,file} = stat;
    subplot(2,3,3)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Fixation (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %CVT new Aligned to Dot on
    stat = lfp_phasecoherence(dot_aligneddata,'all');
    low_freq_coherence_data{6,file} = stat;
    subplot(2,3,6)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Dot on (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:6
        subplot(2,3,sb)
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['AcrossTask-LowFreqCoherence-' listsq_files{file}(1:end-13)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All High Freq Coherence---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    yl = NaN(6,2);
    
    figure
    
    %Saccades during novel images during ListSQ List trials
    stat = lfp_phasecoherence2High(novsac_aligneddata,'all');
    high_freq_coherence_data{1,file} = stat;
    subplot(2,3,1)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    %Saccades during repeat images during ListSQ List trials
    stat = lfp_phasecoherence2High(repsac_aligneddata,'all');
    high_freq_coherence_data{4,file} = stat;
    subplot(2,3,4)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(4,:) = caxis;
    title('List Repeat Images')
    
    %Saccades during ListSQ sequence trials
    stat = lfp_phasecoherence2High(seq_saccade_aligneddata,'all');
    high_freq_coherence_data{2,file} = stat;
    subplot(2,3,2)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    yl(4,:) = caxis;
    xlabel('Time from Saccade (sec)')
    title('Sequence Trials')
    
    %Item Aligned during ListSQ sequence trials
    stat = lfp_phasecoherence2High(item_aligneddata,'all');
    high_freq_coherence_data{5,file} = stat;
    subplot(2,3,5)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time Item On (sec)')
    yl(4,:) = caxis;
    title('Sequence Trials')
    
    %CVT new Aligned to fixation on Cross hair
    stat = lfp_phasecoherence2High(cvtfix_aligneddata,'all');
    high_freq_coherence_data{3,file} = stat;
    subplot(2,3,3)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Fixation (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %CVT new Aligned to Dot on
    stat = lfp_phasecoherence2High(dot_aligneddata,'all');
    high_freq_coherence_data{6,file} = stat;
    subplot(2,3,6)
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Dot on (sec)')
    yl(4,:) = caxis;
    title('CVT Trials')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:6
        subplot(2,3,sb)
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['AcrossTask-highFreqCoherence-' listsq_files{file}(1:end-13)]);
    
end
save('Alldata_across_tasks')
emailme('Done with Across task data analysis')
%%
all_event_time_lock_data = cell(1,6);
all_low_freq_data = cell(1,6);
all_high_freq_data = cell(1,6);
all_low_freq_coherence_data = cell(1,6);
all_high_freq_coherence_data = cell(1,6);

plotsnums  = [1 4 2 5 3 6];
all_titles =  {'List Novel Images','List Repeat Images','Sequence Trials',...
    'Sequence Trials','CVT Trials','CVT Trials'};
all_xlabels = {'Time from Saccade (sec)','Time from Saccade (sec)','Time from Saccade (sec)',...
    'Time Item On (sec)','Time from Fixation (sec)','Time from Dot on (sec)'};

total_sessions = 0;
for file = 1:length(listsq_files)
    if ~isempty(event_time_lock_data{1,file})
        for sb = 1:6
            if file == 1
                all_event_time_lock_data{sb} = event_time_lock_data{sb,file}.avg;
                all_low_freq_data{sb} = low_freq_data{sb,file}.powspctrm;
                all_high_freq_data{sb} = high_freq_data{sb,file}.powspctrm;
                all_low_freq_coherence_data{sb} = low_freq_coherence_data{sb,file}.cohspctrm;
                all_high_freq_coherence_data{sb} = high_freq_coherence_data{sb,file}.cohspctrm;
            else
                all_event_time_lock_data{sb} =  all_event_time_lock_data{sb}+event_time_lock_data{sb,file}.avg;
                all_low_freq_data{sb} = all_low_freq_data{sb}+low_freq_data{sb,file}.powspctrm;
                all_high_freq_data{sb} = all_high_freq_data{sb}+high_freq_data{sb,file}.powspctrm;
                all_low_freq_coherence_data{sb} = all_low_freq_coherence_data{sb}+low_freq_coherence_data{sb,file}.cohspctrm;
                all_high_freq_coherence_data{sb} = all_high_freq_coherence_data{sb}+high_freq_coherence_data{sb,file}.cohspctrm;
            end
        end
        total_sessions = total_sessions + 1;
    end
end
%%
%---Plot All session LFP waveforms---%

yl = NaN(6,2);
figure

for sb = 1:6
    timelock = event_time_lock_data{1,1};
    timelock.avg =  all_event_time_lock_data{plotsnums(sb)}/total_sessions;
    subplot(2,3,plotsnums(sb))
    plot(timelock.time,timelock.avg')
    xlim([-0.75 1.5])
    xlabel(all_xlabels{sb})
    ylabel('Voltage (uV)')
    yl(sb,:) = ylim;
    title(all_titles{sb})
end

%rescale and add line plot
ymin = min(yl(:,1));
ymax = 0.75*max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    line([0 0],[ymin ymax],'Color','k','LineStyle','--')
    hold off
    set(gca,'XMinorTick','on','YMinorTick','on')
    grid on
    grid(gca,'minor')
    ylim([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_AcrossTaskLFPs-' listsq_files{file}(1:end-13)]);


%---Plot All session Low Frequency Power Spectrum---%


cfgfrq = [];
cfgfrq.baseline = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

yl = NaN(6,2);


figure
for sb = 1:6
    freq = low_freq_data{1,1};
    freq.powspctrm = all_low_freq_data{plotsnums(sb)}/total_sessions;
    subplot(2,3,plotsnums(sb))
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 1.5])
    xlabel(all_xlabels{sb})
    yl(sb,:) = caxis;
    title(all_titles{sb})
end

%rescale and add line plot

ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    plot([0 0],[3.5 30.5],'w-')
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_AcrossTask_LowFreqPower-' listsq_files{file}(1:end-13)]);

%---Plot All session high Frequency Power Spectrum---%


cfgfrq = [];
cfgfrq.baseline = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

yl = NaN(6,2);


figure
for sb = 1:6
    freq = high_freq_data{1,1};
    freq.powspctrm = all_high_freq_data{plotsnums(sb)}/total_sessions;
    subplot(2,3,plotsnums(sb))
    ft_singleplotTFR(cfgfrq, freq);
    xlim([-0.75 1.5])
    xlabel(all_xlabels{sb})
    yl(sb,:) = caxis;
    title(all_titles{sb})
end

%rescale and add line plot

ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    plot([0 0],[29 181],'w-')
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_AcrossTask_highFreqPower-' listsq_files{file}(1:end-13)]);


%---Plot All session Low Freq Coherence Spectrum---%

yl = NaN(6,2);

figure
for sb = 1:6
    stat = low_freq_coherence_data{1,1};
    stat.cohspctrm = all_low_freq_coherence_data{plotsnums(sb)}/total_sessions;
    subplot(2,3,plotsnums(sb))
    imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    xlim([-0.75 1.5])
    xlabel(all_xlabels{sb})
    yl(sb,:) = caxis;
    title(all_titles{sb})
end

%rescale and add line plot

ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_AcrossTask_LowFreqCoherence-' listsq_files{file}(1:end-13)]);

%---Plot All session high Freq Coherence Spectrum---%

yl = NaN(6,2);

figure
for sb = 1:6
    stat = high_freq_coherence_data{1,1};
    stat.cohspctrm = all_high_freq_coherence_data{plotsnums(sb)}/total_sessions;
    subplot(2,3,plotsnums(sb))
    imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    xlim([-0.75 1.5])
    xlabel(all_xlabels{sb})
    yl(sb,:) = caxis;
    title(all_titles{sb})
end

%rescale and add line plot

ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_AcrossTask_highFreqCoherence-' listsq_files{file}(1:end-13)]);
