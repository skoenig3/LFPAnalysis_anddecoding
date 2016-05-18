function get_saccade_aligned_LFPs_across_tasks(session_data,data_dir,eye_data_dir)
% written by Seth Konig 4/28/16. Code based on previous LFP analysis code
% and updated to handle new data structure. Code imports LFP data for LFPs
% aligned to saccades in the sequence and list task and fixations and
% dotonset in the coverte attention task. Code also removes any line noise.
% Lastly, code imports LFPs during baseline condition (ITI 250-750 ms
% period) so can do baseline power subtraction. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----Check to make sure there are 2 sessions with tasks with at least 100 trials---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
miniumum_num_trials = 100;

task = 'ListSQ';
[task_file] = get_task_data(session_data,task);
if isempty(task_file)
    return
else
    load([eye_data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg')
    total_rewarded_trials = 0;
    for trl = 1:length(cfg.trl)
        if any(cfg.trl(trl).allval == 3)
            total_rewarded_trials =   total_rewarded_trials+1;
        end
    end
    if total_rewarded_trials < miniumum_num_trials
        return
    end
end
task = 'cvtnew';
[task_file] = get_task_data(session_data,task);
if isempty(task_file)
    return
else
    load([eye_data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg')
    total_rewarded_trials = 0;
    for trl = 1:length(cfg.trl)
        if any(cfg.trl(trl).allval == 3)
            total_rewarded_trials =   total_rewarded_trials+1;
        end
    end
    if total_rewarded_trials < miniumum_num_trials
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set Important Variables & Values---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Important parametrs/values for ListSQ List---%
img_on_code = 23;
img_off_code = 24;
fixspoton = 35;
imageX = 800; %horizontal size of the screen
imageY = 600; %horizontal size of the screen
minfixdur = 100;%how long of fixation is required after a saccade
%1000 gives them a reasonable buffer
twin = 500; %how long to ignore eye data after image onset


%---Important parameters/values for ListSQ Sequence---%
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
reward_code = 3;


%---Important parameters/values for CVTNew---%
dot_on_code = 25;
dot_clrchng_code = 27;
bar_code_response = 4; %monkey made a move
%time is + 300 ms
short = [700 1133]; %short duration trials
mediumm = [1134 1567];
long = [1568 2500];%cap should be 2000 but cortex can have some lag

%---Important parameters/values across all tasks---%
trial_start_code = 15;
fixation_on_cross = 8;
reward_code = 3;
Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity
Fs = 1000;
buffer1 = 2048;
buffer2 = 2048;
baseline_period = [250 750];%into ITI period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%---List portion of ListSQ---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'ListSQ';
[task_file,item_file,cnd_file,~,~,~,~,~,listsq_lfp_quality] = get_task_data(session_data,task);

disp(['Importing ListSQ Data for ' task_file(1:8) ])
load([eye_data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','fixationstats','hdr')

%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

if strcmp('TO160107_3-sorted.nex',task_file)
    %error in cfg.trl.cnd for 1 condition so need to remove these
    cfg.trl(296).allval = NaN;  
    cfg.trl(428).allval = NaN; 
    cfg.trl(429).allval = NaN; 
end                 
%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---process eye data and fixation durations---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_trials = length(cfg.trl);

%---Storage variables for List Images---%
list_saccade_times = NaN(num_trials,75); %start time of saccades
fixation_durations = NaN(num_trials,75); %keep track of fixation duration following saccade because we may want to sort things out
novel_repeat = NaN(num_trials,75);
list_startind = NaN(1,num_trials);%trial start index in all data

%---Storage Variables for Sequence---%
time_to_fixation = NaN(num_trials,4); %reaction times
sequence_saccade_times = NaN(num_trials,4); %start time of saccades
sequence_startind = NaN(1,num_trials);%trial start index in all data

%---Storage Variables for ITI Baseline--%
baseline_start_end = NaN(num_trials,3);

for t = 1:num_trials
    if t == 54
        continue
    end
        
    if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd(1)-1000) > sequence_items(end) %only want image trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        list_startind(t) = cfg.trl(t).begsmpind;%trial start index in all data
        imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code);%to ignore 1st fixation
        imgon = imgon(1);
        imgoff =  cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code);%to ignore last fixation if cutoff
        
        % if monkey isn't paying attention data isn't probably
        % worth much plus have to cut off somewhere
        if imgoff-imgon > 1.5*imgdur-1
            imgoff = imgon+1.5*imgdur-1;
        end
        
        saccadetimes = fixationstats{t}.saccadetimes;
        fixationtimes = fixationstats{t}.fixationtimes;
        if ~isempty(fixationtimes)
            %ignore saccades before the first 250 ms of image on to
            %remove image onset confounds
            
            
            invalid = find(saccadetimes(2,:)) > imgoff-trial_start;
            saccadetimes(:,invalid) = NaN;
            invalid = find(fixationtimes(2,:)) > imgoff-trial_start;
            fixationtimes(:,invalid) = NaN;
            
            invalid = find(saccadetimes(1,:) <= imgon-trial_start+twin);
            saccadetimes(:,invalid) = NaN;
            invalid = find(fixationtimes(1,:) <= imgon-trial_start+twin);
            fixationtimes(:,invalid) = NaN;
            
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
            list_saccade_times(t,1:length(goodsacs)) = saccadetimes(1,goodsacs);
            fixation_durations(t,1:length(goodsacs)) = fixdur(goodsacs);
            novel_repeat(t,1:length(goodsacs)) = novel_vs_repeat(img_cnd == cfg.trl(t).cnd(1));
        end
    elseif any(cfg.trl(t).allval == reward_code) %rewarded sequence trials
        sequence_startind(t) = cfg.trl(t).begsmpind;%trial start index in all data

        locs = sequence_locations{find(itmlist(cfg.trl(t).cnd-1000) == sequence_items)};
        %convert to DVA for this analysis
        locs(1,:) = (locs(1,:)-400)/24;
        locs(2,:) = (locs(2,:)-300)/24;
        fixationstats{t}.fixations(1,:) = (fixationstats{t}.fixations(1,:)-400)/24;
        fixationstats{t}.fixations(2,:) = (fixationstats{t}.fixations(2,:)-300)/24;
        
        events = cfg.trl(t).allval;
        events(events == 100)= 0;
        events(1) = 100;%eye data starts for recording right away
        event_times = cfg.trl(t).alltim;
        
        trialdata = analyze_sequence_trial(fixationstats{t},locs,fixwin,...
            events,event_times);
        
        time_to_fixation(t,:) = trialdata.t2f;
        
        fixation_numbers = trialdata.fixationnums; %fixation number for each item
        fixationtimes = fixationstats{t}.fixationtimes;
        saccadetimes = fixationstats{t}.saccadetimes;
        for item = 1:4
            if ~isnan(fixation_numbers(item))
                saccadeind = find(saccadetimes(2,:)+1 ==  fixationtimes(1,fixation_numbers(item)));
                if ~isempty(saccadeind)
                    sequence_saccade_times(t,item) = saccadetimes(1,saccadeind);
                end
            end
        end
    end
    
    trial_start = cfg.trl(t).begsmpind;%trial start index in all data
    baseline_start_end(t,:) = [trial_start+baseline_period(1)-1-buffer1 ...
        trial_start+baseline_period(2)-1+buffer2 -buffer1];
end
%---remove excess NaNs---%
baseline_start_end = laundry(baseline_start_end);
list_saccade_times = laundry(list_saccade_times);
novel_repeat = laundry(novel_repeat);
time_to_fixation = laundry(time_to_fixation);
sequence_saccade_times = laundry(sequence_saccade_times);
fixation_durations = laundry(fixation_durations);
list_startind = laundry(list_startind);
sequence_startind = laundry(sequence_startind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---convert Data into Field Trip Friendly Format---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% want following organization [event start index-buffer1, event end index+buffer2, -buffer1]

%---For List Trials---%
num_fix_list = sum(sum(~isnan(list_saccade_times)));
list_saccade_aligned = NaN(num_fix_list,3);
stretched_fixation_durations = NaN(1,num_fix_list);
stretched_novel_vs_repeat = NaN(1,num_fix_list);

index = 1;
for trl = 1:size(list_saccade_times,1)
    sacs = list_saccade_times(trl,:);
    fixdur = fixation_durations(trl,:);
    fixdur(isnan(sacs)) = [];
    sacs(isnan(sacs)) = [];
    for sacind = 1:length(sacs)
        if novel_vs_repeat(trl) == 1 %novel image
            stretched_novel_vs_repeat(index) = 1;
        else %repeat image
            stretched_novel_vs_repeat(index) = 2;
        end
        stretched_fixation_durations(index) = fixdur(sacind);
        list_saccade_aligned(index,:) = [sacs(sacind)-1-buffer1+list_startind(trl) ...
            sacs(sacind)-1+buffer2+list_startind(trl) -buffer1];
        index = index+1;
    end
end
%remove trials if data extends past end of file
too_late = find(list_saccade_aligned(:,2) > hdr.nSamples);
list_saccade_aligned(too_late,:) = [];
stretched_fixation_durations(too_late) = [];
stretched_novel_vs_repeat(too_late) = []; 
too_late = find(baseline_start_end(:,2) > hdr.nSamples);
baseline_start_end(too_late,:) = [];


%---For Sequence Trials---%
num_fix_sequence = sum(sum(~isnan(sequence_saccade_times)));
sequence_saccade_aligned = NaN(num_fix_sequence,3);
stretched_time_to_fixation = NaN(1,num_fix_sequence);

index = 1;
for trl = 1:size(sequence_saccade_times)
    sacs = sequence_saccade_times(trl,:);
    t2f = time_to_fixation(trl,:);
    t2f(isnan(sacs)) = [];
    sacs(isnan(sacs)) = [];
    for sacind = 1:length(sacs)
        stretched_time_to_fixation(index) = t2f(sacind);
        sequence_saccade_aligned(index,:) = [sacs(sacind)-1-buffer1+sequence_startind(trl) ...
            sacs(sacind)-1+buffer2+sequence_startind(trl) -buffer1];
        index = index+1;
    end
end
%remove trials if data extends past end of file
too_late = find(sequence_saccade_aligned(:,2) > hdr.nSamples);
sequence_saccade_aligned(too_late,:) = [];
stretched_time_to_fixation(too_late) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import LFP data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%also filters out line noise

%code to find spike/LFP channels
LFPchannels = find_desired_channels(cfg,'LFP');

bad_channels = [];
for channel = 1:4
    if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
        if  listsq_lfp_quality(channel) == 0; %if it is bad
            bad_channels = [bad_channels channel];
        end
    end
end
LFPchannels(bad_channels) = [];

% read the data from file and preprocess them
cfg.channel       = cfg.channel(LFPchannels)';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = Fline;
cfg.padding       = 1;
cfg.continuous    = 'yes';

%LFP data aligned saccades for list
disp(['Importing List LFP Data for ' task_file(1:8) ])
cfg_list_sac = cfg;
cfg_list_sac.trl = list_saccade_aligned;
list_saccade_aligneddata = ft_preprocessing(cfg_list_sac);

%LFP data aligned saccades for sequence
disp(['Importing Sequence LFP Data for ' task_file(1:8) ])
cfg_seq_sac = cfg;
cfg_seq_sac.trl = sequence_saccade_aligned;
sequence_saccade_aligneddata = ft_preprocessing(cfg_seq_sac);

%LFP data aligned to ITI period for baseline
disp(['Importing ListSQ ITI Data for ' task_file(1:8) ])
cfg_listsq_ITI = cfg;
cfg_listsq_ITI.trl =  baseline_start_end;
listsq_ITI_aligneddata = ft_preprocessing(cfg_listsq_ITI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---remove any trials with NaNs, code can't process them well--%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for list saccades
remove_trials = [];
for t = 1:length(list_saccade_aligneddata.trial)
    if sum(sum(isnan(list_saccade_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
list_saccade_aligneddata.time(remove_trials) = [];
list_saccade_aligneddata.trial(remove_trials) = [];
list_saccade_aligneddata.sampleinfo(remove_trials,:) = [];

%for sequence saccades
remove_trials = [];
for t = 1:length(sequence_saccade_aligneddata.trial)
    if sum(sum(isnan(sequence_saccade_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
sequence_saccade_aligneddata.time(remove_trials) = [];
sequence_saccade_aligneddata.trial(remove_trials) = [];
sequence_saccade_aligneddata.sampleinfo(remove_trials,:) = [];

%for listsq ITI period
remove_trials = [];
for t = 1:length(listsq_ITI_aligneddata.trial)
    if sum(sum(isnan(listsq_ITI_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
listsq_ITI_aligneddata.time(remove_trials) = [];
listsq_ITI_aligneddata.trial(remove_trials) = [];
listsq_ITI_aligneddata.sampleinfo(remove_trials,:) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%---Covert Attention Task---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import & Preprocess Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'cvtnew';
[task_file,~,~,~,~,~,~,~,cvtnew_lfp_quality] = get_task_data(session_data,task);
load([eye_data_dir task_file(1:end-11) '-preprocessed'],'cfg','meta','hdr');

%---preallocate space and parallel structure of cfg--%
num_trials = length(cfg.trl);
baseline_start_end = NaN(num_trials,3);
successful_trials = NaN(1,length(cfg.trl));

for t = 1:length(cfg.trl);
    if sum(cfg.trl(t).allval == 3) > 0;
        successful_trials(t) = t;
    end
    trial_start = cfg.trl(t).begsmpind;%trial start index in all data
    baseline_start_end(t,:) = [trial_start+baseline_period(1)-1-buffer1 ...
        trial_start+baseline_period(2)-1+buffer2 -buffer1];
end
successful_trials = laundry(successful_trials);
cfg.trl = cfg.trl(successful_trials);
cfg.trlold = cfg.trl; %store since may need if and we're going to write over

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---convert Data into Field Trip Friendly Format---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
%keep track of sequences and item number for later analysis and store
%by event/eye movement
num_trials = length(successful_trials);
dot_aligned = NaN(num_trials,3); %align to dot turning on
cvtnew_fixation_aligned = NaN(num_trials,3); %aligned to fixation on cross
for t = 1:num_trials
    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
    startind = cfg.trl(t).begsmpind;%trial start index in all data

    all_event_times = NaN(1,2);
    all_event_times(1) = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_on_cross)-trial_start;
    all_event_times(2) = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start;
    
    cvtnew_fixation_aligned(t,:) = [all_event_times(1)-1-buffer1+startind ...
        all_event_times(1)-1+buffer2+startind -buffer1];
    dot_aligned(t,:) = [all_event_times(2)-1-buffer1+startind ...
        all_event_times(2)-1+buffer2+startind -buffer1];
end


%clean up/remove any NaNs (occur when saccades couldn't be detected)
dot_aligned = laundry(dot_aligned);
cvtnew_fixation_aligned = laundry(cvtnew_fixation_aligned);


%remove trials if data extends past end of file
too_late = find(dot_aligned(:,2) > hdr.nSamples);
dot_aligned(too_late,:) = [];
too_late = find(cvtnew_fixation_aligned(:,2) > hdr.nSamples);
cvtnew_fixation_aligned(too_late,:) = [];
too_late = find(baseline_start_end(:,2) > hdr.nSamples);
baseline_start_end(too_late,:) = [];

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

bad_channels = [];
for channel = 1:4
    if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
        if  cvtnew_lfp_quality(channel) == 0; %if it is bad
            bad_channels = [bad_channels channel];
        end
    end
end
LFPchannels(bad_channels) = [];

% read the data from file and preprocess them
cfg.channel       = cfg.channel(LFPchannels)';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = Fline;
cfg.padding       = 1;

%fixation aligned LFP data
disp(['Importing CVTNEW Fixation LFP Data for ' task_file(1:8) ])
cvtfixcfg = cfg;
cvtfixcfg.trl = cvtnew_fixation_aligned;
cvtfix_aligneddata = ft_preprocessing(cvtfixcfg);

%dot aligned LFP data
disp(['Importing CVTNEW Dot LFP Data for ' task_file(1:8) ])
dotcfg = cfg;
dotcfg.trl = dot_aligned;
dot_aligneddata = ft_preprocessing(dotcfg);

%LFP data aligned to ITI period for baseline
cfg_cvtnew_ITI = cfg;
cfg_cvtnew_ITI.trl =  baseline_start_end;
cvtnew_ITI_aligneddata = ft_preprocessing(cfg_cvtnew_ITI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---remove any trials with NaNs, code can't process them well--%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for cvtnew ITI period
remove_trials = [];
for t = 1:length(cvtnew_ITI_aligneddata.trial)
    if sum(sum(isnan(cvtnew_ITI_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
cvtnew_ITI_aligneddata.time(remove_trials) = [];
cvtnew_ITI_aligneddata.trial(remove_trials) = [];
cvtnew_ITI_aligneddata.sampleinfo(remove_trials,:) = [];

%for cvtnew fixation
remove_trials = [];
for t = 1:length(cvtfix_aligneddata.trial)
    if sum(sum(isnan(cvtfix_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
cvtfix_aligneddata.time(remove_trials) = [];
cvtfix_aligneddata.trial(remove_trials) = [];
cvtfix_aligneddata.sampleinfo(remove_trials,:) = [];

%for cvtnew dot aligned data
remove_trials = [];
for t = 1:length(cvtfix_aligneddata.trial)
    if sum(sum(isnan(dot_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
dot_aligneddata.time(remove_trials) = [];
dot_aligneddata.trial(remove_trials) = [];
dot_aligneddata.sampleinfo(remove_trials,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save the Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-preproccessed_saccade_aligned_LFPs.mat'],...
    'sequence_saccade_aligneddata','stretched_time_to_fixation',...
    'buffer1','buffer2','list_saccade_aligneddata','stretched_fixation_durations',...
    'stretched_novel_vs_repeat','listsq_lfp_quality','short_trials',...
    'long_trials','medium_trials','cvtfix_aligneddata','dot_aligneddata',...
     'listsq_ITI_aligneddata','cvtnew_ITI_aligneddata');

end