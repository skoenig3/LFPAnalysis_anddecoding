clar

task = 'ListSQ';
twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trialstart_code = 15;
smval = 50;%gaussian 1/2 width for smoothing
numshuffs = 1000; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;
min_blks = 2;

fc = 90;
fs = 1000;
[b,a] = butter(6,fc/(fs/2),'high');
[b2,a2] = butter(6,12/(fs/2),'low');


%all parameters
LFP_count = 0;%number of LFP/electrodes analyzed
which_monkey = [];


%for sequence trials
seq_HFO_trig_avg = [];
seq_high_HFO_trig_avg = [];
seq_high2_HFO_trig_avg = [];
seq_total_HFO_count = 0;
seq_time_distribution = zeros(1,12500);
seq_total_duration = 0; %total LFP duration analyzed
was_predicted = []; %1 for predicted, 0 for reactive
seq_co_occur_count = 0;
seq_all_trial_power = zeros(100,751);
all_trial_predict = zeros(1,2); %1: predicted count, 2 not predicted count
sequence_proximity = [];

%for image trials
image_HFO_trig_avg = [];
image_high_HFO_trig_avg = [];
image_high2_HFO_trig_avg = [];
image_total_HFO_count = zeros(1,2);
image_time_distribution = zeros(2,12500);
image_total_duration = 0; %total LFP duration analyzed
novel_vs_repeat = [];
image_co_occur_count = 0;
image_all_trial_power = zeros(100,751);
image_proximity = [];

for monkey = 1%2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        %         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess = 1:length(session_data)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        task = 'ListSQ';
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,...
            sorting_quality,waveform_count,lfp_quality,comments] = get_task_data(session_data{sess},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        
        %load trial data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg',...
            'hdr','fixationstats');
        disp(['Session #' num2str(sess)])
        num_trials = length(cfg.trl);
        
        LFPchannels = find_desired_channels(cfg,'LFP');
        %remove bad LFP channels
        bad_channels = [];
        for channel = 1:4
            if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                if  lfp_quality(channel) == 0; %if it is bad
                    bad_channels = [bad_channels channel];
                end
            end
        end
        try
            LFPchannels(bad_channels) = [];
        catch
            continue
        end
        LFP_count = LFP_count+length(LFPchannels);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Get successful trials Information by Task---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);
        
                        
        %set the image duration
        if str2double(task_file(3:8)) < 140805 %first sets done by PW before this data had 7 s images
            imgdur = 7000;
        else %rest of image presentations were 5 seconds
            imgdur = 5000;
        end
        imgdur = imgdur*1.5;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---process eye data locked to trial events---%%%
        num_trials = length(cfg.trl);
        fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
        reaction_times = NaN(length(fixationstats),4);%when did fixation on item start
        which_sequence = NaN(1,length(cfg.trl));
        
        ons_offs = NaN(num_trials,9);
        trialtype = NaN(1,num_trials);%1 for sequence, 3  for novel images, 4 for repeat images
        for trial = 1:num_trials
            if sum(cfg.trl(trial).allval == 3) >= 5; %in which sequence trials were rewarded
                which_sequence(trial) = find(sequence_items == itmlist(cfg.trl(trial).cnd-1000));
                locs = sequence_locations{which_sequence(trial)};
                trialtype(trial) = 1;
                
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
                    event_codes,event_times,1);
                
                fixation_numbers = trialdata.fixationnums; %fixation number for each item
                fixationtimes = fixationstats{trial}.fixationtimes;
                trialstart = cfg.trl(trial).alltim(1);
                reward = cfg.trl(trial).alltim(cfg.trl(trial).allval == 3)-trialstart;
                reward = reward(1);
                ons_offs(trial,9) = reward;
                for item = 1:4
                    
                    item_on = cfg.trl(trial).alltim(cfg.trl(trial).allval == item_event_codes(2*item-1)); %when item turned on
                    item_on = item_on-trialstart;
                    item_off = cfg.trl(trial).alltim(cfg.trl(trial).allval == item_event_codes(2*item));
                    item_off = item_off-trialstart;
                    ons_offs(trial,2*item-1) = item_on;
                    ons_offs(trial,2*item) = item_off;
                    
                    if ~isnan(fixation_numbers(item))
                        fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
                        reaction_times(trial,item) = fixation_start_time(trial,item)-item_on;
                    end
                end
                if any(reaction_times(trial,:) < predict_rt)
                    all_trial_predict(1) = all_trial_predict(1)+1; %had a predictive eye movement
                else 
                     all_trial_predict(2) = all_trial_predict(2)+1; %didnt' ahve a predictive eye movement
                end
            elseif any(cfg.trl(trial).allval == 23) && itmlist(cfg.trl(trial).cnd-1000) > 19 %successful image trial

                %image trial info
                img_index = find(cfg.trl(trial).cnd == img_cnd); %image index
                if any(isnan(which_img(img_index)))
                    continue
                end
                if novel_vs_repeat(img_index) == 1 %novel
                    trialtype(trial) = 3;
                else%repeat
                    trialtype(trial) = 4;
                end
                
                trialstart = cfg.trl(trial).alltim(1);
                crosson = cfg.trl(trial).alltim(cfg.trl(trial).allval == 35)-trialstart; %cross on
                imgon = cfg.trl(trial).alltim(cfg.trl(trial).allval == 23)-trialstart;%image turned on
                imgoff = cfg.trl(trial).alltim(cfg.trl(trial).allval == 24)-trialstart;%image turned on
                
                ons_offs(trial,1) = crosson;
                ons_offs(trial,2) = imgon;
                ons_offs(trial,3) = imgoff;
                
                ons_offs(trial,9) = imgon+imgdur; %cap time
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Calculate Firing Rate Locked to Eye Data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %also locked to "trial start" item 1 on or saccade to item 1
        trial_LFPs = cell(num_trials,length(LFPchannels));
        for trial = 2:num_trials
            if all(isnan(ons_offs(trial,:)))
                continue
            end
            if trialtype(trial) == 1 %sequence trial
                seq_total_duration = seq_total_duration+ons_offs(trial,9);%only anlzying to reward start
            elseif trialtype(trial) >= 3 %image trial
                image_total_duration = image_total_duration+ons_offs(trial,9);%only anlzying to 1.5*image duration afeter image is turned on
            end
            for chan = 1:length(LFPchannels)
                trial_LFPs{trial,chan} =  [data(LFPchannels(chan)).values{trial-1}(end-250:end) ...
                    data(LFPchannels(chan)).values{trial}];
                %add 250 ms from previous trial, will not anlayze but will
                %use as buffer
            end
        end
        
        %High Pass filter LFPs
        filtered_trial_LFPs = cell(size(trial_LFPs));
        for trial = 1:num_trials
            for chan = 1:size(trial_LFPs,2)
                if ~isempty(trial_LFPs{trial,chan})
                    filtered_trial_LFPs{trial,chan} = abs(filtfilt(b,a,trial_LFPs{trial,chan}));
                end
            end
        end
        
        %Low Pass filter rectified LFPs
        envelope_trial_LFPs = cell(size(trial_LFPs));
        for trial = 1:num_trials
            for chan = 1:size(trial_LFPs,2)
                if ~isempty(trial_LFPs{trial,chan})
                    envelope_trial_LFPs{trial,chan} = filtfilt(b2,a2,filtered_trial_LFPs{trial,chan});
                end
            end
        end
        
        % Calculate Mean and ~standad deviation
        chan_mean = zeros(1,length(LFPchannels));
        chan_std = zeros(1,length(LFPchannels));
        count = 0;
        for tr = 1:num_trials
            for chan = 1:size(trial_LFPs,2)
                if ~isempty(trial_LFPs{trial,chan})
                    chan_mean(chan) = chan_mean(chan)+mean(envelope_trial_LFPs{trial,chan});
                    chan_std(chan) = chan_std(chan)+std(envelope_trial_LFPs{trial,chan});
                    count = count+1;
                end
            end
        end
        chan_mean = chan_mean/(count/length(LFPchannels));
        chan_std = chan_std/(count/length(LFPchannels)); %approximate, not exact but should be ok as mean(std(i)) ~= std(:)

        
        proximity = cell(1,length(LFPchannels)); %how close to fixations start
        these_types = cell(1,length(LFPchannels)); %trial types
        HFOs = cell(1,length(LFPchannels));
        high_HFOs = cell(1,length(LFPchannels));
        high_HFOs2 = cell(1,length(LFPchannels));
        predicted = cell(1,length(LFPchannels));
        
        %---detect HFO events---%
        for trial = 1:num_trials
            for chan = 1:size(trial_LFPs,2)
                if isempty(trial_LFPs{trial,chan})
                    continue
                end
                [peaks,locs] = findpeaks(envelope_trial_LFPs{trial,chan},'MinPeakWidth',10);
                these_peaks = find(peaks > chan_mean(chan)+3*chan_std(chan));
                if isempty(these_peaks)
                    continue
                end
                peaks = peaks(these_peaks);
                locs = locs(these_peaks);
                
                for p = 1:length(peaks)
                    if locs(p) < 250 %started in previous trial, though could still be ITI
                        continue
                    elseif locs(p) > ons_offs(trial,9) %started after reward start
                        continue
                    elseif locs(p)+25 > length(trial_LFPs{trial,chan}) %too close to end of trial
                        continue
                    end
                    if any(abs(locs(p)-ons_offs(trial,:)) < 50) %cooccurs with screen flash
                        if trialtype(trial) == 1%sequence trials
                            seq_co_occur_count = seq_co_occur_count+1;
                        else %image trials
                            image_co_occur_count = image_co_occur_count+1;
                        end
                        continue
                    end
                    if all(envelope_trial_LFPs{trial,chan}(locs(p)-25:locs(p)+25)...
                            > (chan_mean(chan)+1*chan_std(chan)))
                        
                        less_than_thresh = find(envelope_trial_LFPs{trial,chan}...
                            < (chan_mean(chan)+1*chan_std(chan)));
                        onset = locs(p)-less_than_thresh;
                        onset(onset < 0) = [];
                        onset = min(onset);
                        locs(p) = locs(p)-onset;
                        if locs(p) < 251
                            continue
                        elseif locs(p)+500 > length(trial_LFPs{trial,chan}) %still too close to end of trial
                            continue
                        end
                                                                 
                        if trialtype(trial) == 1%sequence trials
                            seq_total_HFO_count =  seq_total_HFO_count+1;
                            seq_time_distribution(locs(p)) = seq_time_distribution(locs(p))+1;
                        else %image trials
                            if trialtype(trial) == 3 %novel image
                                image_total_HFO_count(1) =  image_total_HFO_count(1)+1;
                                image_time_distribution(1,locs(p)) = image_time_distribution(1,locs(p))+1;
                            elseif trialtype(trial) == 4 %repeat image
                                image_total_HFO_count(2) = image_total_HFO_count(2)+1;
                                image_time_distribution(2,locs(p)) = image_time_distribution(2,locs(p))+1;
                            end
                        end
                        
                        fixation_times = fixationstats{trial}.fixationtimes;
                        [~,closest_fixation] = min(abs(fixation_times(1,:)-locs(p)));
                        proximity{chan} = [proximity{chan} (fixation_times(1,closest_fixation)-locs(p))];
                        these_types{chan} = [these_types{chan} trialtype(trial)];
                        if trialtype(trial) == 1
                            if any(reaction_times(trial,:) < predict_rt)
                                predicted{chan} = [predicted{chan}  1];
                            else
                                predicted{chan} = [predicted{chan}  0];
                            end
                        else
                           predicted{chan} = [predicted{chan} NaN]; 
                        end
                        
                        temp = filtfilt(b,a,trial_LFPs{trial,chan}); %high pass filter
                        
                        HFOs{chan} = [HFOs{chan}; trial_LFPs{trial,chan}(locs(p)-250:locs(p)+500)]; %raw
                        high_HFOs{chan} = [high_HFOs{chan}; filtered_trial_LFPs{trial,chan}(locs(p)-250:locs(p)+500)]; %high pass
                        high_HFOs2{chan} = [high_HFOs2{chan}; temp(locs(p)-250:locs(p)+500)];
                    end
                end
            end
        end
        
        for chan = 1:length(HFOs)
            which_monkey = [which_monkey monkey];
            
            seq_trials = find(these_types{chan} == 1);
            image_trials = find(these_types{chan} > 1);
            
            %for sequence trials
            temp = HFOs{chan}(seq_trials,:);
            if ~isempty(temp)
                [~,~,wfq,meanpower,~,~,~] = waveletanalysis(temp);
                seq_HFO_trig_avg = [seq_HFO_trig_avg; nanmean(temp)];
                seq_high_HFO_trig_avg = [seq_high_HFO_trig_avg; nanmean(high_HFOs{chan}(seq_trials,:))];
                seq_high2_HFO_trig_avg = [seq_high2_HFO_trig_avg; nanmean(high_HFOs2{chan}(seq_trials,:))];
                was_predicted = [was_predicted predicted{chan}(seq_trials)];
                seq_all_trial_power = seq_all_trial_power+meanpower;
                sequence_proximity = [sequence_proximity  proximity{chan}(seq_trials)];
            end
            
            %for image trials
            temp = HFOs{chan}(image_trials,:);
            if ~isempty(temp)
                [~,~,wfq,meanpower,~,~,~] = waveletanalysis(temp);
                image_HFO_trig_avg = [image_HFO_trig_avg; nanmean(temp)];
                image_high_HFO_trig_avg = [image_high_HFO_trig_avg; nanmean(high_HFOs{chan}(image_trials,:))];
                image_high2_HFO_trig_avg = [image_high2_HFO_trig_avg; nanmean(high_HFOs2{chan}(image_trials,:))];
                image_all_trial_power = image_all_trial_power+meanpower;
                image_proximity = [image_proximity  proximity{chan}(image_trials)];
            end
        end
    end
end
%%
tm = -250:500;
figure
hold on
dofill(tm,seq_HFO_trig_avg,'black',1,1);
dofill(tm,image_HFO_trig_avg,'green',1,1);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
hold off
xlabel('Time From Ripple Onset (ms)')
ylabel('LFP Amplitude (uV)')
legend('Sequences','Images')
xlim([-250 500])
%%
figure
hold on
dofill(tm,seq_high_HFO_trig_avg,'black',1,1);
dofill(tm,image_high_HFO_trig_avg,'green',1,1);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
hold off
xlabel('Time From Ripple Onset (ms)')
ylabel('LFP Amplitude (uV)')
legend('Sequences','Images')
title('High Passed (> 100 Hz) and Rectified')
xlim([-250 500])
%%

figure
hold on
dofill(tm,seq_high2_HFO_trig_avg,'black',1,1);
dofill(tm,image_high2_HFO_trig_avg,'green',1,1);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
hold off
xlabel('Time From Ripple Onset (ms)')
ylabel('LFP Amplitude (uV)')
legend('Sequences','Images')
title('High Passed (> 100 Hz)')
xlim([-250 500])

%%
figure
subplot(2,2,1)
low = find(wfq < 30);
imagesc(tm,wfq(low),seq_all_trial_power(low,:)/size(seq_high2_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subplot(2,2,4)
mid = find(wfq >= 30 & wfq < 80);
imagesc(tm,wfq(mid),seq_all_trial_power(mid,:)/size(seq_high2_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subplot(2,2,2)
high = find(wfq > 70);
imagesc(tm,wfq(high),seq_all_trial_power(high,:)/size(seq_high2_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subtitle('Sequences')

figure
subplot(2,2,1)
low = find(wfq < 30);
imagesc(tm,wfq(low),image_all_trial_power(low,:)/size(image_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subplot(2,2,4)
mid = find(wfq >= 30 & wfq < 80);
imagesc(tm,wfq(mid),image_all_trial_power(mid,:)/size(image_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')


subplot(2,2,2)
high = find(wfq > 70);
imagesc(tm,wfq(high),image_all_trial_power(high,:)/size(image_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subtitle('Images')
%%
figure
subplot(2,2,1)
low = find(wfq < 30);
imagesc(tm,wfq(low),image_all_trial_power(low,:)/size(image_HFO_trig_avg,1)-seq_all_trial_power(low,:)/size(seq_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subplot(2,2,4)
mid = find(wfq >= 30 & wfq < 80);
imagesc(tm,wfq(mid),image_all_trial_power(mid,:)/size(image_HFO_trig_avg,1)-seq_all_trial_power(mid,:)/size(seq_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')


subplot(2,2,2)
high = find(wfq > 70);
imagesc(tm,wfq(high),image_all_trial_power(high,:)/size(image_HFO_trig_avg,1)-seq_all_trial_power(high,:)/size(seq_HFO_trig_avg,1))
xlabel('Time from HFO start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subtitle('Images-Sequences')

%%
figure
hold on
dofill((1:length(seq_time_distribution))-200,seq_time_distribution/seq_total_HFO_count,'black',1,30)
dofill((1:length(image_time_distribution(1,:)))-200,image_time_distribution(1,:)/(image_total_HFO_count(1)),'blue',1,30)
dofill((1:length(image_time_distribution(2,:)))-200,image_time_distribution(2,:)/(image_total_HFO_count(2)),'red',1,30)
legend('Sequences','Novel Images','Repeat Images')
xlabel('Time From Trial Start (ms)')
ylabel('Normalized Occurence Rate')
xlim([0 6500]) 

%%
figure
subplot(1,2,1)
hist(image_proximity,250)
xlim([-250 250])
title('Images')
xlabel('Lag from Closest Fixation Start (ms)')
ylabel('HFO Count')
title(['Median Lag: ' num2str(median(image_proximity)) ' ms'])
box off

subplot(1,2,2)
hist(sequence_proximity,250)
title('Sequence')
xlim([-250 250])
xlabel('Lag from Closest Fixation Start (ms)')
ylabel('HFO Count')
title(['Median Lag: ' num2str(median(sequence_proximity)) ' ms'])
box off
