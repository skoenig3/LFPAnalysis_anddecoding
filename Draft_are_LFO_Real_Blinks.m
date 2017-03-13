clar

task = 'ListSQ';
twin = 750;% how much time to take before and after saccade.
twin2 = 250;
image_on_twin = 500;
fixwin = 5;%size of fixation window on each crosshair
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
imgon_code = 23;
imgoff_code = 24;
Fs = 1000;
min_blks = 2;
dotsize=6; %marker area; point width squared
min_blink_dur = 25;

%create 30 Hz low pass filter
fs = 1000;
[blow,alow] = butter(6,12/(fs/2),'low');
[bhigh,ahigh] = butter(6,30/(fs/2),'high');


LFP_count = 0;
avg_blink_LFP = [];
avg_blink_pupil = [];

%
for monkey = 2%1:2
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
    
    for sess =1:length(session_data)
        
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
        
        pupilchannel= find_desired_channels(cfg,'pupil');
        if isempty(pupilchannel)
            continue
        end
        
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
        
        blink_LFP = cell(1,length(LFPchannels));
        for chan = 1:length(LFPchannels)
            blink_LFP{chan} = NaN(2000,twin+twin2+1);
        end
        blink_pupil = NaN(2000,twin+twin2+1);
        index = 1;
        for trial = 1:num_trials
            if sum(cfg.trl(trial).allval == 3) >= 5; %in which sequence trials were rewarded
                last_code = 3;
                trial_type = 1;
                %continue
            elseif any(cfg.trl(trial).allval == 23) && itmlist(cfg.trl(trial).cnd-1000) > 19 %successful image trial
                last_code = imgoff_code;
                trial_type = 2;
                %continue
            else
                continue
            end
            
            trial_start = cfg.trl(trial).alltim(cfg.trl(trial).allval == trial_start_code); %time at start of trial
            imgon = cfg.trl(trial).alltim(cfg.trl(trial).allval == imgon_code)-trial_start; %time of image on relative to start of trial
            imgoff = cfg.trl(trial).alltim(cfg.trl(trial).allval == last_code)-trial_start; %time of image off relative to start of trial
            imgoff = imgoff(1); %if sequence then last code is reward code so many of them
            
            % if monkey isn't paying attention data isn't probably
            % worth much plus have to cut off somewhere
            % will also cap sequence trial durations, which should be
            % shorter anyway
            if imgoff-imgon > 1.5*imgdur-1
                imgoff = imgon+1.5*imgdur-1;
            end
            imgoff = imgoff-twin;
            
            
            pupil = data(pupilchannel).values{trial};
            thresh = mean(pupil)-3*std(pupil);
            blinks = find(pupil < thresh);
            blinks = findgaps(blinks);
            if size(blinks,2) == 1;
                blinks = blinks';
            end
            
            blink_ind = [];%upsampled ind
            if ~isempty(blinks)
                for b = 1:size(blinks,1)
                    ind = blinks(b,:);
                    ind(ind == 0) = [];
                    if ind(1) <= 1000%twin2 %too early in the trial
                        continue
                    end
                    if ind(end) > imgoff %too late
                        continue
                    end
                    if length(ind) > min_blink_dur
                        continue
                    end
                    blink_ind = [blink_ind ind(1)];%center on up sampled bins
                end
            end
            
            for b = 1:length(blink_ind)
                blink_pupil(index+b-1,:) = pupil(blink_ind(b)-twin2:blink_ind(b)+twin);
            end
            
            for b = 1:length(blink_ind)
                for chan = 1:length(LFPchannels);
                    trial_LFP = data(LFPchannels(chan)).values{trial};
                    %trial_LFP = filtfilt(blow,alow,trial_LFP); %low pass filter
                    %trial_LFP = filtfilt(bhigh,ahigh,trial_LFP); %high pass filter
                    blink_LFP{chan}(index,:) = trial_LFP(blink_ind(b)-twin2:blink_ind(b)+twin);
                end
                index = index+1;
            end
        end
        for chan = 1:length(LFPchannels);
            avg_blink_LFP = [avg_blink_LFP; nanmean(blink_LFP{chan})];
            avg_blink_pupil = [avg_blink_pupil; nanmean(blink_pupil)];
        end
    end
end
%%
tm = -twin2:twin;
figure
subplot(1,2,1)
plot(tm,nanmean(avg_blink_pupil))
xlabel('Time from Blink Start (ms)')
ylabel('Pupil Diameter (a.u.)')
xlim([-twin2 twin])

subplot(1,2,2)
plot(tm,nanmean(avg_blink_LFP))
xlabel('Time from Blink Start (ms)')
ylabel('LFP (uV)')
xlim([-twin2 twin])



if monkey == 1
    subtitle('Vivian')
else
    subtitle('Tobii')
end