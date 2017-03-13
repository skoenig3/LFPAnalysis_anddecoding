clar

task = 'ListSQ';
twin = 500;% how much time to take before and after saccade.
image_on_twin = 1000;
fixwin = 5;%size of fixation window on each crosshair
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
imgon_code = 23;
imgoff_code = 24;
Fs = 1000;
min_blks = 2;
dotsize=6; %marker area; point width squared
%sequence PCA, 5 dva min amp, 400 ms min fix dur!
min_fix_dur = 250;
min_sac_amp = 2*24; %in pixels

imageX = 800;
imageY = 600;

directs = [0 90 180 -90];
fs = 1000;
[blow,alow] = butter(6,20/(fs/2),'low');%create 20 Hz low pass filter
[bhigh,ahigh] = butter(6,30/(fs/2),'high');%create 30 Hz high pass filter



all_directions = [];
all_amplitudes = [];
all_positions = [];
all_fixdurs = [];
all_LFPs = [];

for monkey = 1%1:2
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
        
        fixation_aligned_LFP =cell(1,length(LFPchannels));
        all_saccade_amplitudes = NaN(length(LFPchannels),num_trials*15);
        all_saccade_directions = NaN(length(LFPchannels),num_trials*15);
        all_fixation_durations = NaN(length(LFPchannels),num_trials*15);
        all_fixation_locations = cell(1,length(LFPchannels));
        for chan = 1:length(LFPchannels)
            fixation_aligned_LFP{chan} = NaN(num_trials*15,2*twin+1);
            all_fixation_locations{chan} = NaN(num_trials*15,2);
        end
        
        fixind = ones(1,length(LFPchannels));
        for trial = 1:num_trials
            if sum(cfg.trl(trial).allval == 3) >= 5; %in which sequence trials were rewarded
                last_code = 3;
                trial_type = 1;
                continue
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
            
            fixationtimes = fixationstats{trial}.fixationtimes; %start and end times of fixations
            saccadetimes = fixationstats{trial}.saccadetimes; % start and end times of saccades
            xy = fixationstats{trial}.XY; %x and y eye data
            fixations = round(fixationstats{trial}.fixations);
            fixations(1,fixations(1,:) > imageX) = imageX;
            fixations(2,fixations(2,:) > imageY) = imageY;
            fixations(fixations < 1) = 1;
            
            
            %find fixations and saccades that did not occur during the image period;
            %should also take care of the 1st fixation on the crosshair
            
            %fixation started before image turned on
            invalid= find(fixationtimes(1,:) < imgon);
            fixationtimes(:,invalid) = [];
            fixations(:,invalid) = [];
            
            %fixation ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(fixationtimes(2,:) > imgoff);
            fixationtimes(:,invalid) = [];
            fixations(:,invalid) = [];
            
            %saccade started before image turned on
            invalid= find(saccadetimes(1,:) < imgon);
            saccadetimes(:,invalid) = [];
            
            %Saccaded ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(saccadetimes(2,:) > imgoff);
            saccadetimes(:,invalid) = [];
            
            saccade_amplitude = NaN(1,size(fixationtimes,2));
            fixation_duration = NaN(1,size(fixationtimes,2));
            fixation_start_time = NaN(1,size(fixationtimes,2));
            saccade_direction = NaN(1,size(fixationtimes,2));
            fixation_location = NaN(2,size(fixationtimes,2));
            for f = 1:size(fixationtimes,2);
                prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
                if isempty(prior_sac) %trial ended or eye outside of image
                    continue %try next one
                end
                sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%this fixation duration
                if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough
                    saccade_amplitude(f) = sacamp;
                    fixation_start_time(f) =fixationtimes(1,f);%fixationtimes(1,f);
                    fixation_duration(f) = fixationtimes(2,f)-fixationtimes(1,f)+1;
                    saccade_direction(f) = atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                        xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%saccade_direction
                    fixation_location(:,f) = fixations(:,f);
                end
            end
            
            for chan = 1:length(LFPchannels);
                trial_LFP = data(LFPchannels(chan)).values{trial};
                trial_LFP = filtfilt(blow,alow,trial_LFP); %low pass filter
                %trial_LFP = filtfilt(bhigh,ahigh,trial_LFP); %high pass filter
                for f = 1:length(fixation_start_time)
                    if ~isnan(fixation_start_time(f))
                        this_fixation_LFP = trial_LFP(fixation_start_time(f)-twin:fixation_start_time(f)+twin);
                        fixation_aligned_LFP{chan}(fixind(chan),:) = this_fixation_LFP;
                  
                        all_saccade_amplitudes(chan,fixind(chan)) = saccade_amplitude(f);
                        all_saccade_directions(chan,fixind(chan)) = saccade_direction(f);
                        all_fixation_durations(chan,fixind(chan)) =  fixation_duration(f);
                        all_fixation_locations{chan}(fixind(chan),:) =  fixation_location(:,f)';
                        
                        fixind(chan) = fixind(chan)+1;
                        
                    end
                end
            end
        end
        
        %---remove excess nans---%
        fixation_aligned_LFP = laundry(fixation_aligned_LFP);
        all_saccade_amplitudes = laundry(all_saccade_amplitudes);
        all_saccade_directions = laundry(all_saccade_directions);
        all_fixation_durations = laundry(all_fixation_durations);
        all_fixation_locations = laundry(all_fixation_locations);
        
        for chan = 1:length(LFPchannels);
            all_directions = [all_directions; all_saccade_directions(chan,:)'];
            all_amplitudes = [all_amplitudes; all_saccade_amplitudes(chan,:)'];
            all_positions = [all_positions; all_fixation_locations{chan}];
            all_fixdurs = [all_fixdurs; all_fixation_durations(chan,:)'];
            all_LFPs = [all_LFPs; fixation_aligned_LFP{chan}];
            
        end
       
    end
end

%% Plot Some Statistics
tm = -twin:twin;
figure
plot(tm,nanmean(all_LFPs));
yl = ylim;
hold on
plot([0 0],[yl(1) yl(2)],'k--')
plot([min_fix_dur min_fix_dur],[yl(1) yl(2)],'k--')
hold off
xlabel('Time From Fixation Start (ms)')
ylabel('LFP (uV)')
title('Avg Across Session/Electrode')
%%
afd = all_fixdurs;
afd(afd > 1000) = [];
figure
hist(afd,50);
xlabel('Fixation Duration (ms)')
ylabel('Count')

figure
subplot(1,2,1)
hist((all_positions(:,1)-400)/24,32)
xlabel('Horizontal Postion (dva)')
ylabel('Fixation Count')

subplot(1,2,2)
hist((all_positions(:,2)-300)/24,32)
xlabel('Vertical Postion (dva)')
ylabel('Fixation Count')

eccentricity = all_positions;
eccentricity(:,1) = (eccentricity(:,1)-400)/24;
eccentricity(:,2) = (eccentricity(:,2)-300)/24;
eccentricity = sqrt(sum(eccentricity.^2,2));

figure
hist(eccentricity,20)
xlabel('Eccentricity from Center of Screen (dva)')
ylabel('Fixation Count')


figure
hist(all_amplitudes/24,25)
xlabel('Preceding Saccade Amplitude (dva)')
ylabel('Fixation Count')

figure
hist(all_directions,360)
xlabel('Preceding Saccade Angle')
ylabel('Fixation Count')

%% Compute the PCAs
[U,S,V] = pca(all_LFPs(:,twin-50:twin+min_fix_dur),6);
svar = S/sum(S(:));
%% Visualize the PCAs


figure
for p = 1:size(U,2)
   low_ind = find(U(:,p) < prctile(U(:,p),25)); 
   high_ind = find(U(:,p) >  prctile(U(:,p),75)); 
   
   subplot(2,3,p)
   hold on
   plot(tm,nanmean(all_LFPs(low_ind,:)));
   plot(tm,nanmean(all_LFPs(high_ind,:)));
   hold off
   xlabel('Time From Fixation Start (ms)')
   ylabel('LFP (uV)')
   title(['PCA#' num2str(p)])
end
