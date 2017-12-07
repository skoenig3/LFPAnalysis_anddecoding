% clar

task = 'ListSQ';
twin = 1000;% how much time to take before and after saccade.
image_on_twin = 500;
fixwin = 5;%size of fixation window on each crosshair
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
imgon_code = 23;
imgoff_code = 24;
Fs = 1000;
min_blks = 2;


min_fix_dur = 100;
min_sac_amp = 2;
max_amplitude = 16;
min_amplitude = 2;
bin_amplitude = 4;
% amps = [0:bin_amplitude:18]+2;
% amps(end) = 24;
amps = [2 4 6 10 16 24];

directs = [0 90 180 -90];
%create 30 Hz low pass filter
fs = 1000;
[blow,alow] = butter(6,30/(fs/2),'low');
[bhigh,ahigh] = butter(6,30/(fs/2),'high');


LFP_count = 0;

all_amplitude_LFP = zeros(length(amps)-1,2*twin+1);
all_direction_LFP = zeros(4,2*twin+1);
all_LFPs = zeros(1,2*twin+1);
up_vs_down = zeros(1,2*twin+1);
left_vs_right = zeros(1,2*twin+1);


all_power = zeros(27,2*twin+1);
left_vs_right_power = zeros(27,2*twin+1);
up_vs_down_power = zeros(27,2*twin+1);
small_vs_large_saccade_power =zeros(27,2*twin+1);


for monkey = 2%1%1:2
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
        for chan = 1:length(LFPchannels)
            fixation_aligned_LFP{chan} = NaN(num_trials*15,2*twin+1);
        end
        
        fixind = ones(1,length(LFPchannels));
        for trial = 1:num_trials
            if sum(cfg.trl(trial).allval == 3) >= 5; %in which sequence trials were rewarded
                last_code = 3;
                trial_type = 1;
            elseif any(cfg.trl(trial).allval == 23) && itmlist(cfg.trl(trial).cnd-1000) > 19 %successful image trial
                last_code = imgoff_code;
                trial_type = 2;
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
            
            
            %find fixations and saccades that did not occur during the image period;
            %should also take care of the 1st fixation on the crosshair
            
            %fixation started before image turned on
            invalid= find(fixationtimes(1,:) < imgon);
            fixationtimes(:,invalid) = [];
            
            %fixation ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(fixationtimes(2,:) > imgoff);
            fixationtimes(:,invalid) = [];
            
            %saccade started before image turned on
            invalid= find(saccadetimes(1,:) < imgon);
            saccadetimes(:,invalid) = [];
            
            %Saccaded ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(saccadetimes(2,:) > imgoff);
            saccadetimes(:,invalid) = [];
            
            saccade_amplitude = NaN(1,size(fixationtimes,2));
            fixation_duration = NaN(1,size(fixationtimes,2));
            saccade_start_time = NaN(1,size(fixationtimes,2));
            saccade_direction = NaN(1,size(fixationtimes,2));
            for f = 1:size(fixationtimes,2);
                prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
                if isempty(prior_sac) %trial ended or eye outside of image
                    continue %try next one
                end
                sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%this fixation duration
                if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough
                    saccade_amplitude(f) = sacamp;
                    saccade_start_time(f) = saccadetimes(1,prior_sac);
                    fixation_duration(f) = fixationtimes(2,f)-fixationtimes(1,f)+1;
                    saccade_direction(f) = atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                        xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%saccade_direction
                end
            end
            
            for chan = 1:length(LFPchannels);
                trial_LFP = data(LFPchannels(chan)).values{trial};
                %trial_LFP = filtfilt(blow,alow,trial_LFP); %low pass filter
                trial_LFP = filtfilt(bhigh,ahigh,trial_LFP); %high pass filter
                for f = 1:length(saccade_start_time)
                    if ~isnan(saccade_start_time(f))
                        fixation_aligned_LFP{chan}(fixind(chan),:) = trial_LFP(saccade_start_time(f)-twin:saccade_start_time(f)+twin);
                        all_saccade_amplitudes(chan,fixind(chan)) = saccade_amplitude(f);
                        all_saccade_directions(chan,fixind(chan)) = saccade_direction(f);
                        all_fixation_durations(chan,fixind(chan)) =  fixation_duration(f);
                        fixind(chan) = fixind(chan)+1;
                    end
                end
            end
        end
        
        amplitude_LFP = cell(1,length(LFPchannels));
        direction_LFP = cell(1,length(LFPchannels));
        maximum_value = NaN(1,length(LFPchannels));
        channel_power = cell(1,length(LFPchannels));
        for chan = 1:length(LFPchannels)
            fixation_aligned_LFP{chan}(fixind:end,:) = [];
            amplitude_LFP{chan} = NaN(length(amps)-1,2*twin+1);
            direction_LFP{chan} = NaN(4,2*twin+1);
            maximum_value(chan) = max(abs(mean(fixation_aligned_LFP{chan})));
            all_LFPs = all_LFPs+mean(fixation_aligned_LFP{chan})/maximum_value(chan);
            temp = fixation_aligned_LFP{chan}/maximum_value(chan);
             [~,~,wfq,meanpower,~,~,~] = waveletanalysis(temp(1:4:end,:)); %downsample
             all_power = all_power+meanpower;
        end
        all_saccade_amplitudes(:,fixind(1):end,:) = [];
        all_saccade_directions(:,fixind(1):end,:) = [];
        all_saccade_amplitudes = all_saccade_amplitudes/24; %convert to dva
        all_fixation_durations(:,fixind(1):end,:) = [];
        %%

        for chan = 1:length(LFPchannels)
            for bin = 1:length(amps)-1
                these_amplitudes = (all_saccade_amplitudes(1,:) < amps(bin+1) & all_saccade_amplitudes(1,:) >= amps(bin));
                amplitude_LFP{chan}(bin,:) = nanmean(fixation_aligned_LFP{chan}(these_amplitudes,:));
            end
            all_amplitude_LFP = all_amplitude_LFP+amplitude_LFP{chan}/maximum_value(chan);
            prct25 = prctile(all_saccade_amplitudes(1,:),25);
            prct75 = prctile(all_saccade_amplitudes(1,:),75);
            small = find(all_saccade_amplitudes(1,:) <= prct25);
            large = find(all_saccade_amplitudes(1,:) >= prct75);
            
            temp = fixation_aligned_LFP{chan}(small,:);
            temp = temp/maximum_value(chan);
            [~,~,wfq,smallmeanpower,~,~,~] = waveletanalysis(temp); 
            
            temp = fixation_aligned_LFP{chan}(large,:);
            temp = temp/maximum_value(chan);
            [~,~,wfq,largemeanpower,~,~,~] = waveletanalysis(temp); 
            small_vs_large_saccade_power = small_vs_large_saccade_power+(largemeanpower-smallmeanpower);
            
            dr_power = cell(1,4);
            for dr =1:4
                dirs = all_saccade_directions(1,:);
                if dr == 3 %180
                    these_directions = ((dirs > 150) | (dirs < -150));
                else
                    these_directions = ((dirs >= directs(dr)-30) & (dirs <=  directs(dr)+30));
                end
                direction_LFP{chan}(dr,:) = nanmean(fixation_aligned_LFP{chan}(these_directions,:));
                
                 temp = fixation_aligned_LFP{chan}(these_directions,:);
                 temp = temp/maximum_value(chan);
                 [~,~,wfq,dr_power{dr},~,~,~] = waveletanalysis(temp); 

            end
            all_direction_LFP = all_direction_LFP+direction_LFP{chan}/maximum_value(chan);
            
            left_vs_right = left_vs_right+...
                (direction_LFP{chan}(1,:)-direction_LFP{chan}(3,:))/maximum_value(chan);
            up_vs_down = up_vs_down +...
                (direction_LFP{chan}(2,:)-direction_LFP{chan}(4,:))/maximum_value(chan);
            
            left_vs_right_power = left_vs_right_power + (dr_power{1}-dr_power{3});
            up_vs_down_power = up_vs_down_power +(dr_power{2}-dr_power{4});
        end
        LFP_count = LFP_count+length(LFPchannels);
    end
end

%%
tm = -twin:twin;

figure
subplot(2,2,1)
plot(tm,all_LFPs/LFP_count);
xlabel('Time From Saccade Start (ms)')
ylabel('Normalized LFP')
xlim([-500 500])

subplot(2,2,2)
plot(tm,all_amplitude_LFP'/LFP_count)
legend('2-4','4-6','6-10','10-16','16-24')
xlabel('Time From Saccade Start (ms)')
ylabel('Normalized LFP')
xlim([-500 500])

subplot(2,2,3)
plot(tm,all_direction_LFP'/LFP_count)
legend('Right','Up','Left','Down')
xlabel('Time From Saccade Start (ms)')
ylabel('Normalized LFP')
xlim([-500 500])

subplot(2,2,4)
hold on
plot(tm,left_vs_right/LFP_count)
plot(tm,up_vs_down/LFP_count)
plot([-twin twin],[0 0],'k--')
hold off
xlabel('Time From Saccade Start (ms)')
ylabel('Normalized LFP')
legend('Left-Right','Up-Down')
xlim([-500 500])

if monkey == 1
    subtitle(['Vivian, n = ' num2str(LFP_count)])
else
    subtitle(['Tobii, n = ' num2str(LFP_count)])
end
%%
figure
subplot(2,2,1)
imagesc(tm,wfq,all_power/LFP_count);
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
colormap('jet')
title('All Saccades')
xlim([-500 500])

subplot(2,2,2)
imagesc(tm,wfq,small_vs_large_saccade_power/LFP_count);
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
colormap('jet')
title('25% Largest-25% Smallest Saccades')
xlim([-500 500])

subplot(2,2,3)
imagesc(tm,wfq,up_vs_down_power/LFP_count);
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
colormap('jet')
title('Up vs Down')
xlim([-500 500])

subplot(2,2,4)
imagesc(tm,wfq,left_vs_right_power/LFP_count);
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
colormap('jet')
title('Left vs Right')
xlim([-500 500])

