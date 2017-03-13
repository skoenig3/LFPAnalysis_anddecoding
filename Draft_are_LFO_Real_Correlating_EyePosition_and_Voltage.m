clar

task = 'ListSQ';
twin = 1000;% how much time to take before and after sapos_CCade.
image_on_twin = 500;
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

min_cluster_area = 250; %ms*Hz for multiple comparisons corrections with maris method
numshuffs = 1000;


imageX = 800;
imageY = 600;

directs = [0 90 180 -90];
%create 30 Hz low pass filter
fs = 1000;
[blow,alow] = butter(6,16/(fs/2),'low');
[bhigh,ahigh] = butter(6,30/(fs/2),'high');


LFP_count = 0;
all_LFPs = zeros(1,2*twin+1);


all_left_LFPs = zeros(1,2*twin+1);
all_right_LFPs = zeros(1,2*twin+1);
all_up_LFPs = zeros(1,2*twin+1);
all_down_LFPs = zeros(1,2*twin+1);
all_center_lr_LFPs = zeros(1,2*twin+1);

all_left_count = 0;
all_right_count = 0;
all_up_count = 0;
all_down_count = 0;
all_center_lr_count = 0;


left_center_lr_sapos_CCades = zeros(1,2*twin+1);
left_left_sapos_CCades = zeros(1,2*twin+1);
right_center_lr_sapos_CCades = zeros(1,2*twin+1);
right_right_sapos_CCades = zeros(1,2*twin+1);

allLc = [];
allRc = [];
allLL = [];
allRR = [];

left_center_lr_count = 0;
left_left_count = 0;
right_center_lr_count = 0;
right_right_count = 0;


all_fixation_voltages = zeros(600,800);
all_fixation_voltages_count = zeros(600,800);

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
        %%%---Get supos_CCessful trials Information by Task---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        oveRRap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
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
        all_sapos_CCade_amplitudes = NaN(length(LFPchannels),num_trials*15);
        all_sapos_CCade_directions = NaN(length(LFPchannels),num_trials*15);
        all_fixation_durations = NaN(length(LFPchannels),num_trials*15);
        for chan = 1:length(LFPchannels)
            fixation_aligned_LFP{chan} = NaN(num_trials*15,2*twin+1);
        end
        
        LFPs_by_Location = cell(length(LFPchannels),4);
        fix_ind2 = ones(length(LFPchannels),4);
        for chan = 1:length(LFPchannels);
            for l = 1:4
                LFPs_by_Location{chan,l} = NaN(1000,2*twin+1);
            end
        end
        
        
        fixind = ones(1,length(LFPchannels));
        for trial = 1:num_trials
            if sum(cfg.trl(trial).allval == 3) >= 5; %in which sequence trials were rewarded
                last_code = 3;
                trial_type = 1;
                continue
            elseif any(cfg.trl(trial).allval == 23) && itmlist(cfg.trl(trial).cnd-1000) > 19 %supos_CCessful image trial
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
            saccadetimes = fixationstats{trial}.saccadetimes; % start and end times of sapos_CCades
            xy = fixationstats{trial}.XY; %x and y eye data
            fixations = round(fixationstats{trial}.fixations);
            fixations(1,fixations(1,:) > imageX) = imageX;
            fixations(2,fixations(2,:) > imageY) = imageY;
            fixations(fixations < 1) = 1;
            
            
            %find fixations and sapos_CCades that did not opos_CCur during the image period;
            %should also take care of the 1st fixation on the crosshair
            
            %fixation started before image turned on
            invalid= find(fixationtimes(1,:) < imgon);
            fixationtimes(:,invalid) = [];
            fixations(:,invalid) = [];
            
            %fixation ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(fixationtimes(2,:) > imgoff);
            fixationtimes(:,invalid) = [];
            fixations(:,invalid) = [];
            
            %sapos_CCade started before image turned on
            invalid= find(saccadetimes(1,:) < imgon);
            saccadetimes(:,invalid) = [];
            
            %Sapos_CCaded ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(saccadetimes(2,:) > imgoff);
            saccadetimes(:,invalid) = [];
            
            sapos_CCade_amplitude = NaN(1,size(fixationtimes,2));
            fixation_duration = NaN(1,size(fixationtimes,2));
            fixation_start_time = NaN(1,size(fixationtimes,2));
            sapos_CCade_direction = NaN(1,size(fixationtimes,2));
            fixation_location = NaN(2,size(fixationtimes,2));
            
            for f = 1:size(fixationtimes,2);
                prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
                if isempty(prior_sac) %trial ended or eye outside of image
                    continue %try next one
                end
                sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %sapos_CCade amplitude
                %
                %check where sapos_CCade started
                %                 x_start = xy(1,saccadetimes(1,prior_sac));
                %                 y_start = xy(2,saccadetimes(1,prior_sac));
                %                 if x_start < 266 || x_start > 533 || y_start < 200 || y_start > 400
                %                     continue
                %                 end
                %
                fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%this fixation duration
                if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough
                    sapos_CCade_amplitude(f) = sacamp;
                    fixation_start_time(f) =fixationtimes(1,f);%fixationtimes(1,f);
                    fixation_duration(f) = fixationtimes(2,f)-fixationtimes(1,f)+1;
                    sapos_CCade_direction(f) = atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                        xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%sapos_CCade_direction
                    fixation_location(:,f) = fixations(:,f);
                end
            end
            
            for chan = 1:length(LFPchannels);
                trial_LFP = data(LFPchannels(chan)).values{trial};
                %trial_LFP = filtfilt(blow,alow,trial_LFP); %low pass filter
                %trial_LFP = filtfilt(bhigh,ahigh,trial_LFP); %high pass filter
                for f = 1:length(fixation_start_time)
                    if ~isnan(fixation_start_time(f))
                        this_fixation_LFP = trial_LFP(fixation_start_time(f)-twin:fixation_start_time(f)+twin);
                        fixation_aligned_LFP{chan}(fixind(chan),:) = this_fixation_LFP;
                        
                        all_fixation_voltages(fixation_location(2,f),fixation_location(1,f)) = ...
                            all_fixation_voltages(fixation_location(2,f),fixation_location(1,f))+...
                            mean(trial_LFP(fixation_start_time(f)+100:fixation_start_time(f)+200));
                        all_fixation_voltages_count(fixation_location(2,f),fixation_location(1,f)) = ...
                            all_fixation_voltages_count(fixation_location(2,f),fixation_location(1,f)) +1;
                        
                        all_sapos_CCade_amplitudes(chan,fixind(chan)) = sapos_CCade_amplitude(f);
                        all_sapos_CCade_directions(chan,fixind(chan)) = sapos_CCade_direction(f);
                        all_fixation_durations(chan,fixind(chan)) =  fixation_duration(f);
                        fixind(chan) = fixind(chan)+1;
                        
                        %---Fixation Aligned LFPs by Screen Position and Sapos_CCade Direction---%
                        if fixations(1,f) < 267
                            all_left_LFPs = all_left_LFPs+this_fixation_LFP;
                            all_left_count =all_left_count+1;
                            
                            if sapos_CCade_direction(f) > 150 || sapos_CCade_direction(f) < -150
                                left_left_sapos_CCades = left_left_sapos_CCades+this_fixation_LFP;
                                left_left_count = left_left_count+1;
                                
                                LFPs_by_Location{chan,1}(fix_ind2(chan,1),:) = this_fixation_LFP;
                                fix_ind2(chan,1) = fix_ind2(chan,1)+1;
                            end
                        elseif fixations(1,f) > 534
                            all_right_LFPs = all_right_LFPs+this_fixation_LFP;
                            all_right_count =all_right_count+1;
                            
                            if sapos_CCade_direction(f) < 30 && sapos_CCade_direction(f) > -30
                                right_right_sapos_CCades = right_right_sapos_CCades+this_fixation_LFP;
                                right_right_count = right_right_count+1;
                                LFPs_by_Location{chan,2}(fix_ind2(chan,2),:) = this_fixation_LFP;
                                fix_ind2(chan,2) = fix_ind2(chan,2)+1;
                            end
                        end
                        
                        if fixations(2,f) < 200
                            all_up_LFPs = all_up_LFPs+this_fixation_LFP;
                            all_up_count = all_up_count+1;
                        elseif fixations(2,f) > 400;
                            all_down_LFPs=all_down_LFPs+this_fixation_LFP;
                            all_down_count = all_down_count+1;
                        end
                        
                        if fixations(1,f) > 266 && fixations(1,f) < 534 && ...
                                fixations(2,f) > 200 && fixations(2,f) < 400
                            all_center_lr_LFPs = all_center_lr_LFPs+this_fixation_LFP;
                            all_center_lr_count = all_center_lr_count+1;
                            
                            if sapos_CCade_direction(f) > 150 || sapos_CCade_direction(f) < -150
                                left_center_lr_sapos_CCades = left_center_lr_sapos_CCades+this_fixation_LFP;
                                left_center_lr_count = left_center_lr_count+1;
                                LFPs_by_Location{chan,3}(fix_ind2(chan,3),:) = this_fixation_LFP;
                                fix_ind2(chan,3) = fix_ind2(chan,3)+1;
                            elseif sapos_CCade_direction(f) < 30 && sapos_CCade_direction(f) > -30
                                right_center_lr_sapos_CCades = right_center_lr_sapos_CCades+this_fixation_LFP;
                                right_center_lr_count = right_center_lr_count+1;
                                LFPs_by_Location{chan,4}(fix_ind2(chan,4),:) = this_fixation_LFP;
                                fix_ind2(chan,4) = fix_ind2(chan,4)+1;
                            end
                        end
                    end
                end
            end
        end
        LFPs_by_Location = laundry(LFPs_by_Location);
        
        
        amplitude_LFP = cell(1,length(LFPchannels));
        direction_LFP = cell(1,length(LFPchannels));
        maximum_value = NaN(1,length(LFPchannels));
        for chan = 1:length(LFPchannels)
            fixation_aligned_LFP{chan}(fixind:end,:) = [];
            all_LFPs = all_LFPs+nanmean(fixation_aligned_LFP{chan});
            LFP_count = LFP_count+1;
            
            
            allLc = [allLc; LFPs_by_Location{chan,3}];
            allRc = [allRc; LFPs_by_Location{chan,4}];
            allLL = [allLL; LFPs_by_Location{chan,1}];
            allRR = [allRR; LFPs_by_Location{chan,2}];
        end
    end
end

%%
filter_width = 1;
filter_size = filter_width*6+1;
H = fspecial('gaussian',filter_size,filter_width);
bin_count = bin2(all_fixation_voltages_count,24,24,'upper','sum');
smoothed_count = imfilter(bin_count,H);
bin_all_votlages = bin2(all_fixation_voltages,24,24,'upper','mean');
bin_all_voltages = bin_all_votlages./bin_count;
smooth_all_voltages = imfilter(bin_all_votlages,H);


figure
subplot(2,3,1)
imagesc(smoothed_count)
title('Fixation Density')
axis off

subplot(2,3,4)
imagesc(smooth_all_voltages)
colormap('jet')
colorbar
title('Average Post Fixation Voltage')

subplot(2,3,2)
plot([-16:16],mean(smooth_all_voltages))
xlabel('Horizontal Eye Position (dva)')
ylabel('Avg. Voltage (uV)')
xlim([-16 16])


subplot(2,3,5)
plot([-12:12],mean(smooth_all_voltages'))
xlabel('Vertical Eye Position (dva)')
ylabel('Avg. Voltage (uV)')
xlim([-12 12])

tm = -twin:twin;
subplot(2,3,3)
plot(tm,all_LFPs/LFP_count)
xlabel('Time from Fixation Start (ms)')
ylabel('LFP (uV)')
xlim([-100 500])

subplot(2,3,6)
hold on
plot(tm,all_left_LFPs/all_left_count)
plot(tm,all_right_LFPs/all_right_count);
plot(tm,all_up_LFPs/all_up_count);
plot(tm,all_down_LFPs/all_down_count);
plot(tm,all_center_lr_LFPs/all_center_lr_count,'k')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('LFP (uV)')
xlim([-100 500])
legend('Left','Right','Up','Down','Center')

if monkey == 1
    subtitle(['Vivian'])
else
    subtitle(['Tobii'])
end

%%
figure

hold on
plot(tm,left_left_sapos_CCades/left_left_count)
plot(tm,left_center_lr_sapos_CCades/left_center_lr_count);
plot(tm,right_right_sapos_CCades/right_right_count);
plot(tm,right_center_lr_sapos_CCades/right_center_lr_count);
plot(tm,all_center_lr_LFPs/all_center_lr_count,'k')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('LFP (uV)')
xlim([-100 500])
legend('Left-Left','Left-Center','Right-Right','Right-Center','All Center')


if monkey == 1
    subtitle(['Vivian'])
else
    subtitle(['Tobii'])
end

%% Calculate Observed Phase and Power Differences between Different Postion x Direction Combinations
max_samples = min([size(allLc,1),size(allRc,1),size(allLL,1),size(allRR,1)]);

Lc_ind = randperm(size(allLc,1));
Lc_ind = Lc_ind(1:max_samples);
this_allLc = allLc(Lc_ind,:);
[LC_trialpower,LC_trialphase,wfq,LC_meanpower,LC_meanphase,LC_powervar,LC_phasevar] = waveletanalysis(this_allLc);

RC_ind = randperm(size(allRc,1));
RC_ind = RC_ind(1:max_samples);
this_allRc = allRc(RC_ind,:);
[RC_trialpower,RC_trialphase,wfq,RC_meanpower,RC_meanphase,RC_powervar,RC_phasevar] = waveletanalysis(this_allRc);


LL_ind = randperm(size(allLL,1));
LL_ind = LL_ind(1:max_samples);
this_allLL = allLL(LL_ind,:);
[LL_trialpower,LL_trialphase,wfq,LL_meanpower,LL_meanphase,LL_powervar,LL_phasevar] = waveletanalysis(this_allLL);

RR_ind = randperm(size(allRR,1));
RR_ind = RR_ind(1:max_samples);
this_allRR = allRR(RR_ind,:);
[RR_trialpower,RR_trialphase,wfq,RR_meanpower,RR_meanphase,RR_powervar,RR_phasevar] = waveletanalysis(this_allRR);

[~,~,wfq,avg_meanpower,avg_meanphase,avg_powervar,avg_phasevar] = waveletanalysis(...
    [this_allLc(1:4:end,:); this_allLL(1:4:end,:); this_allRc(1:4:end,:); this_allRR(1:4:end,:)]);

%% Plot Observed Phase and Power Differences between Different Postion x Direction Combinations
figure
subplot(2,2,1)
imagesc(tm,wfq,LC_meanpower)
axis xy
title('Left-Center Fixations')
colorbar
colormap('jet')

subplot(2,2,2)
imagesc(tm,wfq,RC_meanpower)
axis xy
title('Right-Center Fixations')
colorbar
colormap('jet')


subplot(2,2,3)
imagesc(tm,wfq,LL_meanpower)
axis xy
title('Left-Left Fixations')
colorbar
colormap('jet')


subplot(2,2,4)
imagesc(tm,wfq,RR_meanpower)
axis xy
title('Right-Right Fixations')
colorbar
colormap('jet')



figure
subplot(2,2,1)
imagesc(tm,wfq,LC_meanpower-RC_meanpower)
axis xy
title('Left/Center-Right/Center Fixations')
colorbar
colormap('jet')

subplot(2,2,2)
imagesc(tm,wfq,LL_meanpower-RR_meanpower)
axis xy
title('Left/Left-Right/Right Fixations')
colorbar
colormap('jet')

subplot(2,2,3)
imagesc(tm,wfq,LL_meanpower+LC_meanpower-RR_meanpower-RC_meanpower)
axis xy
title('Left-Right Fixations')
colorbar
colormap('jet')

subplot(2,2,4)
imagesc(tm,wfq,LC_meanpower+RC_meanpower-RR_meanpower-LL_meanpower)
axis xy
title('Center-Out Fixations')
colorbar
colormap('jet')


if monkey == 1
    subtitle(['Vivian Power Differences'])
else
    subtitle(['Tobii  Power Differences'])
end


figure
subplot(2,2,1)
imagesc(tm,wfq,LC_phasevar-RC_phasevar)
axis xy
title('Left/Center-Right/Center Fixations')
colorbar
colormap('jet')

subplot(2,2,2)
imagesc(tm,wfq,LL_phasevar-RR_phasevar)
axis xy
title('Left/Left-Right/Right Fixations')
colorbar
colormap('jet')

subplot(2,2,3)
imagesc(tm,wfq,LL_phasevar+LC_phasevar-RR_phasevar-RC_phasevar)
axis xy
title('Left-Right Fixations')
colorbar
colormap('jet')

subplot(2,2,4)
imagesc(tm,wfq,LC_phasevar+RC_phasevar-RR_phasevar-LL_phasevar)
axis xy
title('Center-Out Fixations')
colorbar
colormap('jet')


if monkey == 1
    subtitle(['Vivian Phase Differences'])
else
    subtitle(['Tobii  Phase Differences'])
end

figure
subplot(2,2,1)
imagesc(tm,wfq,1-LC_phasevar)
axis xy
title('Left-Center Fixations')
colorbar
colormap('jet')

subplot(2,2,2)
imagesc(tm,wfq,1-RC_phasevar)
axis xy
title('Right-Center Fixations')
colorbar
colormap('jet')


subplot(2,2,3)
imagesc(tm,wfq,1-LL_phasevar)
axis xy
title('Left-Left Fixations')
colorbar
colormap('jet')


subplot(2,2,4)
imagesc(tm,wfq,1-RR_phasevar)
axis xy
title('Right-Right Fixations')
colorbar
colormap('jet')

if monkey == 1
    subtitle(['Vivian Intersapos_CCadic Phase Consistency'])
else
    subtitle(['Tobii  Intersapos_CCadic Phase Consistency'])
end


figure
imagesc(tm,wfq,1-avg_phasevar)
axis xy
colorbar
colormap('jet')
if monkey == 1
    subtitle(['Vivian Avg Subsampled Intersapos_CCadic Phase Consistency'])
else
    subtitle(['Tobii   Avg Subsampled Intersapos_CCadic Phase Consistency'])
end

%%  Calculated Boostrapped Phase and Power Differences between Different Postion x Direction Combinations

Lc_ind = randperm(size(allLc,1));
Lc_ind = Lc_ind(1:max_samples);
this_allLc = allLc(Lc_ind,:);
[LC_trialpower,LC_trialphase,wfq,LC_meanpower,LC_meanphase,LC_powervar,LC_phasevar] = waveletanalysis(this_allLc);

RC_ind = randperm(size(allRc,1));
RC_ind = RC_ind(1:max_samples);
this_allRc = allRc(RC_ind,:);
[RC_trialpower,RC_trialphase,wfq,RC_meanpower,RC_meanphase,RC_powervar,RC_phasevar] = waveletanalysis(this_allRc);


LL_ind = randperm(size(allLL,1));
LL_ind = LL_ind(1:max_samples);
this_allLL = allLL(LL_ind,:);
[LL_trialpower,LL_trialphase,wfq,LL_meanpower,LL_meanphase,LL_powervar,LL_phasevar] = waveletanalysis(this_allLL);

RR_ind = randperm(size(allRR,1));
RR_ind = RR_ind(1:max_samples);
this_allRR = allRR(RR_ind,:);
[RR_trialpower,RR_trialphase,wfq,RR_meanpower,RR_meanphase,RR_powervar,RR_phasevar] = waveletanalysis(this_allRR);

[~,~,wfq,avg_meanpower,avg_meanphase,avg_powervar,avg_phasevar] = waveletanalysis(...
    [this_allLc(1:4:end,:); this_allLL(1:4:end,:); this_allRc(1:4:end,:); this_allRR(1:4:end,:)]);

%%

%---Left vs Rigth Saccades That were in center 1/3---%
lr_center_all_phase = cat(3,LC_trialphase(:,500:1500,:),RC_trialphase(:,500:1500,:));
lr_center_all_phase = exp(1i*lr_center_all_phase);
lr_center_all_power = cat(3,LC_trialpower(:,500:1500,:),RC_trialpower(:,500:1500,:));

 %---Left/Left vs Right/Right saccades that were in outer 1/3---%
ll_rr_all_phase = cat(3,LL_trialphase(:,500:1500,:),RR_trialphase(:,500:1500,:));
ll_rr_all_phase = exp(1i*ll_rr_all_phase);
ll_rr_all_power = cat(3,LL_trialpower(:,500:1500,:),RR_trialpower(:,500:1500,:));

%---All left vs all right saccades---%
lr_all_phase = cat(3,LL_trialphase(:,500:1500,:),LC_trialphase(:,500:1500,:),RR_trialphase(:,500:1500,:),RC_trialphase(:,500:1500,:));
lr_all_phase = exp(1i*lr_all_phase);
lr_all_power = cat(3,LL_trialpower(:,500:1500,:),LC_trialpower(:,500:1500,:),RR_trialpower(:,500:1500,:),RC_trialpower(:,500:1500,:));

%---All center vs all out---%
center_out_all_phase = cat(3,LC_trialphase(:,500:1500,:),RC_trialphase(:,500:1500,:),LL_trialphase(:,500:1500,:),RR_trialphase(:,500:1500,:));
center_out_all_phase = exp(1i*center_out_all_phase);
center_out_all_power = cat(3,LC_trialpower(:,500:1500,:),RC_trialpower(:,500:1500,:),LL_trialpower(:,500:1500,:),RR_trialpower(:,500:1500,:));

%---Categorical Indeces---%
all_id = [ones(max_samples,1); 2*ones(max_samples,1)];
all_id2 = [ones(2*max_samples,1); 2*ones(2*max_samples,1)];
%% preallocate space

shuffled_power_diff_center_lr = cell(1,numshuffs);
shuffled_phase_diff_center_lr = cell(1,numshuffs);
shuffled_phase_diff_ll_rr = cell(1,numshuffs);
shuffled_power_diff_ll_rr = cell(1,numshuffs);
shuffled_phase_diff_lr = cell(1,numshuffs);
shuffled_power_diff_lr = cell(1,numshuffs);
shuffled_phase_diff_center_out = cell(1,numshuffs);
shuffled_power_diff_center_out = cell(1,numshuffs);

num_frequencies = length(wfq);
num_time_pts = size(lr_center_all_phase,2);
for shuff = 1:numshuffs
    shuffled_power_diff_center_lr{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_phase_diff_center_lr{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_phase_diff_ll_rr{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_power_diff_ll_rr{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_phase_diff_lr{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_power_diff_lr{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_phase_diff_center_out{shuff} = NaN(num_frequencies,num_time_pts);
    shuffled_power_diff_center_out{shuff} = NaN(num_frequencies,num_time_pts);
end
%%

for shuff = 1:numshuffs
    disp(['Shuffled #' num2str(shuff)])
    rand_ind =  randperm(2*max_samples);
    rand_ind2 =  randperm(4*max_samples);

    %---Left vs Rigth Saccades That were in center 1/3---%
    tempsum=sum(lr_center_all_phase(:,:,all_id(rand_ind) == 1),3);
    MLR1=(abs(tempsum)./max_samples);
    tempsum=sum(lr_center_all_phase(:,:,all_id(rand_ind) == 2),3);
    MLR2=(abs(tempsum)./max_samples);
    shuffled_phase_diff_center_lr{shuff} = MLR1-MLR2;

    meanpower1=mean(lr_center_all_power(:,:,all_id(rand_ind) == 1),3);
    meanpower2=mean(lr_center_all_power(:,:,all_id(rand_ind) == 2),3);
    shuffled_power_diff_center_lr{shuff} = meanpower1-meanpower2;
    
    %---Left/Left vs Right/Right saccades that were in outer 1/3---%
    tempsum=sum(ll_rr_all_phase(:,:,all_id(rand_ind) == 1),3);
    MLR1=(abs(tempsum)./max_samples);
    tempsum=sum(ll_rr_all_phase(:,:,all_id(rand_ind) == 2),3);
    MLR2=(abs(tempsum)./max_samples);
    shuffled_phase_diff_ll_rr{shuff} = MLR1-MLR2;

    meanpower1=mean(ll_rr_all_power(:,:,all_id(rand_ind) == 1),3);
    meanpower2=mean(ll_rr_all_power(:,:,all_id(rand_ind) == 2),3);
    shuffled_power_diff_ll_rr{shuff} = meanpower1-meanpower2;
    
    %---All left vs all right saccades---%
    tempsum=sum(lr_all_phase(:,:,all_id(rand_ind) == 1),3);
    MLR1=(abs(tempsum)./max_samples);
    tempsum=sum(lr_all_phase(:,:,all_id(rand_ind) == 2),3);
    MLR2=(abs(tempsum)./max_samples);
    shuffled_phase_diff_lr{shuff} = MLR1-MLR2;

    meanpower1=mean(lr_all_power(:,:,all_id(rand_ind) == 1),3);
    meanpower2=mean(lr_all_power(:,:,all_id(rand_ind) == 2),3);
    shuffled_power_diff_lr{shuff} = meanpower1-meanpower2;
    
    
    %---All center vs all out---%
    tempsum=sum(center_out_all_phase(:,:,all_id(rand_ind) == 1),3);
    MLR1=(abs(tempsum)./max_samples);
    tempsum=sum(center_out_all_phase(:,:,all_id(rand_ind) == 2),3);
    MLR2=(abs(tempsum)./max_samples);
    shuffled_phase_diff_center_out{shuff} = MLR1-MLR2;

    meanpower1=mean(center_out_all_power(:,:,all_id(rand_ind) == 1),3);
    meanpower2=mean(center_out_all_power(:,:,all_id(rand_ind) == 2),3);
    shuffled_power_diff_center_out{shuff} = meanpower1-meanpower2;
    
end


%% Plot Phase Differences Left/Center vs Right/Center
observed_phase_diff_center_lr = LC_phasevar-RC_phasevar;
[all_mc_matrix_phase_lr,pos_big_block_matrix_phase_lr,neg_big_block_matrix_phase_lr,observed_prctile_phase_lr] = ...
    cluster_level_statistic2D(observed_phase_diff_center_lr(:,500:1500),shuffled_phase_diff_center_lr,min_cluster_area);
rectified_phase_observed_prctile =  observed_prctile_phase_lr;
rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50) = 100-rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm,wfq,observed_phase_diff_center_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Phase Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm,wfq,observed_prctile_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm,wfq,rectified_phase_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_phase_lr,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('Phase Differences Left/Center-Right/Center')

observed_power_diff_center_lr = LC_meanpower-RC_meanpower;
[all_mc_matrix_power_lr,pos_big_block_matrix_power_lr,neg_big_block_matrix_power_lr,observed_prctile_power_lr] = ...
    cluster_level_statistic2D(observed_power_diff_center_lr(:,500:1500),shuffled_power_diff_center_lr,min_cluster_area);
rectified_power_observed_prctile =  observed_prctile_power_lr;
rectified_power_observed_prctile(rectified_power_observed_prctile < 50) = 100-rectified_power_observed_prctile(rectified_power_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm,wfq,observed_power_diff_center_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed power Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm,wfq,observed_prctile_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm,wfq,rectified_power_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_power_lr,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('power Differences Left/Center-Right/Center')
%% Plot Phase Differences Left/Left vs Right/Right
observed_phase_diff_ll_rr = LL_phasevar-RR_phasevar;
[all_mc_matrix_phase_lr,pos_big_block_matrix_phase_lr,neg_big_block_matrix_phase_lr,observed_prctile_phase_lr] = ...
    cluster_level_statistic2D(observed_phase_diff_ll_rr(:,500:1500),shuffled_phase_diff_ll_rr,min_cluster_area);
rectified_phase_observed_prctile =  observed_prctile_phase_lr;
rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50) = 100-rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm,wfq,observed_phase_diff_ll_rr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Phase Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm,wfq,observed_prctile_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm,wfq,rectified_phase_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_phase_lr,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('Phase Differences Left/Left vs Right/Right')

observed_power_diff_ll_rr = LL_meanpower-RR_meanpower;
[all_mc_matrix_power_lr,pos_big_block_matrix_power_lr,neg_big_block_matrix_power_lr,observed_prctile_power_lr] = ...
    cluster_level_statistic2D(observed_power_diff_ll_rr(:,500:1500),shuffled_power_diff_ll_rr,min_cluster_area);
rectified_power_observed_prctile =  observed_prctile_power_lr;
rectified_power_observed_prctile(rectified_power_observed_prctile < 50) = 100-rectified_power_observed_prctile(rectified_power_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm,wfq,observed_power_diff_ll_rr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed power Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm,wfq,observed_prctile_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm,wfq,rectified_power_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_power_lr,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('power Differences Left/Left vs Right/Right')

%% Plot Phase Differences All left vs all right
observed_phase_diff_ll_rr = abs(sum(lr_all_phase(:,:,all_id == 1),3))-...
    abs(sum(lr_all_phase(:,:,all_id == 2),3));
observed_phase_diff_ll_rr = observed_phase_diff_ll_rr/max_samples;
[all_mc_matrix_phase_lr,pos_big_block_matrix_phase_lr,neg_big_block_matrix_phase_lr,observed_prctile_phase_lr] = ...
    cluster_level_statistic2D(observed_phase_diff_ll_rr,shuffled_phase_diff_ll_rr,min_cluster_area);
rectified_phase_observed_prctile =  observed_prctile_phase_lr;
rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50) = 100-rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm2,wfq,observed_phase_diff_ll_rr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Phase Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm2,wfq,observed_prctile_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm2,wfq,rectified_phase_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_phase_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_phase_lr,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('Phase Differences All left vs all right')

observed_power_diff_ll_rr = mean(lr_all_power(:,:,all_id == 1),3)-...
    mean(lr_all_power(:,:,all_id == 2),3);
[all_mc_matrix_power_lr,pos_big_block_matrix_power_lr,neg_big_block_matrix_power_lr,observed_prctile_power_lr] = ...
    cluster_level_statistic2D(observed_power_diff_ll_rr,shuffled_power_diff_ll_rr,min_cluster_area);
rectified_power_observed_prctile =  observed_prctile_power_lr;
rectified_power_observed_prctile(rectified_power_observed_prctile < 50) = 100-rectified_power_observed_prctile(rectified_power_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm2,wfq,observed_power_diff_ll_rr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed power Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm2,wfq,observed_prctile_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm2,wfq,rectified_power_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_power_lr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_power_lr,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('Power Differences All left vs all right')

%% Plot Phase Differences All center vs All out
observed_phase_diff_ll_rr = abs(sum(center_out_all_phase(:,:,all_id == 1),3))-...
    abs(sum(center_out_all_phase(:,:,all_id == 2),3));
observed_phase_diff_ll_rr = observed_phase_diff_ll_rr/max_samples;
[all_mc_matrix_phase_center_out,pos_big_block_matrix_phase_center_out,neg_big_block_matrix_phase_center_out,observed_prctile_phase_center_out] = ...
    cluster_level_statistic2D(observed_phase_diff_ll_rr,shuffled_phase_diff_ll_rr,min_cluster_area);
rectified_phase_observed_prctile =  observed_prctile_phase_center_out;
rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50) = 100-rectified_phase_observed_prctile(rectified_phase_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm2,wfq,observed_phase_diff_ll_rr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Phase Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm2,wfq,observed_prctile_phase_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm2,wfq,rectified_phase_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_phase_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_phase_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_phase_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_phase_center_out,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('Phase Differences All center vs All out')

observed_power_diff_ll_rr = mean(center_out_all_power(:,:,all_id == 1),3)-...
    mean(center_out_all_power(:,:,all_id == 2),3);
[all_mc_matrix_power_center_out,pos_big_block_matrix_power_center_out,neg_big_block_matrix_power_center_out,observed_prctile_power_center_out] = ...
    cluster_level_statistic2D(observed_power_diff_ll_rr,shuffled_power_diff_ll_rr,min_cluster_area);
rectified_power_observed_prctile =  observed_prctile_power_center_out;
rectified_power_observed_prctile(rectified_power_observed_prctile < 50) = 100-rectified_power_observed_prctile(rectified_power_observed_prctile < 50);

figure
subplot(2,3,1)
imagesc(tm2,wfq,observed_power_diff_ll_rr)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed power Difference')
xlim([-500 500])

subplot(2,3,2)
imagesc(tm2,wfq,observed_prctile_power_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Observed Percentile')
xlim([-500 500])


subplot(2,3,3)
imagesc(tm2,wfq,rectified_power_observed_prctile)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
colorbar
title('Rectified Percentile')
xlim([-500 500])


tm2 = -500:500;
subplot(2,3,5)
imagesc(tm2,wfq,pos_big_block_matrix_power_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Positive Clusters')

subplot(2,3,6)
imagesc(tm2,wfq,neg_big_block_matrix_power_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Large Negative Clusters')

subplot(2,3,4)
imagesc(tm2,wfq,all_mc_matrix_power_center_out)
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
axis xy
title('Maris Corrected Significant Clusters')

subplot(2,3,1)
[B,L] = bwboundaries(all_mc_matrix_power_center_out,'noholes');
hold on
for K = 1:length(B)
    boundary = B{K};
   plot(boundary(:,2)-twin/2, boundary(:,1), 'w', 'LineWidth', 2)
end
hold off

subtitle('Power Differences All center vs All out')
