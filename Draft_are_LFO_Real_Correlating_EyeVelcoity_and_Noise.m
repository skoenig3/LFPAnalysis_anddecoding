% clar

task = 'ListSQ';
image_on_twin = 500;
twin = 750;
fixwin = 5;%size of fixation window on each crosshair
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
imgon_code = 23;
imgoff_code = 24;
Fs = 1000;
min_blks = 2;


LFP_count = 0;

fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
buffer = 100;

session_avg_xc = [];
all_max_corr = [];
lag_max_corr = [];
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
        count = 0;
        all_xcs = zeros(1,2001);
        max_corr = NaN(1,5000);
        lag_max  = NaN(1,5000);
        
        all_trial_vel = [];
        all_trial_LFP = [];
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
            
            xy = fixationstats{trial}.XY; %x and y eye data
            
            
            % parse the input to remove the breaks in the eye data caused by looking
            % outside as indicated by the presence of NaNs
            svel = [];
            
            
            parsed_eyedat = preparse(xy);
            parsed_vel = cell(size(parsed_eyedat));
            vel_index = cell(size(parsed_eyedat));
            ind = 1;
            for p = 1:length(parsed_eyedat)
                if any(~isnan(parsed_eyedat{p}(1,:)))
                    
                    %raw velocity
                    xss = parsed_eyedat{p}(1,:);
                    yss = parsed_eyedat{p}(2,:);
                    
%                     %low pass filtered velocity 30 Hz cutoff
%                     x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
%                     y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
%                     xss = filtfilt(flt,1,x);
%                     yss = filtfilt(flt,1,y);
%                     xss = xss(101:end-101); %remove buffer after filtering
%                     yss = yss(101:end-101); %remove buffer after filtering
%                     x = x(101:end-101); %remove buffer after filtering
%                     y = y(101:end-101); %remove buffer after filtering
                    
                    svelx = diff(xss);
                    svely = diff(yss);
                    sv = sqrt(svelx.^2+svely.^2);
                    svel = [svel sv];
                    
                    parsed_vel{p} = sv;
                    vel_index{p} = ind:ind+length(sv)-1;
                    %                     ind = ind+length(sv);
                else
                    %                     svel = [svel NaN(1,size(parsed_eyedat{p},2))];
                    %                     parsed_vel{p} = NaN(1,size(parsed_eyedat{p},2));
                    %                     vel_index{p} = ind:ind+size(parsed_eyedat{p},2)-1;
                    ind = ind+size(parsed_eyedat{p},2);
                end
            end
            
            temp = cell(1,length(LFPchannels));
            for p = 1:length(parsed_eyedat)
                for chan = 1:length(LFPchannels);
                    trial_LFP = data(LFPchannels(chan)).values{trial};
                    %                     if any(~isnan(parsed_eyedat{p}(1,:)))
                    %                         if length(parsed_vel{p}) > 200
                    temp{chan}= [temp{chan} trial_LFP(vel_index{p})];
                end
            end
            
            temp = cell2mat(temp');
            all_trial_vel = [all_trial_vel svel];
            all_trial_LFP = [all_trial_LFP temp];
            %                 trial_LFP = filtfilt(flt,1,trial_LFP); %low pass filter
            %                 for p = 1:length(parsed_eyedat)
            %                     if any(~isnan(parsed_eyedat{p}(1,:)))
            %                         if length(parsed_vel{p}) > 2000
            %                             cr = xcorr(trial_LFP(vel_index{p}),parsed_vel{p},'coeff');
            %                             width = length(cr);
            %                             center = (width-1)/2;
            %                             all_xcs =  all_xcs+cr(center-1000:center+1000);
            %                             count = count+1;
            %
            %                             max_corr(count) = max(abs(cr));
            %                             lag_max  = center-find(abs(cr) == max(abs(cr)));
            %                         end
            %                     end
            %                 end
            
        end
        all_trial_vel = filtfilt(flt,1,all_trial_vel); %low pass filter
        for chan = 1:length(LFPchannels)
            all_trial_LFP(chan,:) = filtfilt(flt,1,all_trial_LFP(chan,:)); %low pass filter all_trial_LFP
            temp =  xcorr(all_trial_LFP(chan,200:end-200),all_trial_vel(200:end-200),1000,'coeff');
            session_avg_xc = [session_avg_xc; temp];
        end
%         session_avg_xc = [session_avg_xc; all_xcs/count];
%         all_max_corr = [all_max_corr max_corr(~isnan(max_corr))];
%         lag_max_corr = [lag_max_corr lag_max(~isnan(lag_max))];
    end
end
%%
all_max_corr