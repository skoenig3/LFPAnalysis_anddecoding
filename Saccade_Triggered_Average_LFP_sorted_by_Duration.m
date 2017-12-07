%modified from previous code to calculate intersaccadic coherence on
%October 11, 2016 by Seth Konig
%code also plots Saccade Triggered Average LFP

clar
task = 'ListSQ';

%---Important parametrs/values for ListSQ List---%
image_on_twin = 500;%how much time to ignore eye movements for to remove
%strong visual response though some may last longer
trial_start_code = 15;
trial_end_code = 20;
imgon_code = 23;
imgoff_code = 24;
Fs = 1000;
imageX = 800; %horizontal size of the image
imageY = 600; %horizontal size of the image
min_image_trials = 64; %only analyzes sessions with at least 2+ novel repeat blocks worth of image

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

%fixation durations by category minimum already set at 100 ms
fixdurthresh = [100:5:400 10000];
twin = 750; %time before saccade

all_session_avg_STAs = zeros(length(fixdurthresh),2*twin+1);
LFP_count = 0;
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = 'P:\\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = 'P:\\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        
        
        
        %---Load Important Task Information---%
        [task_file,item_file,cnd_file,~,~,~,~,~,lfp_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        if all(lfp_quality == 0) %all channels dead or something
            continue
        end
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],...
            'cfg','hdr','data','fixationstats');
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Get successful trials---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([task_file ': running saccade modulation anlaysis..'])
        
        %get important task specific information
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        
        %preallocate space and parallel structure of cfg
        num_trials = length(cfg.trl);
        image_trials = zeros(1,num_trials);
        for t = 1:num_trials %only take trials in which image was actually shown
            if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end);
                image_trials(t) = 1;
            end
        end
        
        %remove excess data associated with non image trials
        image_trials = find(image_trials);
        num_trials = length(image_trials);
        
        if num_trials < min_image_trials
            %shouldn't happen since only selected reasonably long recording sessions to analyze
            continue
        end
        
        %only take the eye data from successful image trials
        fixationstats = fixationstats(image_trials);
        cfg.trl = cfg.trl(image_trials);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Import LFP data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %also filters out line noise
        disp(['Importing List LFP Data for ' task_file(1:8) ])
        
        %code to find spike/LFP channels
        LFPchannels = find_desired_channels(cfg,'LFP');
        
        bad_channels = [];
        for channel = 1:4
            if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                if  lfp_quality(channel) == 0; %if it is bad
                    bad_channels = [bad_channels channel];
                end
            end
        end
        LFPchannels(bad_channels) = [];
        lfp_quality(bad_channels) = [];
        
        if isempty(LFPchannels)
            continue
        end
        this_session_avg_STAs = cell(1,length(LFPchannels));
        for channel = 1:length(LFPchannels)
            this_session_avg_STAs{channel} = zeros(length(fixdurthresh),2*twin+1);
        end
        dur_count = zeros(1,length(fixdurthresh));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Get Saccade and Fixation Times---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        list_saccade_times = NaN(num_trials,50); %start time of saccades
        saccade_aligned_fixation_durations = NaN(num_trials,50); %keep track of fixation duration following saccade because we may want to sort things out
        for t = 1:num_trials %will use all trials since may be useful to look at eye movements later for all trials
            fixationtimes = fixationstats{t}.fixationtimes;
            saccadetimes = fixationstats{t}.saccadetimes;
            xy = fixationstats{t}.XY;
            
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
            trial_end = cfg.trl(t).alltim(cfg.trl(t).allval == trial_end_code);
            imgon = cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code)-trial_start;
            imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == imgoff_code)-trial_start;
            max_trial_dur = length(data(LFPchannels(1)).values{image_trials(t)});
            
            
            % if monkey isn't paying attention data isn't probably
            % worth much plus have to cut off somewhere
            if imgoff-imgon > 1.5*imgdur-1
                imgoff = imgon+1.5*imgdur-1;
            end
            imgon = imgon+image_on_twin; %don't use data within 500 ms of image turning on in case of visual response
            
            %find fiations and saccades that did not occur during the image period;
            %should also take care of the 1st fixation on the crosshair
            
            %fixation started before image turned on
            invalid= find(fixationtimes(1,:) < imgon);
            fixationtimes(:,invalid) = [];
            
            %fixation started after the image turned off and/or firing rate could corrupted by image turning off
            invalid= find(fixationtimes(1,:) > imgoff);
            fixationtimes(:,invalid) = [];
            
            %saccade started before image turned on
            invalid= find(saccadetimes(1,:) < imgon);
            saccadetimes(:,invalid) = [];
            
            %saccade started after the image turned off and/or firing rate could corrupted by image turning off
            invalid= find(saccadetimes(1,:) > imgoff);
            saccadetimes(:,invalid) = [];
            
            for s = 1:size(saccadetimes,2);
                next_fix = find(fixationtimes(1,:) == saccadetimes(2,s)+1);%next fixation should start immediately after
                if isempty(next_fix) %trial ended or eye outside of image
                    continue %try next one
                end
                sacamp = sqrt(sum((xy(:,saccadetimes(2,s))-xy(:,saccadetimes(1,s))).^2)); %saccade amplitude
                next_fix_dur = fixationtimes(2,next_fix)-fixationtimes(1,next_fix)+1;%next fixation duration
                if sacamp >= min_sac_amp && next_fix_dur >= min_fix_dur %next fixation has to be long enough & saccade large enough
                    list_saccade_times(t,s) = saccadetimes(1,s);%start time of saccades
                    saccade_aligned_fixation_durations(t,s) =next_fix_dur;%duration of next fixation
                end
            end
            
            sacdur = saccade_aligned_fixation_durations(t,:);
            sacdur(isnan(sacdur)) = [];
            sact = list_saccade_times(t,:);
            sact(isnan(sact)) = [];
            for s = 1:length(sact)
                if sact(s)+twin < max_trial_dur
                    dur_index = find(sacdur(s) >= fixdurthresh(1:end-1) & (sacdur(s) <fixdurthresh(2:end)));
                    dur_count(dur_index) = dur_count(dur_index)+1;
                    for channel = 1:length(LFPchannels)
                        this_session_avg_STAs{channel}(dur_index,:) = this_session_avg_STAs{channel}(dur_index,:)+...
                            data(LFPchannels(channel)).values{image_trials(t)}(sact(s)-twin:sact(s)+twin);
                    end
                end
            end
        end
        
        for channel = 1:length(this_session_avg_STAs)
            for dur = 1:length(fixdurthresh)
                if dur_count(dur) ~= 0
                    this_session_avg_STAs{channel}(dur,:) = this_session_avg_STAs{channel}(dur,:)./dur_count(dur);
                end
            end
            
            all_session_avg_STAs = all_session_avg_STAs+this_session_avg_STAs{channel};
            LFP_count = LFP_count+1;
        end
    end
end

%%
figure
imagesc(-twin:twin,fixdurthresh(1:end-1),all_session_avg_STAs(1:end-1,:))
axis xy
colormap('viridis')
xlabel('Time from Saccade Onsset(ms)')
ylabel('Fixation Duration')
xlim([-250 750])