%modified from previous code to calculate intersaccadic coherence on
%October 11, 2016 by Seth Konig
%code also plots Saccade Triggered Average LFP

clar
task = 'ListSQ';

locations = [2:17];

save_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';

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

%Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%Fline = [59.9 60 60.1 119.9 120 120.1 179.9 180 180.1]; %more limited line noise frequencies to remove
Fline = [60 120 180]; %even more limited line noise frequencies to remove,averaging over so many trials so less important


buffer1 = 2048; %time before saccade
buffer2 = 2048; %time after saccade

%fixation durations by category minimum already set at 100 ms
fixdurthresh = [100 200;
                200 300;
                300 450;
                450 750];


for monkey = 1%:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Get Saccade and Fixation Times---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        list_saccade_times = NaN(num_trials,50); %start time of saccades
        saccade_aligned_fixation_durations = NaN(num_trials,50); %keep track of fixation duration following saccade because we may want to sort things out
        list_fixation_times = NaN(num_trials,50);%start time of fixation
        fixation_aligned_fixation_durations = NaN(num_trials,50); %keep track of fixation duration following saccade because we may want to sort things out
        list_startind = NaN(1,num_trials);%trial start index in all data
        
        for t = 1:num_trials %will use all trials since may be useful to look at eye movements later for all trials
            fixationtimes = fixationstats{t}.fixationtimes;
            saccadetimes = fixationstats{t}.saccadetimes;
            xy = fixationstats{t}.XY;
            
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
            list_startind(t) = trial_start;
            
            trial_end = cfg.trl(t).alltim(cfg.trl(t).allval == trial_end_code);
            imgon = cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code)-trial_start;
            imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == imgoff_code)-trial_start;
            
            
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
            for f = 1:size(fixationtimes,2);
                fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%next fixation duration
                if  fix_dur >= min_fix_dur %if fixation is long enough
                    prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%find fixaiton before in case want to know
                    if ~isempty(prior_sac);
                        sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))...
                            -xy(:,saccadetimes(1,prior_sac))).^2)); %prior saccade amplitude%saccade amplitude
                        if sacamp >= min_sac_amp
                            list_fixation_times(t,f) = fixationtimes(1,f); %fixation start time from trial start
                            fixation_aligned_fixation_durations(t,f) =  fix_dur;%next fixation duration
                        end
%                     else %was probably looking offscreen before and now on screen so sacamp is probably reasonable
%                         list_fixation_times(t,f) = fixationtimes(1,f); %fixation start time from trial start
%                         fixation_aligned_fixation_durations(t,f) =  fix_dur;%next fixation duration
                    end
                end
            end
        end
        
        %remove excess NaNs
        list_saccade_times = laundry(list_saccade_times); %start time of saccades
        saccade_aligned_fixation_durations = laundry(saccade_aligned_fixation_durations); %keep track of fixation duration following saccade because we may want to sort things out
        list_fixation_times = laundry(list_fixation_times);%start time of fixation
        fixation_aligned_fixation_durations = laundry(fixation_aligned_fixation_durations); %keep track of fixation duration following saccade because we may want to sort things out
        list_startind = laundry(list_startind);%trial start index in all data
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---convert Data into Field Trip Friendly Format---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
        
        %---For Saccades---%
        num_sac_list = sum(sum(~isnan(list_saccade_times)));
        list_saccade_aligned = NaN(num_sac_list,3);
        saccade_stretched_fixation_durations = NaN(1,num_sac_list);
        index = 1;
        for trl = 1:size(list_saccade_times,1)
            sacs = list_saccade_times(trl,:);
            fixdur = saccade_aligned_fixation_durations(trl,:);
            fixdur(isnan(sacs)) = [];
            sacs(isnan(sacs)) = [];
            for sacind = 1:length(sacs)
                saccade_stretched_fixation_durations(index) = fixdur(sacind);
                list_saccade_aligned(index,:) = [sacs(sacind)-1-buffer1+list_startind(trl) ...
                    sacs(sacind)-1+buffer2+list_startind(trl) -buffer1];
                index = index+1;
            end
        end
        %remove trials if data extends past end of file which occasionally
        %will happen if monkey ended on a image trial instead of a sequence trial
        too_late = find(list_saccade_aligned(:,2) > hdr.nSamples);
        list_saccade_aligned(too_late,:) = [];
        saccade_stretched_fixation_durations(too_late) = [];
        
        %---For Fixations---%
        num_fix_list = sum(sum(~isnan(list_fixation_times)));
        list_fixation_aligned = NaN(num_fix_list,3);
        fixation_stretched_fixation_durations = NaN(1,num_fix_list);
        index = 1;
        for trl = 1:size(list_fixation_times,1)
            fix = list_fixation_times(trl,:);
            fixdur = fixation_aligned_fixation_durations(trl,:);
            fixdur(isnan(fix)) = [];
            fix(isnan(fix)) = [];
            for sacind = 1:length(fix)
                fixation_stretched_fixation_durations(index) = fixdur(sacind);
                list_fixation_aligned(index,:) = [fix(sacind)-1-buffer1+list_startind(trl) ...
                    fix(sacind)-1+buffer2+list_startind(trl) -buffer1];
                index = index+1;
            end
        end
        %remove trials if data extends past end of file which occasionally
        %will happen if monkey ended on a image trial instead of a sequence trial
        too_late = find(list_fixation_aligned(:,2) > hdr.nSamples);
        list_fixation_aligned(too_late,:) = [];
        fixation_stretched_fixation_durations(too_late) = [];
        
        
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
        
        % read the data from file and preprocess them
        cfg.channel       = cfg.channel(LFPchannels)';
        cfg.dftfilter     = 'yes';
        cfg.dftfreq       = Fline;
        cfg.padding       = 1;
        cfg.continuous    = 'yes';
        
        %LFP data aligned saccades
        cfg_list_sac = cfg;
        cfg_list_sac.trl = list_saccade_aligned;
        list_saccade_aligneddata = ft_preprocessing(cfg_list_sac);
        
        %LFP data aligned fixations
        cfg_list_fix = cfg;
        cfg_list_fix.trl = list_fixation_aligned;
        list_fixation_aligneddata = ft_preprocessing(cfg_list_fix);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---remove any trials with NaNs, code can't process them well--%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %for saccades
        remove_trials = [];
        for t = 1:length(list_saccade_aligneddata.trial)
            if sum(sum(isnan(list_saccade_aligneddata.trial{t}))) >1
                remove_trials =[remove_trials t];
            end
        end
        list_saccade_aligneddata.time(remove_trials) = [];
        list_saccade_aligneddata.trial(remove_trials) = [];
        list_saccade_aligneddata.sampleinfo(remove_trials,:) = [];
        saccade_stretched_fixation_durations(remove_trials) = [];
        
        %for fixations
        remove_trials = [];
        for t = 1:length(list_fixation_aligneddata.trial)
            if sum(sum(isnan(list_fixation_aligneddata.trial{t}))) >1
                remove_trials =[remove_trials t];
            end
        end
        list_fixation_aligneddata.time(remove_trials) = [];
        list_fixation_aligneddata.trial(remove_trials) = [];
        list_fixation_aligneddata.sampleinfo(remove_trials,:) = [];
        fixation_stretched_fixation_durations(remove_trials) = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Save Data since takes a long time to preprocess---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %--get AP location of recording---%
        APs = session_data{session}.location(1)+chamber_zero;
        
        save([save_dir task_file(1:8) '_Eye_Aligned_LFP.mat'],'APs',...
            'list_saccade_aligneddata','saccade_stretched_fixation_durations',...
            'list_fixation_aligneddata','fixation_stretched_fixation_durations',...
            'buffer1','buffer2','lfp_quality');
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get Average Saccade Aligned Waveform Across Recordings Sessions---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clar
task = 'ListSQ';

locations = [2:17];

save_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';

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

%Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%Fline = [59.9 60 60.1 119.9 120 120.1 179.9 180 180.1]; %more limited line noise frequencies to remove
Fline = [60 120 180]; %even more limited line noise frequencies to remove,averaging over so many trials so less important

Fs = 1000; %sanmpling frequence
time_window = -0.75:.01:1.5; %what window to consider eye data over

buffer1 = 2048; %time before saccade
buffer2 = 2048; %time after saccade

%fixation durations by category minimum already set at 100 ms
fixdurthresh = [100 200;
                200 300;
                300 450;
                450 750];

time = (-buffer1:buffer2)/1000;
colors ={'red','magenta','green','blue','cyan'};

%row 1 saccade-aligned, row 2 fixation aligned
average_waveform = cell(2,length(fixdurthresh)); %average LFP waveform for various duration fixations
all_average_waveforms = cell(1,2);%average LFP across all fixations durations
all_which_monkey = []; %Tobii or Vivian
all_APs = [];%AP location of recording chamber

session_count = 0;
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
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
        disp(['Session # ' num2str(session) ': ' task_file(1:8)])
        
        if  exist([save_dir task_file(1:8) '_Eye_Aligned_LFP.mat'],'file')
            load([save_dir task_file(1:8) '_Eye_Aligned_LFP.mat']);
        else
            continue
        end
        session_count = session_count+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---LFPs For Different Fixation Durations---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for fdt = 1:size(fixdurthresh,1)
            fixind = find(saccade_stretched_fixation_durations >= fixdurthresh(fdt,1) & saccade_stretched_fixation_durations <= fixdurthresh(fdt,2));
            fixind2 = find(fixation_stretched_fixation_durations >= fixdurthresh(fdt,1) & fixation_stretched_fixation_durations <= fixdurthresh(fdt,2));
            if length(fixind) > 10 %need at least 10
                
                %---for saccades---%
                temp_aligned = list_saccade_aligneddata;
                temp_aligned.trial = temp_aligned.trial(fixind);
                
                all_sac = cell2mat(temp_aligned.trial');
                sac_electrode_avg = []; %average across trials for each electrode
                for el = 1:length(lfp_quality)
                    sac_electrode_avg = [sac_electrode_avg; nanmean(all_sac(el:length(lfp_quality):end,:))];
                end
                average_waveform{1,fdt} =[average_waveform{1,fdt}; sac_electrode_avg];

                %---for fixations---%
                temp_aligned = list_fixation_aligneddata;
                temp_aligned.trial = temp_aligned.trial(fixind);
                
                all_fix = cell2mat(temp_aligned.trial');
                fix_electrode_avg = []; %average across trials for each electrode
                for el = 1:length(lfp_quality)
                    fix_electrode_avg = [fix_electrode_avg; nanmean(all_fix(el:length(lfp_quality):end,:))];
                end
                average_waveform{2,fdt} =[average_waveform{2,fdt}; fix_electrode_avg];

            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---LFPs All Fixation Durations---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %---For Saccades---%
        all_sac = cell2mat(list_saccade_aligneddata.trial');
        sac_electrode_avg = [];
        for el = 1:length(lfp_quality)
            sac_electrode_avg = [sac_electrode_avg; nanmean(all_sac(el:length(lfp_quality):end,:))];
        end
        all_average_waveforms{1} = [all_average_waveforms{1}; sac_electrode_avg];

        %---For fixations---%
        all_fix = cell2mat(list_fixation_aligneddata.trial');
        fix_electrode_avg = [];
        for el = 1:length(lfp_quality)
            fix_electrode_avg = [fix_electrode_avg; nanmean(all_fix(el:length(lfp_quality):end,:))];
        end
        all_average_waveforms{2} = [all_average_waveforms{2}; fix_electrode_avg];
        
                
        %---session information---%
        all_which_monkey = [all_which_monkey monkey*ones(1,length(lfp_quality))]; %Tobii or Vivian
        APs =  session_data{session}.location(1)+chamber_zero;
        all_APs = [all_APs APs(1)*ones(1,length(lfp_quality))];%AP location of recording chamber
    end
end

save([save_dir 'ListSQ_List_eye_aligned_waveforms_across_sessions.mat'],...
    'all_average_waveforms','average_waveform','all_APs','all_which_monkey',...
    'buffer1','buffer2');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot average LFP waveforms---%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Saccade and Fixation Aligned Average Waveforms by Various Fixation Durations Different SubPlots---%
figure
subplot(2,3,[1 2])
hold on
dofill(time,all_average_waveforms{1},'blue',1,1);%even trials
dofill(time,all_average_waveforms{2},'red',1,1);%odd trials
plot([-0.25 0.5],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.5])
xlabel('Time From Eye Movement (ms)')
legend('Saccade','Fixation','Locaition','NorthWest')
ylabel('LFP (uV')
title('All Fixation Durations')

for fdt = 1:size(fixdurthresh,1)
    xmax = fixdurthresh(fdt,2)/1000;
    subplot(2,3,fdt+2)
    hold on
    dofill(time,average_waveform{1,fdt},'blue',1,[]);%even trials
    dofill(time,average_waveform{2,fdt},'red',1,[]);%odd trials
    plot([-0.25 xmax],[0 0],'k')
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-0.25 xmax])
    xlabel('Time From Eye Movement (ms)')
    ylabel('LFP (uV')
    title(['Fixation Durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2)) ' ms'])
end

subtitle(['Eye Movement Aligned LFP n_{electrodes} = ' num2str(size(all_average_waveforms{1},1))])

%---Saccade Aligned Average Waveforms by Various Fixation Durations Same Plot---%
figure
hold on
dofill(time,all_average_waveforms{1},'k',1,[]);%even trials
for fdt = 1:size(fixdurthresh,1)
    dofill(time,average_waveform{1,fdt},colors{fdt},1,[]);%even trials
end
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
legend('All','100-200','200-300','300-450','450-750')
ylabel('LFP (uV')
title(['All Fixation Durations n_{electrodes} = ' num2str(size(all_average_waveforms{1},1))])

%% Perform PCA on STA LFP
ind = find(time > 0 & time <= 0.4);
[U,S,V] = pca(all_average_waveforms{1}(:,ind));
T = kmeans(U,4);
figure
hold all
for i = 1:max(T)
    if sum(T == i) > 1
        plot(time,mean(all_average_waveforms{1}(T == i,:)))
    else
        plot(time,all_average_waveforms{1}(T == i,:));
    end
end
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])

figure
subplot(2,2,1)
plot(U(:,1),all_APs,'.k')
subplot(2,2,2)
plot(U(:,2),all_APs,'.k')
subplot(2,2,3)
plot(U(:,3),all_APs,'.k')
subplot(2,2,4)
plot(U(:,4),all_APs,'.k')


figure
hold on
dofill(time,all_average_waveforms{1}(all_which_monkey == 1,:),'r',1,[]);%Vivian
dofill(time,all_average_waveforms{1}(all_which_monkey == 2,:),'b',1,[]);%Tobii
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
legend('Vivian','Tobii')
ylabel('LFP (uV)')
title('Saccade Triggered Average LFP')

[U1,I1] = sort(U(:,2));
[U2,I2] = sort(U(:,2));
[U3,I3] = sort(U(:,3));
[U4,I4] = sort(U(:,4));

third = floor(size(U,1)/3);
figure
subplot(2,2,1)
hold on
dofill(time,all_average_waveforms{1}(I1(1:third),:),'r',1,[]);
dofill(time,all_average_waveforms{1}(I1(end-third:end),:),'b',1,[]);
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
ylabel('LFP (uV')
title('PCA1')

subplot(2,2,2)
hold on
dofill(time,all_average_waveforms{1}(I2(1:third),:),'r',1,[]);
dofill(time,all_average_waveforms{1}(I2(end-third:end),:),'b',1,[]);
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
ylabel('LFP (uV')
title('PCA2')

subplot(2,2,3)
hold on
dofill(time,all_average_waveforms{1}(I3(1:third),:),'r',1,[]);
dofill(time,all_average_waveforms{1}(I3(end-third:end),:),'b',1,[]);
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
ylabel('LFP (uV')
title('PCA3')

subplot(2,2,4)
hold on
dofill(time,all_average_waveforms{1}(I4(1:third),:),'r',1,[]);
dofill(time,all_average_waveforms{1}(I4(end-third:end),:),'b',1,[]);
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
ylabel('LFP (uV')
title('PCA4')
%%
figure
ap = all_APs(all_which_monkey == 1);
hist(ap,sort([locations locations+0.5]))
hold on
ap = all_APs(all_which_monkey == 2);
hist(ap,sort([locations locations+0.5]));
hold off
xlabel('AP Location (mm)')
ylabel('Count')
title('LFP electrode recording locations')
legend('Vivian','Tobii')
%%
[sorted,sort_index] = sort(all_APs);
avg = all_average_waveforms{1};
for a = 1:size(avg,1);
    avg(a,:) = avg(a,:)/max(abs(avg(a,:)));
end
avg = avg(sort_index,:);
sorted_monkey = all_which_monkey(sort_index);

figure
imagesc(time,1:size(avg,1),avg);
hold on
mn = find(sorted_monkey == 1);
plot(0,mn,'k*')
mn = find(sorted_monkey == 2);
plot(0,mn,'m*')
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('LFP ranked by AP')
xlim([-0.5 0.5])
colormap('jet')

subtitle('Saccade Triggered Average LFP by AP Location *black: vivian, *pink: tobii')
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
fdt = 5;
zcol = 1;
frequencies = low_freq{5,1,1}.freq;
times = low_freq{5,1,1}.time;

all_session_low_freq = cell(length(frequencies),length(times));
all_session_low_freq_coh= cell(length(frequencies),length(times));

for sess = 1:size(low_freq,2)
    for chan = 1:size(low_freq{fdt,sess,zcol}.powspctrm,1)
        for t = 1:length(times)
            for freq = 1:length(frequencies)
                all_session_low_freq{freq,t} =[all_session_low_freq{freq,t} low_freq{fdt,sess,zcol}.powspctrm(chan,freq,t)];
                all_session_low_freq_coh{freq,t} =[all_session_low_freq_coh{freq,t} low_freq_coh{fdt,sess,zcol}.cohspctrm(chan,freq,t)];
            end
        end
    end
end
%%
avg_power = cellfun(@mean,all_session_low_freq);
for a = 1:size(avg_power)
    avg_power(a,:) = avg_power(a,:)-mean(avg_power(a,:));
end
figure
imagesc(times,frequencies,avg_power)
axis xy
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
xlim([-0.5 0.5])
colormap('jet')
%%
avg_coh = cellfun(@mean,all_session_low_freq_coh);
figure
imagesc(times,frequencies,avg_coh)
axis xy
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
xlim([-0.5 0.5])
colormap('jet')

%%
fdt = 5;
zcol = 1;
frequencies = low_freq{5,1,1}.freq;
times = low_freq{5,1,1}.time;

all_session_low_freq = cell(2,length(frequencies),length(times));
all_session_low_freq_coh= cell(2,length(frequencies),length(times));

for sess = 1:size(low_freq,2)
    for chan = 1:size(low_freq{fdt,sess,zcol}.powspctrm,1)
        for t = 1:length(times)
            for freq = 1:length(frequencies)
                if all_which_monkey(sess)  == 1
                    all_session_low_freq{1,freq,t} =[all_session_low_freq{1,freq,t} low_freq{fdt,sess,zcol}.powspctrm(chan,freq,t)];
                    all_session_low_freq_coh{1,freq,t} =[all_session_low_freq_coh{1,freq,t} low_freq_coh{fdt,sess,zcol}.cohspctrm(chan,freq,t)];
                else
                    all_session_low_freq{2,freq,t} =[all_session_low_freq{2,freq,t} low_freq{fdt,sess,zcol}.powspctrm(chan,freq,t)];
                    all_session_low_freq_coh{2,freq,t} =[all_session_low_freq_coh{2,freq,t} low_freq_coh{fdt,sess,zcol}.cohspctrm(chan,freq,t)];
                end
            end
        end
    end
end
%%
vivian_avg_coh = cellfun(@mean,all_session_low_freq_coh(1,:,:));
vivian_avg_coh = reshape(vivian_avg_coh,size(vivian_avg_coh,2),size(vivian_avg_coh,3));
tobii_avg_coh = cellfun(@mean,all_session_low_freq_coh(2,:,:));
tobii_avg_coh = reshape(tobii_avg_coh,size(tobii_avg_coh,2),size(tobii_avg_coh,3));


figure
subplot(2,2,1)
imagesc(times,frequencies,avg_coh)
axis xy
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
xlim([-0.5 0.5])
colormap('jet')
title('Both')

subplot(2,2,3)
imagesc(times,frequencies,tobii_avg_coh-vivian_avg_coh)
axis xy
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
xlim([-0.5 0.5])
colormap('jet')
title('Tobii-Vivian')



subplot(2,2,2)
imagesc(times,frequencies,vivian_avg_coh)
axis xy
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
xlim([-0.5 0.5])
colormap('jet')
title('Vivian')

subplot(2,2,4)
imagesc(times,frequencies,tobii_avg_coh)
axis xy
xlabel('Time From Saccade Start (ms)')
ylabel('Frequency (Hz)')
xlim([-0.5 0.5])
colormap('jet')
title('Tobii')

subtitle('InterSaccdic Coherence')

%% get cohernce @ 1 time period
period = 1./frequencies;
line_averaged_coh = [];

line_averaged_coh = cell(2,length(frequencies));

for sess = 1:size(low_freq,2)
    for chan = 1:size(low_freq{fdt,sess,zcol}.powspctrm,1)
            for freq = 1:length(frequencies)
                [~,t] =min(abs(times-period(freq)));
                if all_which_monkey(sess)  == 1
                    line_averaged_coh{1,freq} =[line_averaged_coh{1,freq} low_freq_coh{fdt,sess,zcol}.cohspctrm(chan,freq,t)];
                else
                   line_averaged_coh{2,freq} =[line_averaged_coh{2,freq} low_freq_coh{fdt,sess,zcol}.cohspctrm(chan,freq,t)];
                end
            end
    end
end
%%
vivian_line_avg = cell2mat(line_averaged_coh(1,:)');
tobii_line_avg = cell2mat(line_averaged_coh(2,:)');

figure
hold on
dofill(frequencies,vivian_line_avg','r',1,[])
dofill(frequencies,tobii_line_avg','b',1,[])
hold off
ylabel('Coherence')
xlabel('Frequency(Hz)')
xlim([4 30])
title('Coherence @ 1 Time Period')

%%

vivian_avg_pow = cellfun(@mean,all_session_low_freq(1,:,:));
tobii_avg_pow = cellfun(@mean,all_session_low_freq(2,:,:));
avg_pow = vivian_avg_pow+tobii_avg_pow;
vivian_avg_pow = reshape(vivian_avg_pow,size(vivian_avg_pow,2),size(vivian_avg_pow,3));
tobii_avg_pow = reshape(tobii_avg_pow,size(tobii_avg_pow,2),size(tobii_avg_pow,3));
avg_pow = reshape(avg_pow,size(avg_pow,2),size(avg_pow,3));

%%
viv_avg = [];
tob_avg = [];
for f = 1:length(frequencies)
    viv_avg(f) = nanmean(vivian_avg_pow(f,:));
    tob_avg(f) = nanmean(tobii_avg_pow(f,:));
    vivian_avg_pow(f,:) = vivian_avg_pow(f,:)-viv_avg(f);
    tobii_avg_pow(f,:) = tobii_avg_pow(f,:)-tob_avg(f);
    avg_pow(f,:) = avg_pow(f,:)-nanmean(avg_pow(f,:));
end


figure
plot(frequencies,viv_avg./viv_avg(1),'r')
hold on
plot(frequencies,tob_avg./tob_avg(1),'b')
hold off
xlabel('Frequency (Hz)')
ylabel('Relative Power (uV^2)')
legend('Vivian','Tobii')

figure
subplot(2,2,1)
imagesc(times,frequencies,tobii_avg_pow-vivian_avg_pow)
xlabel('Time from Saccade Start (ms)')
ylabel('Frequency')
title('Tobii-Vivian')
colormap('jet')
axis xy

subplot(2,2,3)
imagesc(times,frequencies,vivian_avg_pow)
xlabel('Time from Saccade Start (ms)')
ylabel('Frequency')
title('Tobb+Vivian')
colormap('jet')
axis xy

subplot(2,2,2)
imagesc(times,frequencies,tobii_avg_pow)
xlabel('Time from Saccade Start (ms)')
ylabel('Frequency')
title('Tobii')
colormap('jet')
axis xy

subplot(2,2,4)
imagesc(times,frequencies,vivian_avg_pow)
xlabel('Time from Saccade Start (ms)')
ylabel('Frequency')
title('Vivian')
colormap('jet')
axis xy

subtitle('LFP Power around the time of Saccade')

%%
figure
subplot(2,2,1)
hold on
plot(frequencies,viv_avg,'r')
plot(frequencies,tob_avg','b')
hold off
ylabel('Power')
xlabel('Frequency(Hz)')
xlim([4 30])
legend('Vivian','Tobii')
title('Raw Power')

subplot(2,2,2)
hold on
plot(frequencies,viv_avg./viv_avg(1),'r')
plot(frequencies,tob_avg./tob_avg(1),'b')
hold off
ylabel('Power')
xlabel('Frequency(Hz)')
xlim([4 30])
legend('Vivian','Tobii')
title('Normalized to 4Hz Power')


F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);


Data1 = [frequencies' viv_avg'];
x01 = [1 1 0 0];
Data2 = [frequencies' tob_avg'];
x02 = [1 1 1 0];

t1 = Data1(:,1);
y1 = Data1(:,2);
[x1,resnorm,~,exitflag,output] = lsqcurvefit(F,x01,t,y);

t2 = Data1(:,1);
y2 = Data1(:,2);
[x2,resnorm,~,exitflag,output] = lsqcurvefit(F,x02,t,y);


subplot(2,2,3)
hold on
plot(frequencies,y1-F(x1,t),'r')
plot(frequencies,y2-F(x2,t),'b')
hold off
ylabel('Power')
xlabel('Frequency(Hz)')
xlim([4 30])
legend('Vivian','Tobii')
title('Subtracted 1/f fit')
%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot average LFP waveforms---%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Saccade and Fixation Aligned Average Waveforms by Various Fixation Durations Different SubPlots---%
figure
subplot(2,3,1)
hold on
[~,~,~,ysmth1,~]=dofill(time,all_average_waveforms{1},'blue',1,1);%even trials
[~,~,~,ysmth2,~]=dofill(time,all_average_waveforms{2},'red',1,1);%odd trials
plot([-0.25 0.5],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.5])
xlabel('Time From Eye Movement (ms)')
legend('Saccade','Fixation')
ylabel('LFP (uV')
title('All Fixation Durations')

for fdt = 1:size(fixdurthresh,1)
    xmax = fixdurthresh(fdt,2)/1000;
    subplot(2,3,fdt+1)
    hold on
    dofill(time,average_waveform{1,fdt},'blue',1,[]);%even trials
    dofill(time,average_waveform{2,fdt},'red',1,[]);%odd trials
    plot([-0.25 xmax],[0 0],'k')
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-0.25 xmax])
    xlabel('Time From Eye Movement (ms)')
    ylabel('LFP (uV')
    title(['Fixation Durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2)) ' ms'])
end

subtitle(['Eye Movement Aligned LFP n_{electrodes} = ' num2str(size(all_average_waveforms{1},1))])

%---Saccade Aligned Average Waveforms by Various Fixation Durations Same Plot---%
figure
hold on
dofill(time,all_average_waveforms{1},'k',1,[]);%even trials
for fdt = 1:size(fixdurthresh,1)
    dofill(time,average_waveform{1,fdt},colors{fdt},1,[]);%even trials
end
plot([-0.25 0.75],[0 0],'k')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.25 0.75])
xlabel('Time From Saccade Start (ms)')
legend('All','100-200','200-300','300-400','400-500','500-750')
ylabel('LFP (uV')
title('All Fixation Durations')
