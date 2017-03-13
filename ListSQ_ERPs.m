%draft trying to understand LFPs aligned to events. Essentially ERP
%analysis.

clar

twin1 = 200;% how much time to take before event cross appears and how much to ignore during fixation
twin2 = 1000;%how much time to look at after stimulus onset for short window
twin3 = 750;%minimum fixation on image cross hair duraiton
twin4 = 5000; %for long window on image on
numshuffs = 10000; %number of shuffles to do for bootstrapping
smval =60; %30 ms std, 1/2 width of gaussian smoothing filter
smval2 = 300;%150 ms std, 1/2 width of gaussian smoothing filter

% time_locked_firing since ITI period is defined by 2 events (15 & 16)
cross_on_code = 35;
fixation_code = 8;
image_on_code = 23;
image_off_code = 24;
reward_code = 3;
trial_start_code = 15; %ITI start
task = 'ListSQ';
min_blks = 2;

task = 'ListSQ';
%set(0,'DefaultFigureVisible','OFF');

all_image_on_short = [];
all_image_on_long = [];
all_cross_on = [];
all_fix_on_cross = [];
all_image_off = [];
all_which_monkey = [];
all_reward = [];
all_reward_end = [];

all_nov = [];
all_rep = [];

for monkey = 2:-1:1
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
        %%%---import task and chan data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [task_file,item_file,cnd_file,multichan,chan_names,chan_confidence,sorting_quality,waveform_count,lfp_quality,comments]...
            = get_task_data(session_data{sess},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');
        
        num_trials = length(cfg.trl); %number of trials
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
        
        %get important task specific information
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to events---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(task_file(1:8))
        
        %preallocate space and parallel structure of cfg
        time_lock_LFP = cell(length(LFPchannels),7);%event aligned spike trains
        for chan = 1:length(LFPchannels)
            time_lock_LFP{chan,1} = NaN(192,twin1+twin3);%cross hair on
            time_lock_LFP{chan,2} = NaN(192,twin1+twin2);%fix crosshair
            time_lock_LFP{chan,3} = NaN(192,twin1+twin2);%imgon
            time_lock_LFP{chan,4} = NaN(192,2*twin3);%imgoff
            time_lock_LFP{chan,5} = NaN(192,twin1+twin4);%img on long window
            time_lock_LFP{chan,6} = NaN(404,twin1+twin2);%reward start
            time_lock_LFP{chan,7} = NaN(404,twin2+twin2);%reward end
        end
        
        nvr = NaN(1,192);
        which_images = NaN(1,192);
        trial_index = 1;
        
        for t = 1:num_trials
            if any(cfg.trl(t).allval == image_on_code); %in which image was displayed or 1st item in sequence was displayed
                if (itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)) %then sequence trial
                    
                    %---image info---%
                    img_index = find(cfg.trl(t).cnd == img_cnd); %image index
                    if any(isnan(which_img(img_index)))
                        continue
                    end
                    nvr(img_index) = novel_vs_repeat(img_index); %whether image was novel or repeated
                    which_images(img_index) = which_img(img_index); %image number
                    
                    
                    %---Trial Event Codes/Times---%
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
                    trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after image off
                    crosson = cfg.trl(t).alltim(cfg.trl(t).allval == cross_on_code)-trial_start; %cross hair on
                    fix_cross = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_code)-trial_start; %fixation on cross hair according to cortex
                    fix_cross = fix_cross(1);%since fixation on image also counts as event 8
                    imgon = cfg.trl(t).alltim(cfg.trl(t).allval == image_on_code)-trial_start; %when image turns on
                    imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == image_off_code)-trial_start; %when image turns off
                    
                    if imgon-fix_cross > 750%cortex sometimes has minimal lag up to ~20 ms
                        imgon2 =  fix_cross+750;
                    else
                        imgon2 = imgon;
                    end
                    
                    
                    
                    for chan = 1:length(LFPchannels)
                        LFPs = data(LFPchannels(chan)).values{t}; %spike trains for this trial
                        
                        %---crosson aligned LFP---%
                        time_lock_LFP{chan,1}(img_index,:) = LFPs(crosson-twin1:crosson+twin3-1);
                        
                        %---fix cross aligned LFP---%
                        %variable duration event
                        time_lock_LFP{chan,2}(img_index,1:(imgon2-fix_cross+twin1)) = LFPs(fix_cross-twin1:imgon2-1);
                        
                        %---image on 1s aligned LFP---%
                        time_lock_LFP{chan,3}(img_index,:) = LFPs(imgon-twin1:imgon+twin2-1);
                        
                        %---image on 5 s aligned LFP---%
                        time_lock_LFP{chan,5}(img_index,:) = LFPs(imgon-twin1:imgon+twin4-1);
                        
                        %---image off aligned LFP---%
                        if t ~= length(cfg.trl)
                            %will need to look in next trial to get rest of data
                            %since usually only get ~200 ms of data
                            rest_time = twin3-(trial_end-imgoff);%twin3-time from image off to trial end
                            LFP2 = data(LFPchannels(chan)).values{t+1}(1:rest_time-2);
                            this_trial_LFP = LFPs(imgoff-twin3:end);
                            time_lock_LFP{chan,4}(img_index,:) = [this_trial_LFP LFP2];
                        end
                        
                        
                        
                    end
                elseif any(cfg.trl(t).allval == reward_code)%rewarded sequence trial
                    
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
                    trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after image off
                    rewards = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start;
                    reward_start = rewards(1);
                    reward_end = rewards(end)+115;%115 ms reward pulse duration
                    
                    LFPs = data(LFPchannels(chan)).values{t}; %spike trains for this trial
                    
                    %reward start aligned
                    time_lock_LFP{chan,6}(trial_index,:) = LFPs(reward_start-twin1:reward_start+twin2-1);
                    
                    %---image off aligned LFP---%
                    if t ~= length(cfg.trl)
                        %will need to look in next trial to get rest of data
                        %since usually only get ~200 ms of data
                        rest_time = twin2-(trial_end-reward_end);%twin3-time from image off to trial end
                        if rest_time > 0 
                        LFP2 = data(LFPchannels(chan)).values{t+1}(1:rest_time-2);
                        this_trial_LFP = LFPs(reward_end-twin2:end);
                        time_lock_LFP{chan,7}(trial_index,:) = [this_trial_LFP LFP2];
                        else
                            time_lock_LFP{chan,7}(trial_index,:) = LFPs(reward_end-twin2:reward_end+twin2-1);
                        end
                    end
                    trial_index = trial_index+1;
                end
            end
        end
        
        %for comparing novel to repeat only want the trials in which both
        %images were shown
        rmv = []; %images to remove
        for img = 1:96
            ind = find(which_images == img);
            if length(ind) ~= 2 %so either novel or repeat but not both
                rmv = [rmv ind];
            end
        end
        which_images(rmv) = NaN; %set to NaN to keep structure
        nvr(rmv) = NaN; %set to NaN to keep structure

        all_image_on_short = [all_image_on_short; nanmean(time_lock_LFP{chan,3})];
        all_image_on_long = [all_image_on_long; nanmean(time_lock_LFP{chan,5})];
        all_image_off = [all_image_off; nanmean(time_lock_LFP{chan,4})];
        all_cross_on = [all_cross_on; nanmean(time_lock_LFP{chan,1})];
        all_fix_on_cross = [all_fix_on_cross; nanmean(time_lock_LFP{chan,2})];
        all_which_monkey = [all_which_monkey monkey];
        
        all_reward = [all_reward; nanmean(time_lock_LFP{chan,6})];
        all_reward_end = [all_reward_end; nanmean(time_lock_LFP{chan,7})];
        %         all_nov = [all_nov; nanmean(time_lock_LFP{chan,3}(nvr == 1,:))];
        %         all_rep = [all_rep; nanmean(time_lock_LFP{chan,3}(nvr == 2,:))];
    end
end
%%
% 
numshuffs = 10000;
all_novrep = [all_nov; all_rep];
index = [ones(1,size(all_nov,1)) 2*ones(1,size(all_rep,1))];
all_curves = NaN(numshuffs,size(all_novrep,2));

for shuff = 1:numshuffs
    nind = randperm(length(index));
    ind = index(nind);
    all_curves(shuff,:) = nanmean(all_novrep(ind == 1,:))- nanmean(all_novrep(ind == 2,:));    
end
%%
observed_diff = nanmean(all_nov)-nanmean(all_rep);
smval = 30;
[~,list_sig_times] = cluster_level_statistic(observed_diff,all_curves,2,smval); %multiple comparision corrected significant indeces