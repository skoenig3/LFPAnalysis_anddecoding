%draft trying to understand LFPs aligned to events. Essentially ERP
%analysis.

clar

twin1 = 200;% how much time to take before event cross appears and how much to ignore during fixation
twin2 = 5000;%how much time to look at after stimulus onset for short window
twin4 = 5000; %for long window on image on
numshuffs = 10000; %number of shuffles to do for bootstrapping
imageX = 800;
imageY = 600;

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

all_LFPs = cell(1,85);
fixation_durations = cell(2,85);
saccade_amplitudes = cell(2,85);
all_recurrence_measures = cell(2,85);
image_pair_count = zeros(1,85);

set = 0;
for monkey = 1:2
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
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr','fixationstats');
        
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
        time_lock_LFP = cell(2,length(LFPchannels));%event aligned spike trains
        for chan = 1:length(LFPchannels)
            time_lock_LFP{1,chan} = NaN(96,twin1+twin2);%novel
            time_lock_LFP{2,chan} = NaN(96,twin1+twin2);%repeat
        end
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        set = set+1;
        %---Prealoccate Memory---%
        fixation_durations{1,set} = NaN(96,40);
        fixation_durations{2,set} = NaN(96,40);
        saccade_amplitudes{1,set} = NaN(96,40);
        saccade_amplitudes{2,set} = NaN(96,40);
        all_recurrence_measures{1,set} = NaN(96,5);
        all_recurrence_measures{2,set} = NaN(96,5);

        
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
                    imgon = cfg.trl(t).alltim(cfg.trl(t).allval == image_on_code)-trial_start; %when image turns on
                    imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == image_off_code)-trial_start; %when image turns off
                    
                    for chan = 1:length(LFPchannels)
                        LFPs = data(LFPchannels(chan)).values{t}; %spike trains for this trial
                        time_lock_LFP{nvr(img_index),chan}(which_img(img_index),:) = LFPs(imgon-twin1:imgon+twin2-1);
                    end
                    
                    %---get fixation/saccaqde information---%
                    fixationtimes = fixationstats{t}.fixationtimes; %fixtaion start and end times
                    fixations = fixationstats{t}.fixations; %fixation locations
                    saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                    xy = fixationstats{t}.XY;
                    
                    %remove fixations that started before image turned on
                    invalid= find(fixationtimes(1,:) < imgon);
                    fixationtimes(:,invalid) = [];
                    fixations(:,invalid) = [];
                    invalid= find(saccadetimes(1,:) < imgon);
                    saccadetimes(:,invalid) = [];
                    
                    %remove fixations started ending after image turned off
                    invalid= find(fixationtimes(2,:) > imgoff);
                    fixationtimes(:,invalid) = [];
                    fixations(:,invalid) = [];
                    invalid= find(saccadetimes(2,:) > imgoff);
                    saccadetimes(:,invalid) = [];
                                        
                    fixdurs = diff(fixationtimes)+1; %calculate fixation duration
                    fixation_durations{nvr(img_index),set}(which_img(img_index),1:size(fixationtimes,2)) = fixdurs;
                    
                    for f = 1:size(fixationtimes,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                        if f == 1
                            sacamp = sqrt(sum((fixations(:,1)-[400;300]).^2));
                            saccade_amplitudes{nvr(img_index),set}(which_img(img_index),f) = sacamp;
                        else
                            prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                            if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                                continue;
                            end
                            sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                            saccade_amplitudes{nvr(img_index),set}(which_img(img_index),f) = sacamp;
                        end
                    end
                    
                    
                    [recurrence_rate,~,corm,laminarity,~,forward_trace,~,reverse_trace,~] = ...
                        calculate_auto_recurrence(fixations,40,48);
                    
                    all_recurrence_measures{nvr(img_index),set}(which_img(img_index),1) = recurrence_rate;
                    all_recurrence_measures{nvr(img_index),set}(which_img(img_index),2) = corm;
                    all_recurrence_measures{nvr(img_index),set}(which_img(img_index),3) = laminarity;
                    all_recurrence_measures{nvr(img_index),set}(which_img(img_index),4) = forward_trace;
                    all_recurrence_measures{nvr(img_index),set}(which_img(img_index),5) = reverse_trace;
                end
            end
        end
        
        %for comparing novel to repeat only want the trials in which both
        %images were shown
        rmv = []; %images to remove
        for img = 1:96
            ind = find(which_images == img);
            if length(ind) ~= 2 %so either novel or repeat but not both
                rmv = [rmv img];
            end
        end
        
        %remove unpaired images
        fixation_durations{1,set}(rmv,:) = NaN;
        fixation_durations{2,set}(rmv,:) = NaN;
        saccade_amplitudes{1,set}(rmv,:) = NaN;
        saccade_amplitudes{2,set}(rmv,:) = NaN;
        all_recurrence_measures{1,set}(rmv,:) = NaN;
        all_recurrence_measures{2,set}(rmv,:) = NaN;
        
        for chan = 1:size(time_lock_LFP,2)
            time_lock_LFP{1,chan}(rmv,:) = NaN;
            time_lock_LFP{2,chan}(rmv,:) = NaN;
        end
        all_LFPs{set} = time_lock_LFP;
        
        image_pair_count(set) = sum(~isnan(fixation_durations{1,set}(:,1)));
     
    end
end

%% Average Fixation Durations by Ordinal Fixation #
all_fix_durs = []; %all fixaiton durations
nov_fix_durs = NaN(size(fixation_durations,2),20); %median by set novel fixation durations
rep_fix_durs = NaN(size(fixation_durations,2),20); %median by set repeat fixation durations
for set = 1:size(fixation_durations,2);
    if ~isempty(fixation_durations{1,set})
        nov_durs = fixation_durations{1,set}; %novel images
        rep_durs = fixation_durations{2,set}; %repeat images
        
        %remove fixations shorter than 100 msin duration since removed from analysis
        nov_durs(nov_durs < 100) = NaN;
        rep_durs(rep_durs < 100) = NaN;
        
        nov_fix_durs(set,:) = nanmedian(nov_durs(:,1:20)); %novel
        rep_fix_durs(set,:) = nanmedian(rep_durs(:,1:20)); %repeat
        
        all_fix_durs = [all_fix_durs; nov_durs; rep_durs]; %all fixation durations
    end
end
nov_fix_durs  = laundry(nov_fix_durs);
rep_fix_durs  = laundry(rep_fix_durs);


%---Stats Test---%
p_wilx = signrank(median(nov_fix_durs),median(rep_fix_durs)); %Wilcoxon signed rank test for zero median between novel and repeated images
%nonparametric, repeated measures (across multiple fixations), analysiss
p_vals = [];
for f = 1:size(nov_fix_durs,2)
    [~,p_vals(f)] = ttest(nov_fix_durs(:,f),rep_fix_durs(:,f)); %paired ttest, not sure if really valid but easy to interpret
end

%---Plot Results---%
figure
hold on
plot(nanmean(nov_fix_durs))
errorb(1:20,nanmean(nov_fix_durs),nanstd(nov_fix_durs)./sqrt(sum(~isnan(nov_fix_durs))),'color','b')
plot(nanmean(rep_fix_durs),'r')
errorb(1:20,nanmean(rep_fix_durs),nanstd(rep_fix_durs)./sqrt(sum(~isnan(rep_fix_durs))),'color','r')
for f = 1:size(nov_fix_durs,2)
    if p_vals(f) < 0.05/size(nov_fix_durs,2) %Bonferroni correction
        plot(f,215,'k*')
    end
end
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
xlim([0 21])
ylim([160 220])
axis square
legend('Novel','Repeat')
title(['PW and TO n_{sessions} = ' num2str(size(fixation_durations,2)) ', p_{Wilcoxon} = ' num2str(p_wilx,3)])
%% Average Image Onset aligned LFPs

all_nov_LFP = [];
all_rep_LFP = [];
for set = 1:size(all_LFPs,2)
    if ~isempty(all_LFPs{set})
        for chan = 1:size(all_LFPs{set},2)
            all_nov_LFP = [all_nov_LFP; nanmean(all_LFPs{set}{1,chan})];
            all_rep_LFP = [all_rep_LFP; nanmean(all_LFPs{set}{2,chan})];
        end
    end
end
%%
tm = -twin1:twin2-1;
figure
hold on
plot(tm,nanmean(all_nov_LFP),'blue')
plot(tm,nanmean(all_rep_LFP),'red')
plot([0 0],[-45,30],'k--')
hold off
xlabel('Time from Image Onset (ms)')
ylabel('LFP (uV)')
xlim([-twin1 twin2])
ylim([-45 30])
title('All Novel/Repeat Images')

%%
high_recog_nov = [];
low_recog_nov = [];
high_recog_rep = [];
low_recog_rep = [];


high_recog_nov_durs = [];
low_recog_nov_durs = [];
high_recog_rep_durs = [];
low_recog_rep_durs = [];


for set = 1:size(fixation_durations,2);
    if ~isempty(all_LFPs{set})
        nov_durs = fixation_durations{1,set}; %novel images
        rep_durs = fixation_durations{2,set}; %repeat images

        mean_fix_durs = nanmean(nov_durs(:,4:end)');
        rep_recurrence_measures = all_recurrence_measures{1,set}(:,:);


%         mean_fix_durs = nanmean(rep_durs(:,1:end)');
%         rep_recurrence_measures = all_recurrence_measures{2,set}(:,:);
%         rep_recurrence_measures(:,1) = [];
%         rep_recurrence_measures(:,3) = [];
        
        nandurs = find(isnan(mean_fix_durs));
        nanmeasures = find(isnan(sum(rep_recurrence_measures',1)));
        nanvals = unique([nandurs,nanmeasures]);
        
        
        group_vals = [mean_fix_durs' rep_recurrence_measures];
        group_vals(nanvals,:) = [];
    
        [U,~,~] = pca(group_vals ,1);        
        
        low_thresh = prctile(U,33);
        high_thresh = prctile(U,66);
        
        high_ind = find(U >= high_thresh);
        low_ind = find(U <= low_thresh);
        
        high_recog_nov_durs = [high_recog_nov_durs; nanmean(nov_durs(high_ind,:))];
        low_recog_nov_durs = [low_recog_nov_durs; nanmean(nov_durs(low_ind,:))];
        high_recog_rep_durs = [high_recog_rep_durs; nanmean(rep_durs(high_ind,:))];
        low_recog_rep_durs = [low_recog_rep_durs; nanmean(rep_durs(low_ind,:))];
        
        for chan = 1:size(all_LFPs{set},2)
            
            nov_LFPs = all_LFPs{set}{1,chan};
            nov_LFPs(nanvals,:) = [];
            rep_LPFs = all_LFPs{set}{2,chan};
            rep_LPFs(nanvals,:) = [];
            
            high_recog_nov = [high_recog_nov; nanmean(nov_LFPs(high_ind,1:2200))];
            high_recog_rep = [high_recog_rep; nanmean(rep_LPFs(high_ind,1:2200))];
            
            low_recog_nov =  [low_recog_nov; nanmean(nov_LFPs(low_ind,1:2200))];
            low_recog_rep  = [low_recog_rep; nanmean(rep_LPFs(low_ind,1:2200))];
        end
    end
end

%---Plot Results---%
figure
hold on
errorbar(nanmean(high_recog_nov_durs(:,1:20)),nanstd(high_recog_nov_durs(:,1:20))./sqrt(85))
errorbar(nanmean(low_recog_nov_durs(:,1:20)),nanstd(low_recog_nov_durs(:,1:20))./sqrt(85))
errorbar(nanmean(high_recog_rep_durs(:,1:20)),nanstd(high_recog_rep_durs(:,1:20))./sqrt(85))
errorbar(nanmean(low_recog_rep_durs(:,1:20)),nanstd(low_recog_rep_durs(:,1:20))./sqrt(85))
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('High Nov','Low Nov','High Rep','Low Rep')
xlim([0 21])
ylim([160 260])



numshuffs = 10000;
min_cluster_size = 30;

all_novrep = [high_recog_nov; high_recog_rep];
index = [ones(1,size(high_recog_nov,1)) 2*ones(1,size(high_recog_rep,1))];
all_curves = NaN(numshuffs,size(high_recog_rep,2));

for shuff = 1:numshuffs
    nind = randperm(length(index));
    ind = index(nind);
    all_curves(shuff,:) = nanmean(all_novrep(ind == 1,:))- nanmean(all_novrep(ind == 2,:));    
end
observed_diff = nanmean(high_recog_nov)-nanmean(high_recog_rep);
[~,high_sig_times] = cluster_level_statistic(observed_diff,all_curves,2,min_cluster_size); %multiple comparision corrected significant indeces

all_novrep = [low_recog_nov; low_recog_rep];
index = [ones(1,size(low_recog_nov,1)) 2*ones(1,size(low_recog_rep,1))];
all_curves = NaN(numshuffs,size(low_recog_rep,2));

parfor shuff = 1:numshuffs
    nind = randperm(length(index));
    ind = index(nind);
    all_curves(shuff,:) = nanmean(all_novrep(ind == 1,:))- nanmean(all_novrep(ind == 2,:));    
end
observed_diff = nanmean(low_recog_nov)-nanmean(low_recog_rep);
min_cluster_size = 30;
[~,low_sig_times] = cluster_level_statistic(observed_diff,all_curves,2,min_cluster_size); %multiple comparision corrected significant indeces

tm = -twin1:2000-1

figure
subplot(1,2,1)
hold on
plot(tm, nanmean(high_recog_nov),'b')
plot(tm, nanmean(high_recog_rep),'r')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
gaps = findgaps(find(high_sig_times));
if ~isempty(gaps)
    for g = 1:size(gaps,1)
        gp = gaps(g,:);
        gp(gp == 0) = [];
        h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
        uistack(h,'bottom')
        %set(h,'FaceAlpha',0.25)
    end
end
hold off
xlabel('Time from Image Onset (ms)')
ylabel('LFP (uV)')
title('High Recognition')
xlim([-twin1 2000])
ylim([-45 30])

subplot(1,2,2)
hold on
plot(tm, nanmean(low_recog_nov),'b')
plot(tm, nanmean(low_recog_rep),'r')
gaps = findgaps(find(low_sig_times));
if ~isempty(gaps)
    for g = 1:size(gaps,1)
        gp = gaps(g,:);
        gp(gp == 0) = [];
        h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
        uistack(h,'bottom')
        %set(h,'FaceAlpha',0.25)
    end
end

hold off
xlabel('Time from Image Onset (ms)')
ylabel('LFP (uV)')
title('Low Recognition')
xlim([-twin1 2000])
ylim([-45 30])

